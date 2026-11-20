/*
 * fgumi glue over htslib's CRAM internals.
 *
 * Compiled inside the vendored htslib tree (so it sees cram/cram.h and
 * hfile_internal.h). Everything here is either a thin accessor or a
 * re-statement of htslib's own single-threaded logic with the thread pool
 * removed, so that fgumi can run the expensive parts (slice decode, container
 * encode) on its own workers:
 *
 *   read:  fgumi_cram_next_slice()  -- serial: containers/slices off the stream
 *          fgumi_cram_decode_slice() -- any thread: decode + serialize as BAM
 *   write: cram_put_bam_seq() via sam_write1() hands full containers to the
 *          container executor (see vendor/container-executor.patch); then
 *          cram_encode_container() on any thread, cram_flush_encoded_container()
 *          serially.
 */
#include <config.h>

#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/hfile.h"
#include "htslib/hts_endian.h"
#include "htslib/hts_log.h"
#include "hfile_internal.h"
#include "cram/cram.h"
#include "cram/cram_decode.h"

/* ---------------------------------------------------------------------- */
/* hFILE backend over a Rust reader                                        */

typedef ssize_t (*fgumi_read_cb)(void *ctx, void *buf, size_t n);
typedef void (*fgumi_close_cb)(void *ctx);

typedef struct {
    hFILE base;
    void *ctx;
    fgumi_read_cb read_cb;
    fgumi_close_cb close_cb;
} hFILE_cb;

static ssize_t cb_read(hFILE *fpv, void *buf, size_t n) {
    hFILE_cb *fp = (hFILE_cb *)fpv;
    return fp->read_cb(fp->ctx, buf, n);
}
static ssize_t cb_write(hFILE *fpv, const void *buf, size_t n) {
    (void)fpv; (void)buf; (void)n;
    errno = EBADF;
    return -1;
}
static off_t cb_seek(hFILE *fpv, off_t off, int whence) {
    (void)fpv; (void)off; (void)whence;
    errno = ESPIPE;
    return -1;
}
static int cb_close(hFILE *fpv) {
    hFILE_cb *fp = (hFILE_cb *)fpv;
    if (fp->close_cb) fp->close_cb(fp->ctx);
    return 0;
}
static const struct hFILE_backend cb_backend = {
    cb_read, cb_write, cb_seek, NULL, cb_close
};

/*
 * Open a CRAM reader over callbacks. Returns NULL on failure (ctx is closed).
 *
 * *cram_version is set to the CRAM file definition's version (major * 256 +
 * minor) when the stream starts with one, else 0. A CRAM older than 2.1 is
 * refused (NULL) before htslib parses anything: 1.x and 2.0 have no EOF
 * container, so a truncated file reads exactly like a complete one (2.0) or
 * every complete file looks truncated (1.x).
 */
htsFile *fgumi_hts_open_callback_reader(void *ctx, fgumi_read_cb rd,
                                        fgumi_close_cb cl, const char *name,
                                        int *cram_version) {
    *cram_version = 0;
    hFILE_cb *hf = (hFILE_cb *)hfile_init(sizeof(hFILE_cb), "r", 1 << 20);
    if (!hf) { if (cl) cl(ctx); return NULL; }
    hf->ctx = ctx;
    hf->read_cb = rd;
    hf->close_cb = cl;
    hf->base.backend = &cb_backend;
    unsigned char def[6];
    if (hpeek(&hf->base, def, sizeof def) == (ssize_t)sizeof def
        && memcmp(def, "CRAM", 4) == 0) {
        *cram_version = def[4] << 8 | def[5];
        if (*cram_version < (2 << 8 | 1)) {
            hclose_abruptly(&hf->base);
            return NULL;
        }
    }
    htsFile *fp = hts_hopen(&hf->base, name, "r");
    if (!fp) { hclose_abruptly(&hf->base); return NULL; }
    return fp;
}

/* The CRAM handle of `fp`, or NULL when `fp` is not CRAM. */
cram_fd *fgumi_hts_cram_fd(htsFile *fp) {
    return fp && fp->format.format == cram ? fp->fp.cram : NULL;
}

/*
 * Keep whole reference sequences shared rather than a per-call partial window,
 * which concurrent slice decodes / container encodes would clobber. htslib sets
 * this itself whenever it is given a thread pool (CRAM_OPT_THREAD_POOL); we run
 * the same concurrency on our own threads, so we must set it too.
 */
void fgumi_cram_set_shared_ref(cram_fd *fd) {
    fd->shared_ref = 1;
}

/* ---------------------------------------------------------------------- */
/* Read: serial slice iteration                                            */

/*
 * Advance to the next slice. `cur` is the container the previous slice came
 * from (NULL at start); it is NOT freed here -- the caller owns container
 * lifetimes, since earlier slices of `cur` may still be decoding elsewhere.
 * On return 1, *c_out is the slice's container (either `cur` or a freshly read
 * one, which the caller now owns) and *s_out the slice (caller owns).
 *
 * Mirrors cram_next_slice() with ranges and the thread pool removed.
 * Returns 1 = slice, 0 = clean EOF, -2 = the stream ends without its EOF
 * container (truncated), -3 = stray bytes after the EOF container, -1 = any
 * other error (a read failure, or a corrupt or cut-short container/slice).
 */
/* EOF is a complete, empty compression-header block, not just its sentinel
 * container header. Validate the payload and checksum before probing physical
 * EOF; in particular never turn an I/O error into a successful end. */
static int fgumi_cram_read_eof(cram_fd *fd, cram_container *c) {
    static const unsigned char empty_maps[] = {1, 0, 1, 0, 1, 0};
    int expected_length = CRAM_MAJOR_VERS(fd->version) >= 3 ? 15 : 11;
    if (c->length != expected_length || c->ref_seq_span != 0 ||
        c->record_counter != 0 || c->num_bases != 0 ||
        c->num_blocks != 1 || c->num_landmarks != 0)
        return -1;
    cram_block *b = cram_read_block(fd);
    if (!b) return -1;
    int valid = b->method == RAW && b->content_type == COMPRESSION_HEADER &&
        b->content_id == 0 && b->comp_size == sizeof empty_maps &&
        b->uncomp_size == sizeof empty_maps &&
        cram_uncompress_block(b) == 0 &&
        memcmp(b->data, empty_maps, sizeof empty_maps) == 0;
    cram_free_block(b);
    if (!valid) return -1;
    unsigned char byte;
    ssize_t n = hpeek(fd->fp, &byte, 1);
    if (n < 0 || herrno(fd->fp)) return -1;
    if (n != 0) {
        hts_log_error("fgumi: CRAM input has stray bytes after its EOF container");
        return -3;
    }
    fd->eof = 1;
    return 0;
}

int fgumi_cram_next_slice(cram_fd *fd, cram_container *cur,
                          cram_container **c_out, cram_slice **s_out) {
    cram_container *c = cur;
    *c_out = NULL;
    *s_out = NULL;

    for (;;) {
        if (!c || c->curr_slice_mt == c->max_slice) {
            for (;;) {
                c = cram_read_container(fd);
                if (!c) {
                    if (herrno(fd->fp)) return -1;
                    if (fd->eof == 2) {
                        hts_log_error("fgumi: CRAM input ends without its EOF "
                                      "container; the file is truncated");
                        return -2;
                    }
                    return -1;
                }
                if (fd->empty_container) {
                    int ret = fgumi_cram_read_eof(fd, c);
                    cram_free_container(c);
                    return ret;
                }
                c->curr_slice_mt = c->curr_slice;
                if (c->length != 0)
                    break;
                cram_free_container(c);
            }
            if (!(c->comp_hdr_block = cram_read_block(fd))
                || c->comp_hdr_block->content_type != COMPRESSION_HEADER) {
                cram_free_container(c);
                return -1;
            }
            c->comp_hdr = cram_decode_compression_header(fd, c->comp_hdr_block);
            if (!c->comp_hdr) {
                cram_free_container(c);
                return -1;
            }
            if (!c->comp_hdr->AP_delta &&
                sam_hrecs_sort_order(fd->header->hrecs) != ORDER_COORD) {
                pthread_mutex_lock(&fd->ref_lock);
                fd->unsorted = 1;
                pthread_mutex_unlock(&fd->ref_lock);
            }
            if (c->num_records == 0) {
                cram_free_container(c);
                c = NULL;
                continue;
            }
            *c_out = c;  /* newly read: caller takes ownership */
        } else {
            *c_out = c;
        }

        cram_slice *s = cram_read_slice(fd);
        if (!s) {
            if (c != cur) {  /* newly read: still ours */
                cram_free_container(c);
                *c_out = NULL;
            }
            return -1;
        }
        s->slice_num = ++c->curr_slice_mt;
        s->curr_rec = 0;
        s->max_rec = s->hdr->num_records;
        s->last_apos = s->hdr->ref_seq_start;
        *s_out = s;
        return 1;
    }
}

/* ---------------------------------------------------------------------- */
/* Read: decode one slice straight to framed BAM records                   */

static int grow(uint8_t **buf, size_t *cap, size_t need) {
    if (need <= *cap) return 0;
    size_t n = *cap ? *cap : 1 << 16;
    while (n < need) n *= 2;
    uint8_t *p = realloc(*buf, n);
    if (!p) return -1;
    *buf = p;
    *cap = n;
    return 0;
}

static inline void put_u32(uint8_t *p, uint32_t v) { memcpy(p, &v, 4); }

/* Append `b` as one BAM record (block_size prefix included), byte for byte as
 * bam_write1() would emit it, including its long-CIGAR form: more than 65535
 * operations are written as a `<qlen>S<rlen>N` placeholder plus a CG:B,I tag
 * holding the real CIGAR (SAM spec 4.2.2). Little-endian hosts only (checked
 * at build time). */
static int append_bam(const bam1_t *b, uint8_t **buf, size_t *len, size_t *cap) {
    const bam1_core_t *c = &b->core;
    uint32_t l_qname = c->l_qname - c->l_extranul;
    if (l_qname > 255) {
        hts_log_error("fgumi: QNAME %s is longer than 254 characters", bam_get_qname(b));
        return -1;
    }
    if (c->pos > INT32_MAX || c->mpos > INT32_MAX ||
        c->isize > INT32_MAX || c->isize < INT32_MIN) {
        hts_log_error("fgumi: positional data of %s is too large for BAM", bam_get_qname(b));
        return -1;
    }
    int long_cigar = c->n_cigar > 0xffff;
    hts_pos_t cigreflen = 0;
    if (long_cigar) {
        cigreflen = bam_cigar2rlen(c->n_cigar, bam_get_cigar(b));
        if (cigreflen >= (1 << 28)) {
            hts_log_error("fgumi: record %s with %u CIGAR ops spans too much reference for BAM",
                          bam_get_qname(b), c->n_cigar);
            return -1;
        }
    }
    uint32_t block_len = 32 + b->l_data - c->l_extranul + (long_cigar ? 16 : 0);
    if (grow(buf, cap, *len + 4 + block_len) < 0) return -1;
    uint8_t *x = *buf + *len;
    put_u32(x, block_len);
    put_u32(x + 4, (uint32_t)c->tid);
    put_u32(x + 8, (uint32_t)(int32_t)c->pos);
    put_u32(x + 12, (uint32_t)c->bin << 16 | (uint32_t)c->qual << 8 | l_qname);
    put_u32(x + 16, (uint32_t)c->flag << 16 | (long_cigar ? 2u : (uint32_t)c->n_cigar));
    put_u32(x + 20, (uint32_t)c->l_qseq);
    put_u32(x + 24, (uint32_t)c->mtid);
    put_u32(x + 28, (uint32_t)(int32_t)c->mpos);
    put_u32(x + 32, (uint32_t)(int32_t)c->isize);
    uint8_t *p = x + 36;
    memcpy(p, b->data, l_qname);
    p += l_qname;
    if (!long_cigar) {
        memcpy(p, b->data + c->l_qname, b->l_data - c->l_qname);
    } else {
        uint32_t cigar_st = (uint8_t *)bam_get_cigar(b) - b->data;
        uint32_t cigar_en = cigar_st + c->n_cigar * 4;
        put_u32(p, (uint32_t)c->l_qseq << 4 | BAM_CSOFT_CLIP);
        put_u32(p + 4, (uint32_t)cigreflen << 4 | BAM_CREF_SKIP);
        p += 8;
        memcpy(p, b->data + cigar_en, b->l_data - cigar_en);
        p += b->l_data - cigar_en;
        memcpy(p, "CGBI", 4);
        put_u32(p + 4, c->n_cigar);
        memcpy(p + 8, b->data + cigar_st, (size_t)c->n_cigar * 4);
    }
    *len += 4 + block_len;
    return 0;
}

/* Upper bound on a slice's presized decode buffer (see fgumi_cram_decode_slice). */
#define FGUMI_PRESIZE_MAX ((size_t)1 << 28)

/*
 * Decode slice `s` of container `c` and append its records, BAM-encoded, to
 * a malloc'd buffer (*out, *out_len; free with free()). Safe to call
 * concurrently for different slices (htslib's own thread pool does the same).
 * Returns the number of records, or -1.
 */
int64_t fgumi_cram_decode_slice(cram_fd *fd, cram_container *c, cram_slice *s,
                                uint8_t **out, size_t *out_len) {
    *out = NULL;
    *out_len = 0;
    if (cram_decode_slice(fd, c, s, fd->header) != 0)
        return -1;

    int64_t n = s->hdr->num_records;
    /* Presize the output from the decoded records' sizes, so a typical slice
     * needs no realloc: each record's fixed fields, name, CIGAR, sequence,
     * qualities and stored tags, plus slack for the tags decoding adds
     * (RG, and MD/NM on mapped reads) and the long-CIGAR form. Only a hint:
     * append_bam() still grows the buffer as needed, so the estimate is
     * capped rather than trusted. */
    size_t cap = 0;
    for (int64_t i = 0; i < n && cap < FGUMI_PRESIZE_MAX; i++) {
        const cram_record *r = &s->crecs[i];
        size_t name = r->name_len > 0 ? (size_t)r->name_len : 0;
        size_t ncig = r->ncigar > 0 ? (size_t)r->ncigar : 0;
        size_t rlen = r->len > 0 ? (size_t)r->len : 0;
        cap += 36 + name + 1 + 4 * ncig + (rlen + 1) / 2 + rlen + r->aux_size + 64;
    }
    if (cap > FGUMI_PRESIZE_MAX) cap = FGUMI_PRESIZE_MAX;
    uint8_t *buf = cap ? malloc(cap) : NULL;
    if (cap && !buf) return -1;
    size_t len = 0;

    bam1_t *b = bam_init1();
    if (!b) { free(buf); return -1; }
    for (int64_t i = 0; i < n; i++) {
        if (cram_to_bam(fd->header, fd, s, &s->crecs[i], (int)i, &b) < 0
            || append_bam(b, &buf, &len, &cap) < 0) {
            bam_destroy1(b);
            free(buf);
            return -1;
        }
    }
    bam_destroy1(b);
    *out = buf;
    *out_len = len;
    return n;
}

/* Compressed bytes held by a read (undecoded) slice: its header and data blocks. */
size_t fgumi_cram_slice_bytes(const cram_slice *s) {
    size_t n = s->hdr_block ? (size_t)s->hdr_block->comp_size : 0;
    for (int j = 0; s->hdr && j < s->hdr->num_blocks; j++)
        if (s->block && s->block[j])
            n += (size_t)s->block[j]->comp_size;
    return n;
}

/* Admission estimate, not an adversarial allocation limit: account for the
 * whole container's bases (slice base counts are not stored), decoded blocks,
 * record structs/names and repeated tag-dictionary labels. The multiplier
 * covers native scratch, growing CIGAR/aux buffers and native + Rust BAM copies.
 * Reference caches are shared file state and are outside this per-slice charge. */
int fgumi_cram_slice_decoded_size(const cram_container *c, const cram_slice *s,
                                 uint64_t *out) {
    if (!s->hdr || c->num_bases < 0 || s->hdr->num_records < 0 ||
        s->hdr->num_blocks < 0 || s->hdr->ref_seq_span < 0)
        return -1;
    uint64_t n = (uint64_t)c->num_bases;
    uint64_t records = (uint64_t)s->hdr->num_records;
    uint64_t per_record = sizeof(cram_record) + 512;
    if (c->comp_hdr && c->comp_hdr->TD_blk) {
        uint64_t dictionary_bytes = (uint64_t)BLOCK_SIZE(c->comp_hdr->TD_blk);
        if (dictionary_bytes > UINT64_MAX - per_record) return -1;
        per_record += dictionary_bytes;
    }
    if (records && per_record > (UINT64_MAX - n) / records) return -1;
    n += records * per_record;
    uint64_t base_and_block_bytes = (uint64_t)c->num_bases;
    for (int j = 0; j < s->hdr->num_blocks; j++) {
        const cram_block *b = s->block ? s->block[j] : NULL;
        if (!b || b->uncomp_size < 0 || (uint64_t)b->uncomp_size > UINT64_MAX - n)
            return -1;
        n += (uint64_t)b->uncomp_size;
        if ((uint64_t)b->uncomp_size > UINT64_MAX - base_and_block_bytes) return -1;
        base_and_block_bytes += (uint64_t)b->uncomp_size;
    }
    /* A worst-case span * records MD allowance serializes normal mapped slices.
     * Reserve four times the bases/blocks instead. Repeated exceptional long
     * deletions may exceed this heuristic; this is not a hard allocation cap. */
    if (base_and_block_bytes > (UINT64_MAX - n) / 4) return -1;
    n += 4 * base_and_block_bytes;
    if (n > (UINT64_MAX - (1 << 20)) / 16) return -1;
    *out = n * 16 + (1 << 20);
    return 0;
}

void fgumi_free(void *p) { free(p); }

/* ---------------------------------------------------------------------- */
/* Write: raw BAM record -> bam1_t                                          */

/*
 * Move a long CIGAR stored in a CG:B,I tag back into the CIGAR field, as
 * bam_read1() does. A copy of htslib's static bam_tag2cigar() (sam.c), with
 * its buffer growth inlined. Returns 1 if moved, 0 if not applicable, -1 on
 * error.
 */
static int fgumi_tag2cigar(bam1_t *b) {
    bam1_core_t *c = &b->core;
    uint32_t test_CG = BAM_CSOFT_CLIP | (c->l_qseq << BAM_CIGAR_SHIFT);
    if (c->n_cigar == 0 || test_CG != *bam_get_cigar(b))
        return 0;
    if (c->tid < 0 || c->pos < 0)
        return 0;
    uint8_t *CG = bam_aux_get(b, "CG");
    if (!CG)
        return errno == ENOENT ? 0 : -1;
    if (CG[0] != 'B' || !(CG[1] == 'I' || CG[1] == 'i'))
        return 0;
    uint32_t CG_len = le_to_u32(CG + 2);
    if (CG_len < c->n_cigar || CG_len >= 1U << 29)
        return 0;

    uint32_t ori_len = b->l_data;
    uint32_t cigar_st = (uint8_t *)bam_get_cigar(b) - b->data;
    uint32_t fake_bytes = c->n_cigar * 4;
    uint32_t n_cigar4 = CG_len * 4;
    size_t need = (size_t)ori_len + n_cigar4 - fake_bytes;
    if (need > (size_t)b->m_data) {
        ptrdiff_t cg_off = CG - b->data;
        uint8_t *d = realloc(b->data, need);
        if (!d) return -1;
        b->data = d;
        b->m_data = need;
        CG = b->data + cg_off;
    }
    uint32_t CG_st = CG - b->data - 2;
    uint32_t CG_en = CG_st + 8 + n_cigar4;
    c->n_cigar = CG_len;
    b->l_data = b->l_data - fake_bytes + n_cigar4;
    memmove(b->data + cigar_st + n_cigar4, b->data + cigar_st + fake_bytes,
            ori_len - (cigar_st + fake_bytes));
    memcpy(b->data + cigar_st, b->data + (n_cigar4 - fake_bytes) + CG_st + 8, n_cigar4);
    if (ori_len > CG_en)
        memmove(b->data + CG_st + n_cigar4 - fake_bytes,
                b->data + CG_en + n_cigar4 - fake_bytes, ori_len - CG_en);
    b->l_data -= n_cigar4 + 8;
    return 1;
}

/*
 * Fill `b` from one raw BAM record (WITHOUT its 4-byte block_size prefix),
 * as bam_read1() would. Returns 0 or a negative error.
 */
int fgumi_bam_from_raw(bam1_t *b, const uint8_t *x, size_t len) {
    bam1_core_t *c = &b->core;
    uint32_t v;
    if (len < 32) return -4;
    memcpy(&v, x, 4);       c->tid = (int32_t)v;
    memcpy(&v, x + 4, 4);   c->pos = (int32_t)v;
    memcpy(&v, x + 8, 4);
    c->bin = v >> 16; c->qual = v >> 8 & 0xff; c->l_qname = v & 0xff;
    c->l_extranul = (c->l_qname % 4 != 0) ? (4 - c->l_qname % 4) : 0;
    memcpy(&v, x + 12, 4);  c->flag = v >> 16; c->n_cigar = v & 0xffff;
    memcpy(&v, x + 16, 4);  c->l_qseq = (int32_t)v;
    memcpy(&v, x + 20, 4);  c->mtid = (int32_t)v;
    memcpy(&v, x + 24, 4);  c->mpos = (int32_t)v;
    memcpy(&v, x + 28, 4);  c->isize = (int32_t)v;

    size_t body = len - 32;
    if (c->l_qname < 1 || c->l_qseq < 0 || body < c->l_qname) return -4;
    if (x[32 + c->l_qname - 1] != '\0') return -4;  /* fgumi always NUL-terminates */
    uint64_t new_l_data = body + c->l_extranul;
    if (new_l_data > INT_MAX) return -4;
    if (((uint64_t)c->n_cigar << 2) + c->l_qname + c->l_extranul
        + (((uint64_t)c->l_qseq + 1) >> 1) + c->l_qseq > new_l_data)
        return -4;
    if ((uint64_t)b->m_data < new_l_data) {
        uint8_t *p = realloc(b->data, new_l_data);
        if (!p) return -4;
        b->data = p;
        b->m_data = new_l_data;
    }
    b->l_data = (int)new_l_data;
    memcpy(b->data, x + 32, c->l_qname);
    memset(b->data + c->l_qname, 0, c->l_extranul);
    memcpy(b->data + c->l_qname + c->l_extranul, x + 32 + c->l_qname,
           body - c->l_qname);
    c->l_qname += c->l_extranul;
    return fgumi_tag2cigar(b) < 0 ? -4 : 0;
}

/*
 * Writer options, applied right after opening (before the header).
 * `no_ref`: store bases verbatim instead of against a reference (valid for any
 * data; mapped reads compress worse). `use_bzip2`/`use_lzma`: allow those
 * block codecs. Returns 0 or -1.
 */
int fgumi_cram_writer_options(htsFile *fp, int no_ref, int use_bzip2, int use_lzma) {
    if (no_ref && hts_set_opt(fp, CRAM_OPT_NO_REF, 1) != 0) return -1;
    if (use_bzip2 && hts_set_opt(fp, CRAM_OPT_USE_BZIP2, 1) != 0) return -1;
    if (use_lzma && hts_set_opt(fp, CRAM_OPT_USE_LZMA, 1) != 0) return -1;
    return 0;
}

/*
 * Set the 20-byte file-definition ID (htslib uses the name the file was opened
 * under; a staged output wants its final name). Must precede the header.
 */
void fgumi_cram_set_file_id(cram_fd *fd, const char *id) {
    memset(fd->file_def->file_id, 0, 20);
    strncpy(fd->file_def->file_id, id, 20);
}

/*
 * Write one raw BAM record (no block_size prefix) through sam_write1, with the
 * checks bam_read1/sam_read1 apply to a BAM record: reference ids within the
 * header, and a mapped CIGAR's query length equal to the sequence length.
 * Unlike bam_read1, the BIN is not recomputed: CRAM does not store it, and
 * neither its encoder nor its index (CRAI) reads it.
 */
int fgumi_hts_write_raw(htsFile *fp, const sam_hdr_t *h, bam1_t *scratch,
                        const uint8_t *x, size_t len) {
    int r = fgumi_bam_from_raw(scratch, x, len);
    if (r < 0) return r;
    bam1_core_t *c = &scratch->core;
    int n_targets = sam_hdr_nref(h);
    if (c->tid < -1 || c->tid >= n_targets || c->mtid < -1 || c->mtid >= n_targets) {
        hts_log_error("fgumi: record %s has a reference id outside the header",
                      bam_get_qname(scratch));
        return -5;
    }
    if (c->n_cigar > 0) {
        hts_pos_t qlen = bam_cigar2qlen(c->n_cigar, bam_get_cigar(scratch));
        if (c->l_qseq > 0 && !(c->flag & BAM_FUNMAP) && qlen != c->l_qseq) {
            hts_log_error("fgumi: CIGAR and query sequence lengths differ for %s",
                          bam_get_qname(scratch));
            return -6;
        }
    }
    return sam_write1(fp, h, scratch) < 0 ? -1 : 0;
}

/*
 * Abandon a CRAM being written (an error path): drop the open container and
 * mark the file abandoned, so cram_close writes no EOF block and the file
 * cannot pass for complete (the container-executor patch refuses the EOF while
 * any handed-over container is unwritten, or once the file is abandoned).
 */
void fgumi_cram_abandon(cram_fd *fd) {
    if (fd->ctr_mt && fd->ctr_mt != fd->ctr)
        cram_free_container(fd->ctr_mt);
    if (fd->ctr)
        cram_free_container(fd->ctr);
    fd->ctr = fd->ctr_mt = NULL;
    fd->executor_abandoned = 1;
}

/* Header text accessors (sam_hdr_str is not const-correct for Rust). */
const char *fgumi_hdr_text(sam_hdr_t *h, size_t *len) {
    const char *s = sam_hdr_str(h);
    *len = s ? sam_hdr_length(h) : 0;
    return s;
}

#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ != __ORDER_LITTLE_ENDIAN__
#error "fgumi_hts.c assumes a little-endian host"
#endif
