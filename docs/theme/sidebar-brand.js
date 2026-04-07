/**
 * Fulcrum Genomics — fgumi sidebar branding, breadcrumbs, footer injection.
 * Logo SVG is inlined to avoid mdBook asset-path resolution issues.
 */
(function () {
    // Selectively clear only sidebar visibility state (not theme preference).
    // Bumping this version resets fold/sidebar state for all returning visitors.
    var SIDEBAR_VERSION = '3';
    if (localStorage.getItem('fg-sidebar-version') !== SIDEBAR_VERSION) {
        localStorage.removeItem('sidebar');
        localStorage.setItem('fg-sidebar-version', SIDEBAR_VERSION);
    }

    // ── Official Fulcrum Genomics logo (inline SVG) ───────────────────────────
    var LOGO_SVG =
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 102.33 31.81" class="fg-brand-logo" aria-label="Fulcrum Genomics">' +
        '<defs><style>.fg-w{fill:#fff}.fg-b{fill:#26a8e0}.fg-g{fill:#38b44a}</style></defs>' +
        '<path class="fg-w" d="M43.41,17.49h1.43v-4.03h3.92v-1.27h-3.92v-2.78h4.59v-1.35h-6.02v9.42z"/>' +
        '<path class="fg-w" d="M51.34,14.29c0,.99.34,1.78.98,2.36.63.6,1.42.9,2.37.91.96,0,1.76-.31,2.38-.91.62-.58.95-1.37.96-2.36v-6.23h-1.43v6.07c0,.64-.19,1.13-.54,1.47-.36.35-.82.53-1.38.53s-1.01-.18-1.36-.53c-.36-.34-.55-.82-.56-1.47v-6.07h-1.43v6.23z"/>' +
        '<path class="fg-w" d="M60.67,17.49h6.02v-1.35h-4.59v-8.07h-1.43v9.42z"/>' +
        '<path class="fg-w" d="M73.54,14.85c-.4.86-1,1.29-1.8,1.29-.34,0-.62-.07-.87-.21-.25-.12-.44-.28-.59-.47-.19-.2-.31-.47-.37-.79-.07-.33-.1-.95-.1-1.88s.03-1.56.1-1.89c.06-.32.18-.58.37-.78.15-.19.35-.36.59-.48.24-.12.53-.19.87-.2.46,0,.84.14,1.16.39.31.27.52.59.63.97h1.51c-.15-.79-.51-1.45-1.09-1.98-.57-.53-1.31-.8-2.22-.81-.74,0-1.36.19-1.85.53-.5.34-.87.73-1.1,1.16-.14.23-.25.53-.32.9-.06.37-.1,1.1-.1,2.2s.03,1.8.1,2.18c.03.2.08.37.13.5.06.13.12.26.19.41.23.44.59.82,1.1,1.15.5.34,1.11.53,1.85.54.82,0,1.53-.23,2.13-.7.58-.47.98-1.14,1.18-2.02h-1.51z"/>' +
        '<path class="fg-w" d="M77.98,9.34h2.24c.46,0,.81.1,1.05.29.31.22.46.57.47,1.07,0,.41-.13.75-.39,1.03-.27.3-.67.46-1.2.47h-2.16v-2.86zM76.55,17.49h1.43v-4.03h1.82l1.94,4.03h1.7l-2.18-4.18c1.2-.46,1.8-1.33,1.82-2.61-.03-.87-.34-1.54-.94-2.01-.5-.41-1.13-.62-1.92-.62h-3.68v9.42z"/>' +
        '<path class="fg-w" d="M84.95,14.29c0,.99.34,1.78.98,2.36.63.6,1.42.9,2.37.91.96,0,1.76-.31,2.38-.91.62-.58.95-1.37.96-2.36v-6.23h-1.43v6.07c0,.64-.19,1.13-.54,1.47-.36.35-.82.53-1.38.53s-1.01-.18-1.36-.53c-.36-.34-.55-.82-.56-1.47v-6.07h-1.43v6.23z"/>' +
        '<path class="fg-w" d="M94.28,17.49h1.43v-5.87h.03l1.97,4.52h1.19l1.97-4.52h.03v5.87h1.43v-9.42h-1.35l-2.65,6.14-2.7-6.14h-1.34v9.42z"/>' +
        '<path class="fg-w" d="M45.06,22.21h.99v.25c0,.29-.1.53-.29.71-.19.19-.42.28-.71.28-.17,0-.32-.04-.45-.11-.13-.06-.23-.14-.31-.24-.1-.1-.16-.23-.19-.4-.04-.16-.05-.48-.05-.94s.02-.78.05-.95c.03-.16.09-.29.19-.39.08-.1.18-.18.31-.24.12-.06.27-.1.45-.1.24,0,.43.07.6.2.16.13.27.29.33.48h.78c-.08-.39-.26-.72-.56-.99-.3-.26-.68-.4-1.14-.4-.38,0-.7.09-.96.26-.26.17-.45.36-.57.58-.07.11-.13.26-.16.45-.03.18-.05.55-.05,1.1s.02.9.05,1.09c.02.1.04.19.07.25.03.06.06.13.1.2.12.22.31.41.57.58.26.17.57.26.96.27.49,0,.9-.17,1.23-.49.32-.32.49-.71.5-1.19v-.96h-1.72v.68z"/>' +
        '<path class="fg-w" d="M51.67,24.12h3.1v-.68h-2.36v-1.38h2.02v-.63h-2.02v-1.34h2.36v-.68h-3.1v4.71z"/>' +
        '<path class="fg-w" d="M59.49,24.12h.74v-3.35h.01l2.19,3.35h.7v-4.71h-.74v3.35h-.01l-2.2-3.35h-.69v4.71z"/>' +
        '<path class="fg-w" d="M68.02,21.77c0,.54.02.9.05,1.09.02.1.04.19.07.25.03.06.06.13.1.2.12.22.31.41.57.58.26.17.57.26.96.27.39,0,.71-.1.96-.27.25-.17.44-.36.55-.58.08-.11.14-.27.17-.46.03-.19.04-.55.04-1.09s-.01-.91-.04-1.1c-.03-.19-.09-.33-.17-.45-.11-.22-.3-.41-.55-.58-.26-.17-.58-.26-.96-.26-.38,0-.7.09-.96.26-.26.17-.45.36-.57.58-.07.11-.13.26-.16.45-.03.18-.05.55-.05,1.1zM68.75,21.77c0-.46.02-.78.05-.95.03-.16.09-.29.19-.39.08-.1.18-.18.31-.24.12-.06.27-.1.45-.1.18,0,.33.04.46.1.12.06.22.15.29.24.1.1.16.23.2.39.03.17.05.48.05.95s-.02.78-.05.94c-.04.16-.1.3-.2.4-.07.1-.17.18-.29.24-.13.07-.28.11-.46.11s-.32-.04-.45-.11c-.13-.06-.23-.14-.31-.24-.1-.1-.16-.23-.19-.4-.04-.16-.05-.48-.05-.94z"/>' +
        '<path class="fg-w" d="M76.37,24.12h.74v-2.94h.01l1.01,2.26h.61l1.01-2.26h.02v2.94h.74v-4.71h-.7l-1.36,3.07-1.39-3.07h-.69v4.71z"/>' +
        '<rect class="fg-w" x="85.56" y="19.41" width=".74" height="4.71"/>' +
        '<path class="fg-w" d="M93.84,22.81c-.2.43-.51.64-.93.64-.17,0-.32-.04-.45-.11-.13-.06-.23-.14-.31-.24-.1-.1-.16-.23-.19-.4-.04-.16-.05-.48-.05-.94s.02-.78.05-.95c.03-.16.09-.29.19-.39.08-.1.18-.18.31-.24.12-.06.27-.1.45-.1.24,0,.43.07.6.2.16.13.27.29.33.48h.78c-.08-.39-.26-.72-.56-.99-.3-.26-.68-.4-1.14-.4-.38,0-.7.09-.95.26-.26.17-.45.36-.57.58-.07.11-.13.26-.16.45-.03.18-.05.55-.05,1.1s.02.9.05,1.09c.02.1.04.19.07.25.03.06.06.13.1.2.12.22.31.41.57.58.25.17.57.26.95.27.42,0,.79-.12,1.09-.35.3-.23.5-.57.61-1.01h-.78z"/>' +
        '<path class="fg-w" d="M99.22,22.98l-.48.54c.52.43,1.13.65,1.85.65,1.11-.01,1.68-.47,1.7-1.37,0-.33-.11-.63-.32-.88-.22-.26-.55-.41-1.01-.47-.23-.03-.41-.05-.55-.07-.24-.04-.41-.12-.52-.23-.11-.11-.16-.23-.16-.37,0-.23.09-.4.24-.51.15-.11.34-.16.57-.16.44,0,.84.13,1.2.36l.41-.59c-.45-.31-.97-.47-1.57-.49-.5,0-.89.13-1.16.38-.28.25-.42.58-.42,1,0,.34.11.63.34.87.22.23.53.38.95.45.23.03.45.06.64.09.43.07.64.28.63.63,0,.43-.33.65-.96.66-.53,0-.99-.16-1.38-.47z"/>' +
        '<path class="fg-g" d="M21.21,12.15h-10.13v5.54h10.13v1.88c0,1.39.12,3.66-1.3,3.66h-8.83v7.51h2.37c2.03,0,3.87-.25,5.52-.75,1.64-.5,3.04-1.24,4.19-2.21,1.15-.97,2.03-2.14,2.66-3.51.15-.33.29-.68.4-1.04.36-1.11.54-2.33.54-3.66v-7.42h-5.54z"/>' +
        '<path class="fg-b" d="M13.31,1.07c-2.08,0-3.94.25-5.59.74-1.64.49-3.04,1.22-4.19,2.18-1.15.96-2.02,2.13-2.63,3.51-.6,1.38-.91,2.96-.91,4.74v5.46h5.54v-5.54c0-.82.11-1.69.44-2.44.26-.59.61-1.06,1.1-1.47.64-.54,1.44-.89,2.23-1.12,1.28-.38,2.67-.51,4-.51h13.45V1.07h-13.45zM0,23.23v7.51h5.54v-7.51H0z"/>' +
        '<rect class="fg-w" x="34.9" width=".42" height="31.81"/>' +
        '</svg>';

    // ── Footer HTML (injected on all non-index pages) ─────────────────────────
    var FOOTER_HTML =
        '<div class="fg-visit-us fg-page-footer">' +
        '<a href="https://www.fulcrumgenomics.com" target="_blank" rel="noopener">' +
        '<picture>' +
        '<source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/fulcrumgenomics/fgumi/main/.github/logos/fulcrumgenomics-light.svg">' +
        '<img alt="Fulcrum Genomics" src="https://raw.githubusercontent.com/fulcrumgenomics/fgumi/main/.github/logos/fulcrumgenomics-dark.svg" height="50">' +
        '</picture>' +
        '</a>' +
        '<p>Visit us at <a href="https://www.fulcrumgenomics.com" target="_blank" rel="noopener">Fulcrum Genomics</a> to learn more about how we can power your Bioinformatics with fgumi and beyond.</p>' +
        '<div class="fg-visit-buttons">' +
        '<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]" class="fg-btn fg-btn-green">&#9993; Email Us</a>' +
        '<a href="https://www.fulcrumgenomics.com" target="_blank" rel="noopener" class="fg-btn fg-btn-blue">&#127760; Visit Us</a>' +
        '</div>' +
        '<div class="fg-footer-links">' +
        '<a href="https://github.com/fulcrumgenomics/fgumi" target="_blank" rel="noopener">GitHub</a>' +
        '<span class="fg-footer-sep">·</span>' +
        '<a href="https://docs.rs/fgumi" target="_blank" rel="noopener">API Docs</a>' +
        '<span class="fg-footer-sep">·</span>' +
        '<a href="https://github.com/fulcrumgenomics/fgumi/issues" target="_blank" rel="noopener">Issues</a>' +
        '<span class="fg-footer-sep">·</span>' +
        '<a href="https://github.com/fulcrumgenomics/fgumi/discussions" target="_blank" rel="noopener">Discussions</a>' +
        '</div>' +
        '</div>';

    // Compute the path prefix to reach the docs root from the current page.
    // All non-index content lives exactly one directory below the docs root
    // (e.g., guide/page.html, tools/page.html, metrics/page.html), so we
    // need "../" from those pages and "" from the root level.
    var docsRoot = (function () {
        var segs = window.location.pathname.split('/').filter(Boolean);
        if (segs.length >= 2) {
            var parent = segs[segs.length - 2];
            if (parent === 'guide' || parent === 'tools' || parent === 'metrics') {
                return '../';
            }
        }
        return '';
    }());

    // ── Sidebar brand injection ───────────────────────────────────────────────
    function injectBrand() {
        var scrollbox = document.querySelector('.sidebar-scrollbox');
        if (!scrollbox || document.querySelector('.fg-brand-header')) return;

        var header = document.createElement('div');
        header.className = 'fg-brand-header';
        // Logo links home; corporate link is separate below
        header.innerHTML =
            '<a href="' + docsRoot + '" class="fg-brand-link">' + LOGO_SVG + '</a>' +
            '<div class="fg-brand-divider"></div>' +
            '<span class="fg-brand-product">' +
            '<span style="color:#26a8e0">fg</span>' +
            '<span style="color:#38b44a">umi</span>' +
            '</span>' +
            '<span class="fg-brand-tagline">UMI Consensus Toolkit</span>' +
            '<a href="https://www.fulcrumgenomics.com" target="_blank" rel="noopener" class="fg-brand-corp">' +
            'Fulcrum Genomics &#8599;' +
            '</a>';

        scrollbox.insertBefore(header, scrollbox.firstChild);
    }

    // ── Footer injection ──────────────────────────────────────────────────────
    function injectFooter() {
        if (document.querySelector('.fg-page-footer')) return;
        var main = document.querySelector('.content main');
        if (!main) return;
        var wrapper = document.createElement('div');
        wrapper.innerHTML = FOOTER_HTML;
        main.appendChild(wrapper.firstChild);
    }

    // ── Breadcrumb injection ──────────────────────────────────────────────────
    var SECTION_MAP = {
        'guide':   { label: 'User Guide',        url: null },
        'tools':   { label: 'Tool Reference',    url: 'tools/index.html' },
        'metrics': { label: 'Metrics Reference', url: 'metrics/index.html' }
    };

    function injectBreadcrumb() {
        var main = document.querySelector('.content main');
        if (!main || document.querySelector('.fg-breadcrumb')) return;

        var path = window.location.pathname;
        var segments = path.replace(/\.html$/, '').split('/').filter(Boolean);
        if (segments.length < 2) return; // root or single-segment: no breadcrumb

        var dir = segments[segments.length - 2];
        var section = SECTION_MAP[dir];
        if (!section) return;

        var h1 = main.querySelector('h1');
        var pageTitle = h1 ? h1.textContent.trim() : segments[segments.length - 1].replace(/-/g, ' ');

        var sep = '<span class="fg-breadcrumb-sep">&#8250;</span>';
        var html = '<a href="' + docsRoot + '">Home</a>' + sep;
        if (section.url) {
            html += '<a href="' + docsRoot + section.url + '">' + section.label + '</a>';
        } else {
            html += '<span>' + section.label + '</span>';
        }
        var currentSpan = document.createElement('span');
        currentSpan.className = 'fg-breadcrumb-current';
        currentSpan.textContent = pageTitle;

        var crumb = document.createElement('nav');
        crumb.className = 'fg-breadcrumb';
        crumb.setAttribute('aria-label', 'breadcrumb');
        crumb.innerHTML = html;
        crumb.appendChild(currentSpan);

        if (h1) {
            main.insertBefore(crumb, h1);
        } else {
            main.insertBefore(crumb, main.firstChild);
        }
    }

    // ── Menu bar title coloring ───────────────────────────────────────────────
    function colorMenuTitle() {
        var title = document.querySelector('.menu-title');
        if (!title || title.querySelector('span')) return;
        var text = title.textContent;
        if (text.startsWith('fg')) {
            title.innerHTML =
                '<span style="color:#26a8e0">fg</span>' +
                '<span style="color:#38b44a">' + text.slice(2) + '</span>' +
                '<span class="fg-byline-inline">by <a href="https://www.fulcrumgenomics.com" target="_blank" rel="noopener">Fulcrum Genomics &#8599;</a></span>';
        }
    }

    // ── Init ──────────────────────────────────────────────────────────────────
    function init() {
        injectBrand();
        injectFooter();
        injectBreadcrumb();
        colorMenuTitle();
    }

    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', init);
    } else {
        init();
    }
})();
