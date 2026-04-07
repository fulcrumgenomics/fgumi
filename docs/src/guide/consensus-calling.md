# Calling Consensus Reads

## Overview

Reads with the same molecular identifier (MI tag) are examined base-by-base to determine the most likely base in the original source molecule. The consensus calling model has three steps:

1. Adjusting input base qualities
2. Computing the maximum posterior probability base
3. Adjusting the output consensus base quality

## Glossary

| Symbol | Description |
|--------|-------------|
| Q | Phred-scaled base quality for a single base (measures sequencing error) |
| S_Q | Value subtracted from input base qualities (prior to capping) |
| M_Q | Maximum base quality cap (applied after shifting) |
| Err_pre | Phred-scaled error rate for errors *before* UMI integration (e.g. deamination, oxidation during library prep) |
| Err_post | Phred-scaled error rate for errors *after* UMI integration but before sequencing (e.g. amplification, target capture) |
| B_i | The base of the i-th read at a given position |

## Step 1: Adjusting Input Base Qualities

Base qualities are assumed to represent the probability of a sequencing error. Two optional adjustments are applied:

1. **Shift**: Subtract a fixed value from the phred-scaled qualities (e.g., Q30 with shift of 10 becomes Q20)
2. **Cap**: Limit to a maximum phred-scaled value

```text
Q' = min(Q - S_Q, M_Q)
```

These adjustments should only be used if input base qualities are systematically over-estimated.

The adjusted quality is converted to an error probability:

```text
P_Q' = 10^(-Q'/10)
```

Then combined with the post-UMI error rate to produce a compound error probability covering all processes from UMI integration through sequencing:

```text
P_Q'' = Err_post * (1 - P_Q') + (1 - Err_post) * P_Q' + (Err_post * P_Q' * 2/3)
```

This formula sums three terms:
1. Error in post-UMI processes, no sequencing error
2. No post-UMI error, but sequencing error
3. Both errors occur, but the second doesn't reverse the first (probability 2/3 for DNA with 4 bases)

## Step 2: Computing the Consensus Base

For each position, the likelihood that the true base is A, C, G, or T is computed by multiplying across all reads:

```
L(Call=B) = ∏_i { P_Q''/3  if B ≠ B_i
                 { 1 - P_Q'' if B = B_i
```

The likelihoods are normalized to posterior probabilities (assuming a uniform prior):

```
Post(Call=B) = L(Call=B) / Σ L(Call=C) for C in {A, C, G, T}
```

The base with the maximum posterior probability becomes the consensus call.

## Step 3: Adjusting Output Quality

The consensus posterior is converted to an error probability and then modified to incorporate the pre-UMI error rate (errors before UMI integration, such as deamination or oxidation):

```
Pr_err = 1 - Post(Call)
Pr_err' = Err_pre * (1 - Pr_err) + (1 - Err_pre) * Pr_err + (Err_pre * Pr_err * 2/3)
Q_call = -10 * log10(Pr_err')
```

The final consensus base quality represents the probability of error across the entire process: from sample extraction through library preparation, UMI integration, amplification, and sequencing.

Any consensus base with quality below the minimum threshold is masked to `N`.

## Caveats

- Each end of a pair is treated independently; overlapping bases within a pair are jointly called by default (disable with `--consensus-call-overlapping-bases false`)
- Indel errors in the reads are not considered in the consensus model
- `simplex` and `codec` do not accept a `--sort-order` flag; consensus reads are emitted as
  unmapped and should be sorted by the downstream pipeline (`fgumi zipper` + `fgumi sort`)
