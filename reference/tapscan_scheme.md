# Data frame of TAPscan's classification scheme

TAPscan v4 classification rules, with one row per rule. Genes are
assigned to a family (and subfamily, if any) if they have all domains in
variable **Should** and none of the domains in variable **Should_not**.
For families with more than one rule (e.g., bZIP), genes are assigned to
the family if they satisfy any of the rules.

## Usage

``` r
data(tapscan_scheme)
```

## Format

A data frame with the following variables:

- Family:

  TAPscan family, as in variable **TAPscan_family** of the output of
  [`classify_tfs()`](classify_tfs.md).

- Subfamily:

  TAPscan subfamily, as in variable **TAPscan_subfamily** of the output
  of [`classify_tfs()`](classify_tfs.md). Families without subfamilies
  have the family name as subfamily.

- Rule:

  Name of the rule in TAPscan's rules file.

- Should:

  Domains that genes should have, separated by ' and '. Numbers in
  parentheses (e.g., 'Myb_DNA-binding (2)') indicate the required number
  of copies of a domain.

- Should_not:

  Domains that genes should not have, separated by ' or '.

- TAP_class:

  TAP class. One of 'TF' (transcription factor), 'TR' (transcriptional
  regulator), or 'PT' (putative TAP). NA for families without a TAP
  class in TAPscan (see `tap_class`).

## References

Petroll, R., Varshney, D., Hiltemann, S., Finke, H., Schreiber, M., de
Vries, J., & Rensing, S. A. (2025). Enhanced sensitivity of TAPscan v4
enables comprehensive analysis of streptophyte transcription factor
evolution. The Plant Journal, 121(1), e17184.

## Examples

``` r
data(tapscan_scheme)
```
