# Classes of transcription-associated proteins (TAPs) in TAPscan

Families of transcription-associated proteins (TAPs) in TAPscan v4,
their classes, descriptions, and references. Data were obtained from
TAPscan's website (file *\_data/import-tapinfo/tapinfo_v4.csv* in the
GitHub repository Rensing-Lab/TAPscan-v4-website). HTML tags in
descriptions and references were removed. Family MADS_MIKC, which is not
included in the original file, was manually added as a TF family (as
family MADS), without description and references.

## Usage

``` r
data(tap_class)
```

## Format

A data frame with the following variables:

- TAP_family:

  TAP family or subfamily, as in the variable **TAPscan_subfamily** of
  the output of [`classify_tfs()`](classify_tfs.md).

- TAP_class:

  TAP class. One of 'TF' (transcription factor), 'TR' (transcriptional
  regulator), or 'PT' (putative TAP).

- Description:

  Description of the TAP family.

- References:

  References for the TAP family, separated by ' \| '.

## References

Petroll, R., Varshney, D., Hiltemann, S., Finke, H., Schreiber, M., de
Vries, J., & Rensing, S. A. (2025). Enhanced sensitivity of TAPscan v4
enables comprehensive analysis of streptophyte transcription factor
evolution. The Plant Journal, 121(1), e17184.

## Examples

``` r
data(tap_class)
```
