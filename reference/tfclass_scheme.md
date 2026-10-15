# Correspondence table between CIS-BP and TFClass

Table used by [`classify_tfs()`](classify_tfs.md) to infer TFClass
superclasses, classes, and families from CIS-BP families. The table was
created by comparing the CIS-BP families of human TFs in CIS-BP (build
3.10) with their classifications in TFClass. CIS-BP families were
assigned to the TFClass class of most of their TFs, and to a TFClass
family only if= 90\\

## Usage

``` r
data(tfclass_scheme)
```

## Format

A data frame with the following variables:

- CISBP_family:

  CIS-BP family.

- TFClass_superclass:

  TFClass superclass.

- TFClass_class:

  TFClass class.

- TFClass_family:

  TFClass family (NA if the CIS-BP family corresponds to more than one
  TFClass family).

- TFClass_ID:

  TFClass identifier of the most specific classification level (class or
  family).

## References

Wingender, E., Schoeps, T., Haubrock, M., Krull, M., & Dönitz, J.
(2018). TFClass: expanding the classification of human transcription
factors to their mammalian orthologs. Nucleic Acids Research, 46(D1),
D343-D347.

## Examples

``` r
data(tfclass_scheme)
```
