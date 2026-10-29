# Correspondence table between CIS-BP and TFClass

Table used by [`classify_tfs()`](classify_tfs.md) to infer TFClass
superclasses, classes, and families from CIS-BP families. The table was
created by comparing the CIS-BP families of human TFs in CIS-BP (build
3.10) with their classifications in TFClass. CIS-BP families were
assigned to the TFClass class of most of their TFs, and to the TFClass
family of at least 90\\ of their TFs (in this class) are listed,
separated by ';' (e.g., CIS-BP's Homeodomain corresponds to TFClass
families HOX, NK, and others). TFs in these CIS-BP families are assigned
to a single TFClass family by [`classify_tfs()`](classify_tfs.md) when
possible (see [`classify_tfs()`](classify_tfs.md) for details). Some
entries were manually curated, including CIS-BP families with multiple
types of DNA-binding domains (e.g., 'Homeodomain,POU') and CIS-BP's AT
hook, whose TFClass family is the TFClass class 'A.T hook factors' (as
CIS-BP's AT hook TFs include proteins other than those in HMGA, TFClass'
only family in this class). CIS-BP families that are not in this table
(e.g., families that do not occur in mammals, or with DNA-binding
domains that are not classified by TFClass) are not classified by
TFClass, except for families with multiple types of DNA-binding domains
(see Details in [`?classify_tfs`](classify_tfs.md)).

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

  TFClass family, or all possible families separated by ';' if the
  CIS-BP family corresponds to more than one TFClass family. For CIS-BP
  families corresponding to a TFClass class but not to a specific family
  (i.e., AT hook), the TFClass class.

- TFClass_ID:

  TFClass identifier(s) of the most specific classification level (class
  or family), separated by ';' if there are multiple families.

## References

Wingender, E., Schoeps, T., Haubrock, M., Krull, M., & Dönitz, J.
(2018). TFClass: expanding the classification of human transcription
factors to their mammalian orthologs. Nucleic Acids Research, 46(D1),
D343-D347.

## Examples

``` r
data(tfclass_scheme)
```
