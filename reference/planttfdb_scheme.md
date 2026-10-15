# Data frame of PlantTFDB's TF family classification scheme

The classification scheme is the same as the one used by PlantTFDB,
except for three family/subfamily names that were changed to match the
output of [`classify_tfs()`](classify_tfs.md) ('M_type' to 'M-type',
'MYB_related' to 'MYB-related', and 'LBD (AS2/LOB)' to 'LBD'). As in
PlantTFDB, families without subfamilies have the family name as
subfamily.

## Usage

``` r
data(planttfdb_scheme)
```

## Format

A data frame with the following variables:

- Family:

  TF family name.

- Subfamily:

  TF subfamily name.

- DBD:

  DNA-binding domain

- Auxiliary:

  Auxiliary domain

- Forbidden:

  Forbidden domain

## References

Jin, J., Tian, F., Yang, D. C., Meng, Y. Q., Kong, L., Luo, J., & Gao,
G. (2016). PlantTFDB 4.0: toward a central hub for transcription factors
and regulatory interactions in plants. Nucleic acids research, gkw982.

## Examples

``` r
data(planttfdb_scheme)
```
