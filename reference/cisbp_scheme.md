# Data frame of CIS-BP's TF family classification scheme

DNA-binding domains (DBDs) used by CIS-BP to identify TFs, and the TF
family corresponding to each DBD. Genes with more than one type of DBD
are assigned to a family whose name is the family names of all DBDs
sorted alphabetically (case-insensitive) and separated by ',' (e.g.,
'Homeodomain,POU'), as in CIS-BP. DBDs and families were obtained from
CIS-BP (build 3.10), and DBDs were matched to Pfam (release 38.2)
accessions. Some Pfam families were renamed since they were used by
CIS-BP (e.g., 'DUF260' is now 'LOB'), and the Pfam family 'FLO_LFY' was
split into 'SAM_LFY' (PF01698) and 'C_LFY_FLO' (PF17538), both of which
are used for family 'Leafy'. The RNA-binding domains 'RRM_1' and
'zf-CCCH' are not included, as CIS-BP classifies TFs in families 'RRM'
and 'CCCH ZF' only if there is experimental evidence of DNA binding (6
TFs in CIS-BP 3.10), and never uses these domains for TFs with other
DBDs.

## Usage

``` r
data(cisbp_scheme)
```

## Format

A data frame with the following variables:

- Family:

  CIS-BP family of TFs with the DBD.

- CISBP_DBD:

  DBD name in CIS-BP.

- Pfam:

  Pfam accession of the DBD.

- Pfam_name:

  Pfam name of the DBD (Pfam 38.2).

- Pfam_description:

  Pfam description of the DBD.

## References

Weirauch, M. T., Yang, A., Albu, M., Cote, A. G., Montenegro-Montero,
A., Drewe, P., ... & Hughes, T. R. (2014). Determination and inference
of eukaryotic transcription factor sequence specificity. Cell, 158(6),
1431-1443.

## Examples

``` r
data(cisbp_scheme)
```
