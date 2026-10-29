# Identify TFs and classify them in families

TFs are classified using the classification schemes of PlantTFDB,
TAPscan, Plant-TFClass, CIS-BP, and TFClass (see
[`list_schemes()`](list_schemes.md)).

## Usage

``` r
classify_tfs(domain_annotation = NULL, scheme = "all", min_shared = 0.1)
```

## Arguments

- domain_annotation:

  A list of data frames with domain annotation, as returned by
  [`annotate_domains()`](annotate_domains.md).

- scheme:

  Character indicating which classification scheme(s) to return. Either
  'all' (all schemes for which domain annotation is available) or one or
  more of 'PlantTFDB', 'TAPscan', 'PlantTFClass', 'CISBP', and
  'TFClass'. Default: 'all'.

- min_shared:

  Numeric indicating the minimum fraction of 5-mers shared with a
  reference TF to assign a C2H2 ZF TF to a single TFClass family (see
  details). Only used for TFClass. Default: 0.1.

## Value

A data frame with one row per gene and the following variables (if
**scheme = "all"**):

- Gene:

  Gene ID.

- PlantTFDB_family:

  TF family according to PlantTFDB.

- PlantTFDB_subfamily:

  TF subfamily according to PlantTFDB.

- TAPscan_family:

  TF family according to TAPscan.

- TAPscan_subfamily:

  TF subfamily according to TAPscan.

- TAP_class:

  Class of transcription-associated protein (TAP) according to TAPscan:
  transcription factor (TF), transcriptional regulator (TR), or putative
  TAP (PT). See `data(tap_class)`.

- PlantTFClass_superclass:

  TF superclass according to Plant-TFClass.

- PlantTFClass_class:

  TF class according to Plant-TFClass.

- PlantTFClass_family:

  TF family according to Plant-TFClass.

- CISBP_family:

  TF family according to CIS-BP.

- TFClass_superclass:

  TF superclass according to TFClass.

- TFClass_class:

  TF class according to TFClass.

- TFClass_family:

  TF family according to TFClass.

Only genes classified by at least one scheme are included. Genes not
classified by a scheme have NA in the respective variables. Genes
assigned to more than one PlantTFDB family have families (and
subfamilies, in the same order) separated by ';'.

If only some schemes are used, only variable **Gene** and the variables
of the chosen schemes are returned, and only genes classified by at
least one of the chosen schemes are included. Note that TAPscan
classifies not only TFs, but also transcriptional regulators and
putative TAPs, which can be distinguished using the variable
**TAP_class**.

## Details

PlantTFDB and TAPscan classifications are based on the domains found in
each gene. Plant-TFClass is a structure-based classification of TFs into
superclasses, classes, and families, and Plant-TFClass classifications
are inferred from TAPscan classifications using a correspondence table
between the two classification systems. When a TAPscan subfamily
corresponds to more than one Plant-TFClass family (e.g., TAPscan's AP2
includes Plant-TFClass's AP2, ERF/DREB, and RAV), PlantTFDB subfamilies
are used to resolve the ambiguity. If the ambiguity cannot be resolved
(e.g., the gene was not classified by PlantTFDB), all possible
classifications are reported, separated by ';'.

CIS-BP classifications are based on the set of DNA-binding domains
(DBDs) found in each gene. TFClass is a classification of TFs into
superclasses, classes, families, and subfamilies based on the structure
of their DBDs, and TFClass classifications are inferred from CIS-BP
families using a correspondence table between the two classification
systems. As TFClass includes only TFs of human and other mammals, it
should be used only for mammals (use CIS-BP for other eukaryotes). TFs
of CIS-BP families that do not occur in mammals (e.g., most plant and
fungal TF families) are not classified by TFClass. CIS-BP families with
multiple types of DBDs that are not in the correspondence table are
classified based on each DBD's family; if these are in different TFClass
classes or families, all possible superclasses, classes, and families
are reported, separated by ';'.

Some CIS-BP families correspond to more than one TFClass family (e.g.,
CIS-BP's Homeodomain corresponds to TFClass families HOX, NK, TALE-type
HD, and others). TFs in these families are assigned to a single TFClass
family using reference TFs of human and mouse with known TFClass
families, as follows:

- **C2H2 ZF**: TFs are assigned to the family of the reference C2H2 ZF
  TF that shares the largest fraction of the TF's amino acid 5-mers
  (ignoring 5-mers found in \> 2% of the reference TFs), if this
  fraction is \>= **min_shared**.

- **Other families**: TFs are scanned with profile HMMs built from the
  DNA-binding domains of the reference TFs in each TFClass family, and
  assigned to the family with the highest score, if this score is \>= 10
  bits higher than the score of any other possible family.

TFs that cannot be assigned to a single family (e.g., C2H2 ZF TFs that
are specific to some lineages, which are not classified by TFClass) have
all possible TFClass families separated by ';'. In validation tests with
TFs of cow, dog, and opossum classified by TFClass, TFs assigned to a
single family were assigned to the correct family in \>99% of the cases.

In all schemes, classification levels without a more specific level are
filled with the name of the more general level. For example, PlantTFDB
family bHLH has no subfamilies, so its subfamily is bHLH.

The classification schemes are available in the data sets
`planttfdb_scheme` (PlantTFDB), `tapscan_scheme` (TAPscan),
`planttfclass_scheme` (correspondence between TAPscan and
Plant-TFClass), `cisbp_scheme` (CIS-BP), and `tfclass_scheme`
(correspondence between CIS-BP and TFClass).

## References

Jin, J., Tian, F., Yang, D. C., Meng, Y. Q., Kong, L., Luo, J., & Gao,
G. (2017). PlantTFDB 4.0: toward a central hub for transcription factors
and regulatory interactions in plants. Nucleic Acids Research, 45(D1),
D1040-D1045.

Petroll, R., Varshney, D., Hiltemann, S., Finke, H., Schreiber, M., de
Vries, J., & Rensing, S. A. (2025). Enhanced sensitivity of TAPscan v4
enables comprehensive analysis of streptophyte transcription factor
evolution. The Plant Journal, 121(1), e17184.

Blanc-Mathieu, R., Dumas, R., Turchi, L., Lucas, J., & Parcy, F. (2024).
Plant-TFClass: a structural classification for plant transcription
factors. Trends in Plant Science, 29(1), 40-51.

Weirauch, M. T., Yang, A., Albu, M., Cote, A. G., Montenegro-Montero,
A., Drewe, P., ... & Hughes, T. R. (2014). Determination and inference
of eukaryotic transcription factor sequence specificity. Cell, 158(6),
1431-1443.

Wingender, E., Schoeps, T., Haubrock, M., Krull, M., & Dönitz, J.
(2018). TFClass: expanding the classification of human transcription
factors to their mammalian orthologs. Nucleic Acids Research, 46(D1),
D343-D347.

## Examples

``` r
data(gsu_annotation)
families <- classify_tfs(gsu_annotation)
head(families)
#>         Gene     PlantTFDB_family PlantTFDB_subfamily TAPscan_family
#> 1 AIG92585.1                 <NA>                <NA>            HAT
#> 2 gsu04730.1             Nin-like            Nin-like         RWP-RK
#> 3 gsu05000.1             Nin-like            Nin-like         RWP-RK
#> 4 gsu06140.1 GARP;MYB superfamily G2-like;MYB-related            MYB
#> 5 gsu06730.1                 <NA>                <NA>           bZIP
#> 6 gsu07030.1                 <NA>                <NA>            PHD
#>   TAPscan_subfamily TAP_class  PlantTFClass_superclass
#> 1              GNAT        TR                     <NA>
#> 2               RKD        TF Helix-turn-helix domains
#> 3               RKD        TF Helix-turn-helix domains
#> 4       MYB-related        TF Helix-turn-helix domains
#> 5              bZIP        TF            Basic domains
#> 6               PHD        TR                     <NA>
#>                    PlantTFClass_class PlantTFClass_family CISBP_family
#> 1                                <NA>                <NA>         <NA>
#> 2                              RWP-RK              RWP-RK       RWP-RK
#> 3                              RWP-RK              RWP-RK       RWP-RK
#> 4          Tryptophan cluster factors         MYB-related     Myb/SANT
#> 5 Basic leucine zipper factors (bZIP)                bZIP         bZIP
#> 6                                <NA>                <NA>         <NA>
#>         TFClass_superclass                       TFClass_class
#> 1                     <NA>                                <NA>
#> 2                     <NA>                                <NA>
#> 3                     <NA>                                <NA>
#> 4 Helix-turn-helix domains          Tryptophan cluster factors
#> 5            Basic domains Basic leucine zipper factors (bZIP)
#> 6                     <NA>                                <NA>
#>                                                                                          TFClass_family
#> 1                                                                                                  <NA>
#> 2                                                                                                  <NA>
#> 3                                                                                                  <NA>
#> 4                                                                                              Myb-SANT
#> 5 Jun-related;Fos-related;Maf-related;B-ATF-related;XBP1-related;ATF4-related;CREB-related;CEBP-related
#> 6                                                                                                  <NA>

# Only TAPscan classifications, with TAP classes (TF, TR, or PT)
tapscan_families <- classify_tfs(gsu_annotation, scheme = "TAPscan")
table(tapscan_families$TAP_class)
#> 
#>  PT  TF  TR 
#>   9 143  96 

# Only CIS-BP classifications
cisbp_families <- classify_tfs(gsu_annotation, scheme = "CISBP")
head(cisbp_families)
#>          Gene CISBP_family
#> 1  gsu04730.1       RWP-RK
#> 2  gsu05000.1       RWP-RK
#> 3  gsu06140.1     Myb/SANT
#> 4  gsu06730.1         bZIP
#> 5  gsu08760.1           DM
#> 6 gsu100410.1       RWP-RK

# TFClass classifications (inferred from CIS-BP families)
tfclass_families <- classify_tfs(gsu_annotation, scheme = "TFClass")
head(tfclass_families)
#>          Gene                    TFClass_superclass
#> 1  gsu06140.1              Helix-turn-helix domains
#> 2  gsu06730.1                         Basic domains
#> 3  gsu08760.1 Zinc-coordinating DNA-binding domains
#> 4 gsu100590.1              Helix-turn-helix domains
#> 5 gsu102830.1              Helix-turn-helix domains
#> 6 gsu105330.1                         Basic domains
#>                             TFClass_class
#> 1              Tryptophan cluster factors
#> 2     Basic leucine zipper factors (bZIP)
#> 3 DM-type intertwined zinc finger factors
#> 4              Tryptophan cluster factors
#> 5                      Heat shock factors
#> 6     Basic leucine zipper factors (bZIP)
#>                                                                                          TFClass_family
#> 1                                                                                              Myb-SANT
#> 2 Jun-related;Fos-related;Maf-related;B-ATF-related;XBP1-related;ATF4-related;CREB-related;CEBP-related
#> 3                                                                                                  DMRT
#> 4                                                                                              Myb-SANT
#> 5                                                                                                   HSF
#> 6 Jun-related;Fos-related;Maf-related;B-ATF-related;XBP1-related;ATF4-related;CREB-related;CEBP-related
```
