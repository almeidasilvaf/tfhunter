# Curated human transcription factors and their classifications

Human TFs curated by Lambert et al. (2018), classified with the CIS-BP
and TFClass schemes. As this list was manually curated by experts, it is
more accurate than TFs identified from the human proteome with
[`annotate_domains()`](annotate_domains.md) and
[`classify_tfs()`](classify_tfs.md), and it can be used as the set of
human TFs in comparative analyses (e.g., with
[`get_tf_counts()`](get_tf_counts.md)). CIS-BP families were obtained
from CIS-BP (build 3.10), whose human TFs are the TFs in this list.
TFClass classifications were obtained from TFClass for TFs in TFClass,
and inferred from CIS-BP families (as in
[`classify_tfs()`](classify_tfs.md)) for other TFs.

## Usage

``` r
data(human_tfs)
```

## Format

A data frame with one row per TF and the following variables:

- Gene:

  Ensembl gene ID (or UniProt entry name, for DUX1 and DUX3, which are
  not in Ensembl).

- Symbol:

  HGNC gene symbol.

- CISBP_family:

  TF family according to CIS-BP. TFs without a known DNA-binding domain
  are in family 'Unknown'.

- TFClass_superclass:

  TF superclass according to TFClass.

- TFClass_class:

  TF class according to TFClass.

- TFClass_family:

  TF family according to TFClass.

- TFClass_source:

  Source of TFClass classifications. One of 'TFClass' (TF classified by
  TFClass) or 'Inferred from CIS-BP' (TFClass classifications inferred
  from the CIS-BP family). NA for TFs without TFClass classifications.

- DBD:

  DNA-binding domain(s), as in Lambert et al. (2018).

- TF_assessment:

  Assessment of DNA-binding activity. One of 'Known motif', 'Inferred
  motif' (motif inferred from a closely related TF), or 'Likely to be
  sequence specific TF'.

- Binding_mode:

  Mode of interacting with DNA. One of 'Monomer or homomultimer',
  'Obligate heteromer', or 'Low specificity DNA-binding protein'.

- Motif_status:

  Availability of DNA-binding motifs. One of 'High-throughput in vitro',
  '100 perc ID - in vitro' (motif available for a TF with an identical
  DNA-binding domain), 'In vivo/Misc source', or 'No motif'.

## References

Lambert, S. A., Jolma, A., Campitelli, L. F., Das, P. K., Yin, Y., Albu,
M., ... & Weirauch, M. T. (2018). The human transcription factors. Cell,
172(4), 650-665.

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
data(human_tfs)
head(human_tfs)
#>              Gene     Symbol CISBP_family                    TFClass_superclass
#> 1 ENSG00000267179 AC008770.3      C2H2 ZF Zinc-coordinating DNA-binding domains
#> 2 ENSG00000267281 AC023509.3 bZIP,C2H2 ZF                         Basic domains
#> 3 ENSG00000233757 AC092835.1      C2H2 ZF Zinc-coordinating DNA-binding domains
#> 4 ENSG00000264668 AC138696.1      C2H2 ZF Zinc-coordinating DNA-binding domains
#> 5 ENSG00000101126       ADNP  Homeodomain              Helix-turn-helix domains
#> 6 ENSG00000101544      ADNP2  Homeodomain              Helix-turn-helix domains
#>                         TFClass_class
#> 1            C2H2 zinc finger factors
#> 2 Basic leucine zipper factors (bZIP)
#> 3            C2H2 zinc finger factors
#> 4            C2H2 zinc finger factors
#> 5                Homeo domain factors
#> 6                Homeo domain factors
#>                                                                                                                                                          TFClass_family
#> 1 Three-zinc finger Krüppel-related;Other factors with up to three adjacent zinc fingers;More than 3 adjacent zinc fingers;Factors with multiple dispersed zinc fingers
#> 2                                                                                                                                                           Jun-related
#> 3 Three-zinc finger Krüppel-related;Other factors with up to three adjacent zinc fingers;More than 3 adjacent zinc fingers;Factors with multiple dispersed zinc fingers
#> 4                                                                                                                                     More than 3 adjacent zinc fingers
#> 5                                                                                                                                                                 HD-ZF
#> 6                                                                                                                                                                 HD-ZF
#>         TFClass_source         DBD                     TF_assessment
#> 1 Inferred from CIS-BP     C2H2 ZF Likely to be sequence specific TF
#> 2              TFClass        bZIP                       Known motif
#> 3 Inferred from CIS-BP     C2H2 ZF Likely to be sequence specific TF
#> 4              TFClass     C2H2 ZF                       Known motif
#> 5              TFClass Homeodomain Likely to be sequence specific TF
#> 6              TFClass Homeodomain Likely to be sequence specific TF
#>              Binding_mode           Motif_status
#> 1 Monomer or homomultimer               No motif
#> 2 Monomer or homomultimer 100 perc ID - in vitro
#> 3 Monomer or homomultimer               No motif
#> 4 Monomer or homomultimer 100 perc ID - in vitro
#> 5 Monomer or homomultimer               No motif
#> 6 Monomer or homomultimer               No motif

# Count human TFs per TFClass class
se <- get_tf_counts(
    list(Human = human_tfs), scheme = "TFClass", level = "class"
)
```
