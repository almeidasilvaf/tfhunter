# List available TF classification schemes

List available TF classification schemes

## Usage

``` r
list_schemes()
```

## Value

A data frame with one row per classification scheme and the following
variables:

- Scheme:

  Name of the classification scheme, to be passed to the argument
  **scheme** of [`annotate_domains()`](annotate_domains.md),
  [`classify_tfs()`](classify_tfs.md), and
  [`get_tf_counts()`](get_tf_counts.md).

- Taxa:

  Taxonomic group for which the scheme was designed.

- Levels:

  Classification levels, from the most general to the most specific.

- Annotation:

  Domain annotation(s) required by the scheme (i.e., elements of the
  output of [`annotate_domains()`](annotate_domains.md)), separated by
  ', '.

- Reference:

  Reference of the classification scheme.

## Examples

``` r
list_schemes()
#>         Scheme             Taxa                    Levels         Annotation
#> 1    PlantTFDB           Plants         family, subfamily          PlantTFDB
#> 2      TAPscan Plants and algae         family, subfamily            TAPscan
#> 3 PlantTFClass           Plants superclass, class, family PlantTFDB, TAPscan
#> 4        CISBP       Eukaryotes                    family              CISBP
#> 5      TFClass          Mammals superclass, class, family     CISBP, TFClass
#>                                                        Reference
#> 1                      Jin et al. (2017), doi:10.1093/nar/gkw982
#> 2                   Petroll et al. (2025), doi:10.1111/tpj.17184
#> 3 Blanc-Mathieu et al. (2024), doi:10.1016/j.tplants.2023.06.023
#> 4         Weirauch et al. (2014), doi:10.1016/j.cell.2014.08.009
#> 5                Wingender et al. (2018), doi:10.1093/nar/gkx987
```
