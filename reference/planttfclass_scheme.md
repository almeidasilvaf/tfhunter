# Correspondence table between TAPscan and Plant-TFClass

Table used by [`classify_tfs()`](classify_tfs.md) to infer Plant-TFClass
superclasses, classes, and families from TAPscan families and
subfamilies. TAPscan families/subfamilies with multiple rows correspond
to more than one Plant-TFClass family, and the PlantTFDB subfamily of
each gene (variable **PlantTFDB_subfamily**) is used to choose among
them.

## Usage

``` r
data(planttfclass_scheme)
```

## Format

A data frame with the following variables:

- TAPscan_family:

  TAPscan family.

- TAPscan_subfamily:

  TAPscan subfamily (or family, for families without subfamilies).

- PlantTFClass_superclass:

  Plant-TFClass superclass.

- PlantTFClass_class:

  Plant-TFClass class.

- PlantTFClass_family:

  Plant-TFClass family (or class, for classes without families).

- PlantTFDB_subfamily:

  PlantTFDB subfamily used to choose among multiple Plant-TFClass
  families for the same TAPscan family/subfamily. NA for TAPscan
  families/subfamilies with a single Plant-TFClass family.

## References

Blanc-Mathieu, R., Dumas, R., Turchi, L., Lucas, J., & Parcy, F. (2024).
Plant-TFClass: a structural classification for plant transcription
factors. Trends in Plant Science, 29(1), 40-51.

## Examples

``` r
data(planttfclass_scheme)
```
