# Annotate protein sequences with domains used for TF classification

Protein sequences are scanned with HMMER for the domains used by the
classification schemes of choice (see
[`list_schemes()`](list_schemes.md)). Each domain annotation set has its
own profile HMMs and filtering criteria:

- **PlantTFDB**: PFAM and self-built profile HMMs, with the
  domain-specific score cutoffs defined by PlantTFDB (or an E-value
  threshold, for domains without cutoffs).

- **TAPscan**: TAPscan v4 profile HMMs, with their gathering (GA)
  thresholds, as in the original TAPscan implementation.

- **CISBP**: Pfam profile HMMs of the DNA-binding domains used by
  CIS-BP, with the thresholds used by CIS-BP (per-sequence E-value \<
  0.01 and per-domain conditional E-value \< 0.01). As in CIS-BP,
  overlapping hits of different domains are resolved by keeping only the
  most significant hit (lowest independent E-value).

## Usage

``` r
annotate_domains(seq = NULL, scheme = "all", evalue = 1e-05, threads = 1)
```

## Arguments

- seq:

  An `AAStringSet` object. The sequences in this object must represent
  only the translated sequences of primary (or longest) transcripts. To
  create such an object from a FASTA file, use
  [`Biostrings::readAAStringSet()`](https://rdrr.io/pkg/Biostrings/man/XStringSet-io.html).

- scheme:

  Character vector indicating the classification scheme(s) for which
  domains will be annotated. One or more of the schemes in
  [`list_schemes()`](list_schemes.md), or 'all' to annotate domains for
  all schemes. Default: 'all'.

- evalue:

  Numeric indicating the E-value threshold for `hmmsearch` to be used
  for PlantTFDB domains without pre-defined domain cutoffs. Default:
  1e-05.

- threads:

  Numeric indicating the number of threads to be used by `hmmsearch`.
  Default: 1.

## Value

A list of data frames with domain annotation for each domain annotation
set required by the chosen scheme(s) (**PlantTFDB**, **TAPscan**, and/or
**CISBP**; see variable **Annotation** in
[`list_schemes()`](list_schemes.md)). All data frames have one row per
domain hit and the following variables:

- Gene:

  Character, gene ID.

- Domain:

  Character, domain ID or name.

- qlen:

  Numeric, length of the profile HMM.

- c_evalue:

  Numeric, conditional E-value of the domain hit.

- hmm_from, hmm_to:

  Numeric, start and end coordinates of the alignment in the profile
  HMM.

Rows are kept in the same order as in the HMMER output, as this order is
used by the TAPscan classification algorithm.

## Examples

``` r
data(gsu)
seq <- gsu[1:5]
if(hmmer_is_installed()) {
    annotate_domains(seq)

    # Only domains used by CIS-BP's scheme
    annotate_domains(seq, scheme = "CISBP")
}
```
