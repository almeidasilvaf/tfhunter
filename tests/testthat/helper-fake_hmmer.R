# Create a data frame of hmmsearch hits, as returned by .read_hmmsearch()
fake_hits <- function(
        target_name, query_name, full_evalue = 1e-10, full_score = 50,
        c_evalue = 1e-10, i_evalue = 1e-10, dom_score = 50, qlen = 60,
        hmm_from = 1, hmm_to = 60, ali_from = 1, ali_to = 60
) {
    n <- length(target_name)
    hits <- data.frame(
        target_name = target_name, target_accession = "-", tlen = 100,
        query_name = query_name, query_accession = "-", qlen = qlen,
        full_evalue = full_evalue, full_score = full_score, full_bias = 0,
        dom_n = 1, dom_of = 1, c_evalue = c_evalue, i_evalue = i_evalue,
        dom_score = dom_score, dom_bias = 0, hmm_from = hmm_from,
        hmm_to = hmm_to, ali_from = ali_from, ali_to = ali_to,
        env_from = ali_from, env_to = ali_to, acc = 0.9,
        description = rep("-", n)
    )
    return(hits)
}

# Create a fake 'hmmsearch' program that writes a single domain hit to the
# file passed to --domtblout and exits with the given status. Returns the
# directory containing the program.
fake_hmmsearch <- function(status = 0) {
    dir <- tempfile(pattern = "fake_hmmer")
    dir.create(dir)
    script <- c(
        "#!/bin/sh",
        "while [ $# -gt 0 ]; do",
        "  if [ \"$1\" = \"--domtblout\" ]; then out=\"$2\"; fi",
        "  shift",
        "done",
        "echo '# target name' > \"$out\"",
        paste(
            "echo 'g1 - 100 PF00046 PF00046.1 57 1e-20 70.0 0.1 1 1 1e-22",
            "1e-20 69.0 0.1 1 57 10 66 9 67 0.98 a description' >> \"$out\""
        ),
        paste("exit", status)
    )
    path <- file.path(dir, "hmmsearch")
    writeLines(script, path)
    Sys.chmod(path, mode = "0755")
    return(dir)
}
