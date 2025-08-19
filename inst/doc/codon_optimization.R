## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(cubar)

seq <- 'ATGCTACGA'
cf_all <- count_codons(yeast_cds)
optimal_codons <- est_optimal_codons(cf_all)
seq_opt <- codon_optimize(seq, optimal_codons)
print(seq_opt)

## -----------------------------------------------------------------------------
seq_opt <- codon_optimize(seq, cf = cf_all, method = "IDT")
print(seq_opt)

## ----eval=FALSE---------------------------------------------------------------
# seq_opt <- codon_optimize(seq, method = "CodonTransformer", organism = "Saccharomyces cerevisiae")
# print(seq_opt)
# #> 9-letter DNAString object
# #> seq: ATGTTAAGATGA

## ----eval=FALSE---------------------------------------------------------------
# seqs_opt <- codon_optimize(seq, cf = cf_all, method = "IDT", num_sequences = 10)
# print(seqs_opt)
# #> DNAStringSet object of length 6:
# #>     width seq
# #> [1]     9 ATGCTCCGT
# #> [2]     9 ATGCTGCGT
# #> [3]     9 ATGCTTCGT
# #> [4]     9 ATGCTACGT
# #> [5]     9 ATGCTCCGA
# #> [6]     9 ATGCTTCGA
# seqs_opt <- codon_optimize(seq, method = "CodonTransformer", organism = "Saccharomyces cerevisiae",
# num_sequences = 10, deterministic =FALSE, temperature = 0.4)
# print(seqs_opt)
# #> DNAStringSet object of length 4:
# #>     width seq
# #> [1]    12 ATGTTGAGATAA
# #> [2]    12 ATGTTAAGATAA
# #> [3]    12 ATGTTGAGATGA
# #> [4]    12 ATGTTGAGATAG

## ----eval=FALSE---------------------------------------------------------------
# seqs_opt <- codon_optimize(seq, cf = cf_all, method = "IDT", num_sequences = 10, spliceai = TRUE)
# print(seqs_opt)
# #>     Candidate_optimized_sequence Possible_splice_junction
# #>                           <char>                   <lgcl>
# #>  1:                    ATGCTACGC                    FALSE
# #>  2:                    ATGCTGCGA                    FALSE
# #>  3:                    ATGCTGCGT                    FALSE
# #>  4:                    ATGCTTCGC                    FALSE
# #>  5:                    ATGCTACGT                    FALSE
# #>  6:                    ATGCTTCGG                    FALSE
# #>  7:                    ATGCTCCGT                    FALSE
# #>  8:                    ATGCTTCGT                    FALSE
# #>  9:                    ATGCTCCGA                    FALSE
# #> 10:                    ATGCTTCGA                    FALSE
# seq_opt <- codon_optimize(seq, method = "CodonTransformer", organism = "Saccharomyces cerevisiae", spliceai = TRUE)
# print(seq_opt)
# #>    Candidate_optimized_sequence Possible_splice_junction
# #>                          <char>                   <lgcl>
# #> 1:                 ATGTTAAGATGA                    FALSE

