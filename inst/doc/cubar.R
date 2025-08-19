## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
library(cubar)
library(ggplot2)

## -----------------------------------------------------------------------------
# example data
yeast_cds

# qc
yeast_cds_qc <- check_cds(yeast_cds)
yeast_cds_qc

## -----------------------------------------------------------------------------
# convert a CDS to codon sequence
seq_to_codons(yeast_cds_qc[['YDR320W-B']])

# convert a CDS to amino acid sequence
Biostrings::translate(yeast_cds_qc[['YDR320W-B']])

## -----------------------------------------------------------------------------
# get codon frequency
yeast_cf <- count_codons(yeast_cds_qc)

## -----------------------------------------------------------------------------
yeast_cf[1:3, 1:3]

## -----------------------------------------------------------------------------
# get codon table for the standard genetic code
ctab <- get_codon_table(gcid = '1')

# plot possible codon and anticodon pairings
pairing <- ca_pairs(ctab, plot = TRUE)
plot_ca_pairs(ctab, pairing)

## -----------------------------------------------------------------------------
# example of a custom mapping
head(aa2codon)

# create a custom codon table
custom_ctab <- create_codon_table(aa2codon)
head(custom_ctab)

## -----------------------------------------------------------------------------
# get enc
enc <- get_enc(yeast_cf)
head(enc)

plot_dist <- function(x, xlab = 'values'){
    x <- stack(x)
    ggplot(x, aes(x = values)) +
        geom_histogram(bins = 40, fill = '#88CCEE') +
        labs(x = xlab, y = 'Number of genes') +
        theme_classic(base_size = 12) +
        theme(axis.text = element_text(color = 'black'))
}

plot_dist(enc, 'ENC')

## -----------------------------------------------------------------------------
# get fop
fop <- get_fop(yeast_cf)
plot_dist(fop, 'Fop')

## -----------------------------------------------------------------------------
optimal_codons <- est_optimal_codons(yeast_cf, codon_table = ctab)
head(optimal_codons[optimal == TRUE])

## -----------------------------------------------------------------------------
# estimate RSCU of highly expressed genes
yeast_exp <- as.data.table(yeast_exp)
yeast_exp <- yeast_exp[gene_id %in% rownames(yeast_cf), ]
yeast_heg <- head(yeast_exp[order(-fpkm), ], n = 500)
rscu_heg <- est_rscu(yeast_cf[yeast_heg$gene_id, ], codon_table = ctab)
head(rscu_heg) # RSCU values are shown in the column `rscu`

# calculate CAI of all genes
# note: CAI values are usually calculated based RSCU of highly expressed genes.
cai <- get_cai(yeast_cf, rscu = rscu_heg)
plot_dist(cai, xlab = 'CAI')

## -----------------------------------------------------------------------------
# get tRNA gene copy number from GtRNADB
trna_gcn <- extract_trna_gcn(yeast_trna)

# calculate tRNA weight for each codon
trna_w <- est_trna_weight(trna_level = trna_gcn, codon_table = ctab)

# get tAI
tai <- get_tai(yeast_cf, trna_w = trna_w)
plot_dist(tai, 'tAI')

## -----------------------------------------------------------------------------
# path_gtrnadb <- 'http://gtrnadb.ucsc.edu/genomes/eukaryota/Scere3/sacCer3-mature-tRNAs.fa'
# yeast_trna <- Biostrings::readRNAStringSet(path_gtrnadb)

## -----------------------------------------------------------------------------
pairs(cbind(CAI = cai, ENC = enc, Fop = fop, TAI = tai),
      cex = 0.5, col = alpha('black', 0.2))

## -----------------------------------------------------------------------------
# get lowly expressed genes
yeast_leg <- head(yeast_exp[order(fpkm), ], n = 500)
yeast_leg <- yeast_leg[gene_id %in% rownames(yeast_cf), ]

# differetial usage test
du_test <- codon_diff(yeast_cds_qc[yeast_heg$gene_id], yeast_cds_qc[yeast_leg$gene_id])
du_test <- du_test[amino_acid != '*', ]

## -----------------------------------------------------------------------------
du_test$codon <- factor(du_test$codon, levels = du_test[order(-global_or), codon])

ggplot(du_test, aes(x = codon, y = log2(global_or))) +
    geom_col() +
    labs(x = NULL, y = 'log2(OR)') +
    theme_classic(base_size = 12) +
    theme(axis.text = element_text(color = 'black'),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## -----------------------------------------------------------------------------
du_test2 <- du_test[!amino_acid %in% c('Met', 'Trp'), ]
du_test2$codon <- factor(du_test2$codon, levels = du_test2[order(-fam_or), codon])

ggplot(du_test2, aes(x = codon, y = log2(fam_or))) +
    geom_col() +
    labs(x = NULL, y = 'log2(OR)') +
    facet_grid(cols = vars(amino_acid), space = 'free', scales = 'free', drop = TRUE) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## -----------------------------------------------------------------------------
# optimize codon usage to the optimal codon of each amino acid
opc_aa <- est_optimal_codons(yeast_cf, codon_table = ctab, level = 'amino_acid')
seq_optimized <- codon_optimize(yeast_cds_qc[['YFR025C']], optimal_codons, level = 'amino_acid')

# CAI before and after optimization
plot_dist(cai, 'CAI') + 
    geom_vline(xintercept = cai['YFR025C'], linetype = 'dashed', color = 'red') +  # before
    geom_vline(xintercept = get_cai(count_codons(seq_optimized), rscu_heg),
               linetype = 'dashed', color = 'black')  # after

## -----------------------------------------------------------------------------
swa <- slide_apply(yeast_cds_qc[['YHR099W']], .f = get_cai,
                   step = 30, before = 20, after = 20, rscu = rscu_heg)

# plot the results
slide_plot(swa, 'CAI')

