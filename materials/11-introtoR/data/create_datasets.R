# John Blischak
# 2013-09-16

fetch_biomart_data <- function(genes = NULL) {
  # Fetches Ensembl Biomart data for translated human genes.
  #
  # Args
  #  genes: character vector of ensembl gene IDs.
  #         If NULL, all genes are fetched.
  #
  # Results
  #  A data.frame where each row is the Ensembl transcipt with the
  #  longest coding sequence for a given Ensembl gene ID.
  #
  require(biomaRt, quietly = TRUE)
  require(plyr, quietly = TRUE)
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#   listFilters(ensembl)
#   listAttributes(ensembl, "feature_page")
#   attributePages(ensembl)
  stopifnot(is.null(genes) || is.vector(genes, mode = "character"))
  atts <- c("ensembl_transcript_id", "ensembl_gene_id", "chromosome_name", 
            "transcript_start", "transcript_end",
            "external_gene_id", "percentage_gc_content", "gene_biotype",
            "source", "name_1006")
  if (is.null(genes)) {
    gene_info <- getBM(attributes = atts, mart = ensembl)
  } else {
    gene_info <- getBM(attributes = atts, filters = "ensembl_gene_id",
                       values = genes, mart = ensembl)
  }
  gene_coding_len <- getBM(attributes = c("ensembl_transcript_id", "cds_length"),
                           filters = "ensembl_transcript_id",
                           values = gene_info$ensembl_transcript_id,
                           mart = ensembl)
  gene_coding_len <- na.omit(gene_coding_len)
  gene_merged <- merge(gene_info, gene_coding_len, by = "ensembl_transcript_id")
  gene_final <- ddply(gene_merged, .(ensembl_gene_id), 
                      function(df) df[which.max(df$cds_length), ])
  return(gene_final)
}  

# Testing:
# g <- c("ENSG00000198888", "ENSG00000198763", "ENSG00000198804",
#        "ENSG00000198712", "ENSG00000228253", "ENSG00000198899")
# test_data <- fetch_biomart_data(g)

# Create dataset for bootcamp
gene_data <- fetch_biomart_data()
write.table(gene_data, file = "ensembl_human_genes.txt", sep = "\t",
            quote = FALSE, row.names = FALSE)
