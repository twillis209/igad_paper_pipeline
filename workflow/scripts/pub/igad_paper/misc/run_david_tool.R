library(data.table)
library(magrittr)
library(httr)

dat <- fread(snakemake@input[[1]], select = 'ensemblId')

base_url <- "https://davidbioinformatics.nih.gov/api.jsp"

# annot tags
go_term_tags <- c("GOTERM_BP_1",
"GOTERM_BP_2",
"GOTERM_BP_3",
"GOTERM_BP_4",
"GOTERM_BP_5",
"GOTERM_BP_ALL",
"GOTERM_BP_FAT",
"GOTERM_CC_1",
"GOTERM_CC_2",
"GOTERM_CC_3",
"GOTERM_CC_4",
"GOTERM_CC_5",
"GOTERM_CC_ALL",
"GOTERM_CC_FAT",
"GOTERM_MF_1",
"GOTERM_MF_2",
"GOTERM_MF_3",
"GOTERM_MF_4",
"GOTERM_MF_5",
"GOTERM_MF_ALL",
"GOTERM_MF_FAT")

pathway_tags <- c("BBID",
"BIOCARTA",
"EC_NUMBER",
"KEGG_COMPOUND",
"KEGG_PATHWAY",
"KEGG_REACTION")

functional_categories_tags <- c("CGAP_EST_QUARTILE",
"CGAP_EST_RANK",
"COG_ONTOLOGY",
"PIR_SEQ_FEATURE",
"SP_COMMENT_TYPE",
"SP_PIR_KEYWORDS",
"UP_SEQ_FEATURE")

params <- list(
  'type' = 'ENSEMBL_GENE_ID',
  'ids' = paste(dat$ensemblId, collapse = ','),
  'tool' = 'term2term',
  'annot' = paste(c(go_term_tags, pathway_tags, functional_categories_tags), collapse = ',')
  )

#r <- POST(url = base_url, body = params, encode = 'form')
#res <- tryCatch(content(r)$data$search$genes[[1]], error = function(e) return(NULL))

writeLines(modify_url(base_url, query = params), snakemake@output[[1]])
