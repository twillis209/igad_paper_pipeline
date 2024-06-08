library(data.table)
library(httr)

dat <- fread(snakemake@input[[1]], select = snakemake@params$gene_column_label)


# NB: Based on the sample script provided by Open Targets Genetics here: https://genetics-docs.opentargets.org/data-access/graphql-api#available-endpoints
ensembl_id_query = "
query idQuery($geneName: String!) {
    search(queryString: $geneName) {
        genes {
            id
            symbol
        }
    }
}
"

base_url <- "https://api.genetics.opentargets.org/graphql"

for(i in 1:nrow(dat)) {
  variables <- list("geneName" = dat[i, get(snakemake@params$gene_column_label)])
  post_body <- list(query = ensembl_id_query, variables = variables)
  r <- POST(url = base_url, body = post_body, encode = 'json')
  res <- tryCatch(content(r)$data$search$genes[[1]], error = function(e) return(NULL))

  if(!is.null(res)) {
    dat[i, ensemblId := res$id]
  }
}

cols <- c(snakemake@params$gene_column_label, 'ensemblId')

fwrite(dat[, ..cols], file = snakemake@output[[1]], sep = '\t')
