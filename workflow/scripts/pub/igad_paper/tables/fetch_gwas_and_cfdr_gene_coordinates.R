library(data.table)
library(httr)

dat <- fread(snakemake@input[[1]])

flank <- snakemake@params$window/2

# NB: Based on the sample script provided by Open Targets Genetics here: https://genetics-docs.opentargets.org/data-access/graphql-api#available-endpoints
query = "
query geneQuery($geneName: String!) {
    search(queryString: $geneName) {
        genes {
            id
            symbol
            start
            end
        }
    }
}
"

base_url <- "https://api.genetics.opentargets.org/graphql"

for(i in 1:nrow(dat)) {
  variables <- list("geneName" = dat[i, chosen_gene])
  post_body <- list(query = query, variables = variables)
  r <- POST(url = base_url, body = post_body, encode = 'json')
  res <- tryCatch(content(r)$data$search$genes[[1]], error = function(e) return(NULL))

  if(!is.null(res)) {
    print(i)
    dat[i, `:=` (start = res$start, end = res$end)]
  }
}

dat[, `:=` (start_with_flank = start - flank, end_with_flank = end + flank)]

fwrite(dat[dataset%in% c('gwas', 'both'), .(chosen_gene, CHR, start, end, start_with_flank, end_with_flank)], file = snakemake@output[['gwas']], sep = '\t')
fwrite(dat[dataset%in% c('cfdr', 'both'), .(chosen_gene, CHR, start, end, start_with_flank, end_with_flank)], file = snakemake@output[['cfdr']], sep = '\t')
