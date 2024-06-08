library(data.table)

non_imd_files <- snakemake@input[names(snakemake@input) %like% 'with|sans']
imd_files <- snakemake@input$imd_gif_files

dats <- list()

for(i in seq_along(non_imd_files)) {
  dat <- fread(non_imd_files[[i]])

  if(basename(non_imd_files[[i]]) %like% 'with_mhc') {
    dat[, variant_set := 'with_mhc']
  } else if(basename(non_imd_files[[i]]) %like% 'sans_mhc')  {
    dat[, variant_set := 'sans_mhc']
  } else {
    stop('Cannot determine if with/sans MHC')
  }

  dat[, dataset := gsub('_', '-', gsub('(_?(with|sans)_mhc_?)|_gif|\\.tsv', '', names(non_imd_files[i])))]

  dats[[i]] <- dat
}

non_imd_dat <- rbindlist(dats)

dats <- list()

for(i in seq_along(imd_files)) {
  dat <- fread(imd_files[[i]])

  if(basename(imd_files[[i]]) %like% 'with_mhc') {
    dat[, variant_set := 'with_mhc']
  } else if(basename(imd_files[[i]]) %like% 'sans_mhc')  {
    dat[, variant_set := 'sans_mhc']
  } else {
    stop('Cannot determine if with/sans MHC')
  }

  dat[, dataset := gsub('_', '-', gsub('(_?(with|sans)_mhc_?)|_gif|\\.tsv', '', basename(imd_files[i])))]

  dats[[i]] <- dat
}

imd_dat <- rbindlist(dats)

fwrite(rbindlist(list(non_imd_dat, imd_dat)), sep = '\t', file = snakemake@output[[1]])
