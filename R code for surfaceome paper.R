### ---- preliminaries ----
library(Seurat)
library(tidyverse)
library(correlation)
library(psych)
library(descr)
library(DescTools)
library(writexl)

## load custom functions
source('./funcs.R')

### ---- get TMS data ----
## raw data: https://figshare.com/articles/dataset/Processed_files_to_use_with_scanpy_/8273102/3?file=23936684 (version 3)
## github: https://github.com/czbiohub/tabula-muris-senis

## select droplet data 
sel_dataset <- 'droplet'

## raw data directory
path_data <- file.path('data', 'tabula_muris', sel_dataset, 'version3')
#fid <- file.path(path_data, 'tabula-muris-senis-droplet-processed-official-annotations.h5ad')
fid <- paste0(path_data, '/', list.files(path_data))   

## list of genes of interest
genes_interest <- c('PLXNA1', 'PLXNA3', 'PTK7', 'CYB5R1', 'CDKN2A')


### ---- data info ----
## get meta info on genes (ID, mean, dispersion, variablity)
info_genes <- h5read(file = fid, name = 'var')   
num_genes <- nrow(info_genes)

## capitalise gene names
info_genes$index <- toupper(info_genes$index)           

## index of senescence gene
idx_cdkn2a <- which(info_genes$index == 'CDKN2A')

## indices of surface marker genes of interest
idx_markers <- which(info_genes$index %in% genes_interest)

## get total number of senescent cells (all tissues)
num_sen_total <- info_genes %>% filter(index == 'CDKN2A') %>% select(n_cells) %>% unlist()


### ---- load TMS data ----
## get expression data (normalised)
dd_data <- h5read(file = fid, name = 'X/data')
dd_indices <- h5read(file = fid, name = 'X/indices') + 1   #convert from python to R indexing
dd_indptr <- h5read(file = fid, name = 'X/indptr') + 1     #convert from python to R indexing

## get meta info on cells
info_cells <- h5read(file = fid, name = 'obs') 
info_cells$tissue_name <- NA            #pre-allocation
info_cells$cell_ontology_name <- NA     #pre-allocation
num_cells <- nrow(info_cells)           #total number of cells

## get meta info on various donor categories
info_meta <- h5read(file = fid, name = 'uns')  

## get names of tissue types 
tissue_names <- info_meta$tissue_categories
for(t in seq_along(tissue_names)){
  idx_class <- which(info_cells$tissue+1 == t)
  info_cells$tissue_name[idx_class] <- tissue_names[t]
}

## get names of cell types 
cell_names <- info_meta$cell_ontology_class_categories
for(c in seq_along(cell_names)){
  idx_class <- which(info_cells$cell_ontology_class+1 == c)
  info_cells$cell_ontology_name[idx_class] <- cell_names[c]
}

## get meaningful age 
curr_age_cats <- info_meta$age_categories
tmp <- info_cells$age+1
for(a in seq_along(curr_age_cats)){
  idx_age <- which(tmp == a)
  info_cells$age[idx_age] <- curr_age_cats[a]
}
info_cells$age <- as.integer(gsub(pattern = 'm', replacement = '', info_cells$age))

## check for discrepancies
tmp <- h5read(file = fid, name = 'var')   
tmp$index <- toupper(tmp$index)
stopifnot(all.equal(info_genes$index, tmp$index))

## find surface marker genes (genes of interest)
idx_cols <- na.omit(match(genes_interest, info_genes$index))

## reconstruct (sparse) expression matrix
xx <- matrix(0, nrow = num_cells, ncol = length(idx_cols))
info_cells$cdkn2a <- NA
chunk_size <- 5e4                           #note: reduce if not enough RAM
num_chunks <- ceiling(num_cells/chunk_size)

for(j in seq(num_chunks)){                 
  cat('processing chunk ', j, '/', num_chunks, '... \n')
  idx_chunk_start <- 1 + (j-1)*chunk_size
  idx_chunk_end <- min(j*chunk_size, num_cells)
  chunk_length <- 1 + idx_chunk_end - idx_chunk_start
  xx_chunk <- matrix(0, nrow = chunk_length, ncol = num_genes)
  
  ## construct gene expression matrix
  for(i in seq(idx_chunk_start, idx_chunk_end)){
    curr_values <- dd_data[dd_indptr[i]:(dd_indptr[i+1]-1)]    #excluding right end of interval, i.e. using [a,b)
    curr_genes <- dd_indices[dd_indptr[i]:(dd_indptr[i+1]-1)]
    xx_chunk[i%%chunk_length, curr_genes] <- curr_values
    
    ## select only surface marker genes (genes of interest)
    xx[idx_chunk_start:idx_chunk_end, ] <- xx_chunk[, idx_cols]
  }
}

## clean-up
rm(xx_chunk); rm(dd_data); rm(dd_indices)

## name columns of gene expression matrix xx
colnames(xx) <- genes_interest

## create composite variable: tissue-cell-type
info_cells$tissue_cell_name <- character(num_cells)
for(j in seq_along(tissue_names)){
  curr_idx <- which(info_cells$tissue_name == tissue_names[j])
  info_cells$tissue_cell_name[curr_idx] <- paste0(info_cells$tissue_name[curr_idx],'-',
                                                  info_cells$cell_ontology_name[curr_idx])
}

## get list of all tissue-cell-types
tissue_cell_names <- unique(info_cells$tissue_cell_name)

## create binary age variable: young = <=3m, old = >=18m
info_cells$age_young <- info_cells$age < 18

## gene expression matrix as dataframe (tibble) and add relevant variables
xx <- as_tibble(xx) 
xx$tissue_cell_name <- info_cells$tissue_cell_name
xx$age_young <- info_cells$age_young


### ---- intermediary output: save processed data ----
fid_out <- file.path('data', paste0('processed_data_', sel_dataset, '_v3.rds'))
if(!file.exists(fid_out)){
  saveRDS(list('expr_mat' = xx, 'info' = info_cells), file = fid_out)
}else{
  warning('File already exists!')
}


### ---- data exploration ----
cat('\n number of sen+ cells per cell type: \n')
table(info_cells$tissue_cell_name[xx$CDKN2A > 0])

cat('\n number of sen+ cells per tissue type: \n')
table(info_cells$tissue_name[xx$CDKN2A > 0])

cat('\n fraction of sen+ cells per tissue type: \n')
table(info_cells$tissue_name[xx$CDKN2A > 0]) %>% 
  `/`(table(info_cells$tissue_name)) %>% round(., digits = 3)

## initialise output matrices 
out_all <- out_young <- out_old <- list('count_mat' = matrix(NA, length(tissue_cell_names), length(genes_interest)),
                                        'prop_mat' = matrix(NA, length(tissue_cell_names), length(genes_interest)))

## loop over all tissue-cell-types: count / proportion of (+) cells 
for(j in seq_along(tissue_cell_names)){
  out_all$count_mat[j, ] <- xx %>% filter(tissue_cell_name == tissue_cell_names[j]) %>% 
    select(all_of(genes_interest)) %>% `>`(0) %>% colSums(.)
  out_all$prop_mat[j, ] <- xx %>% filter(tissue_cell_name == tissue_cell_names[j]) %>% 
    select(all_of(genes_interest)) %>% `>`(0) %>% apply(., 2, mean)
  
  out_young$count_mat[j, ] <- xx %>% filter(tissue_cell_name == tissue_cell_names[j]) %>% 
    filter(age_young) %>% select(all_of(genes_interest)) %>% `>`(0) %>% colSums(.)
  out_young$prop_mat[j, ] <- xx %>% filter(tissue_cell_name == tissue_cell_names[j]) %>% 
    filter(age_young) %>%select(all_of(genes_interest)) %>% `>`(0) %>% apply(., 2, mean)
  
  out_old$count_mat[j, ] <- xx %>% filter(tissue_cell_name == tissue_cell_names[j]) %>% 
    filter(!age_young) %>%select(all_of(genes_interest)) %>% `>`(0) %>% colSums(.)
  out_old$prop_mat[j, ] <- xx %>% filter(tissue_cell_name == tissue_cell_names[j]) %>% 
    filter(!age_young) %>%select(all_of(genes_interest)) %>% `>`(0) %>% apply(., 2, mean)
}

## assign column names: genes of interest
colnames(out_all$count_mat) <- colnames(out_all$prop_mat) <- genes_interest
colnames(out_young$count_mat) <- colnames(out_young$prop_mat) <- genes_interest
colnames(out_old$count_mat) <- colnames(out_old$prop_mat) <- genes_interest

## assign row names: tissue-cell-types
rownames(out_all$count_mat) <- rownames(out_all$prop_mat) <- tissue_cell_names
rownames(out_young$count_mat) <- rownames(out_young$prop_mat) <- tissue_cell_names
rownames(out_old$count_mat) <- rownames(out_old$prop_mat) <- tissue_cell_names

## ---- LORs ----
## initialise output matrices
out_odds <- matrix(NA, length(tissue_cell_names), length(genes_interest))
p_odds <- matrix(NA, length(tissue_cell_names), length(genes_interest))

## loop over all tissue-cell-types
for(j in seq_along(tissue_cell_names)){
  ## find indices of current tissue-cell-type
  curr_idx <- which(info_cells$tissue_cell_name == tissue_cell_names[j])
  
  ## compute odds ratios for current tissue-cell-type
  tmp <- xx[curr_idx, ] %>% select(all_of(genes_interest)) %>% 
    my_getLogOR(., y = as.numeric(xx$CDKN2A[curr_idx] > 0))
  out_odds[j, ] <- tmp$logor
  p_odds[j, ] <- tmp$pval
  rm(tmp)
}

## assign column and row names
colnames(out_odds) <- genes_interest 
colnames(p_odds) <- genes_interest
rownames(out_odds) <- tissue_cell_names 
rownames(p_odds) <- tissue_cell_names


### ---- save output to file ----
data_out <- list('gene_mat' = xx, 
                 'info_cells' = info_cells)

results_out <- list('count_all' = out_all$count_mat,
                    'count_young' = out_young$count_mat,
                    'count_old' = out_old$count_mat,
                    'prop_all' = out_all$prop_mat,
                    'prop_young' = out_young$prop_mat,
                    'prop_old' = out_old$prop_mat, 
                    'logodds_ratios' = out_odds[, 1:4],      
                    'logodds_pval' = p_odds[, 1:4])


## export results
if(!dir.exists('results')){
  dir.create('results')
}
fid_data_out <- file.path('results', 'data.rds')
fid_results_out <- file.path('results', 'results.rds')

saveRDS(data_out, file = fid_data_out)
saveRDS(results_out, file = fid_results_out)

## export as excel spreadsheet
library(writexl)
r2 <- lapply(results_out, function(x) as.data.frame(x))
r2 <- lapply(r2, function(x) bind_cols(tissue_cell_names, x))
write_xlsx(r2, path = file.path('results', 'results.xlsx'))

rm(list = ls())


### ---- get Brain data ----
## raw data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129788

dir='GSE129788_RAW/' 
samples=list.files( dir ,pattern = 'gz')
samples 
library(data.table)
ctList = lapply(samples,function(pro){ 
  # pro=samples[1] 
  print(pro)
  ct=fread(file.path( dir ,pro),data.table = F)
  ct[1:4,1:4]
  rownames(ct)=ct[,1]
  colnames(ct) = paste(gsub('.txt.gz','',pro),
                       colnames(ct) ,sep = '_')
  ct=ct[,-1] 
  return(ct)
})

lapply(ctList, dim)
tmp =table(unlist(lapply(ctList, rownames)))
cg = names(tmp)[tmp==length(samples)]
bigct = do.call(cbind,
                lapply(ctList,function(ct){ 
                  ct = ct[cg,] 
                  return(ct)
                }))
sce.all=CreateSeuratObject(counts = bigct, 
                           min.cells = 5,
                           min.features = 300)
sce.all
as.data.frame(sce.all@assays$RNA$counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all@meta.data$orig.ident)

# Load the supplementary metadata
meta_data_path <- "GSE129788_Supplementary_meta_data_Cell_Types_Etc.xlsx"
supplementary_metadata <- read_excel(meta_data_path)

# Inspect the supplementary metadata
head(supplementary_metadata)
supplementary_metadata$cell_id <- gsub(".*_", "", supplementary_metadata$sample_id)

# Merge the supplementary metadata with the Seurat object's metadata
merged_metadata <- cbind(sce.all@meta.data, supplementary_metadata)
sce.all@meta.data <- merged_metadata

info_cells <-sce.all@meta.data

if (any(duplicated(colnames(info_cells)))) {
  colnames(info_cells) <- make.unique(colnames(info_cells))
}

if ("animal_type" %in% colnames(info_cells)) {
  print(str(info_cells))
  if (is.factor(info_cells$animal_type)) {
    info_cells <- info_cells %>%
      mutate(animal_type = as.character(animal_type))  
  }
  print(str(info_cells))
  info_cells <- info_cells %>%
    mutate(age_young = animal_type == "young")
  print(head(info_cells))
  tissue_cell_names <- unique(info_cells$tissue_cell_name)
  genes_interest <- c('Plxna1', 'Plxna3', 'Ptk7', 'Cyb5r1', 'Cdkn2a')  
  out_all <- out_young <- out_old <- list(
    'count_mat' = matrix(NA, length(tissue_cell_names), length(genes_interest)),
    'prop_mat' = matrix(NA, length(tissue_cell_names), length(genes_interest))
  )
  
  for (j in seq_along(tissue_cell_names)) {
    out_all$count_mat[j, ] <- info_cells %>% filter(tissue_cell_name == tissue_cell_names[j]) %>% 
      select(all_of(genes_interest)) %>% `>`(0) %>% colSums(.)
    
    out_all$prop_mat[j, ] <- info_cells %>% filter(tissue_cell_name == tissue_cell_names[j]) %>% 
      select(all_of(genes_interest)) %>% `>`(0) %>% apply(., 2, mean)
    
    out_young$count_mat[j, ] <- info_cells %>% filter(tissue_cell_name == tissue_cell_names[j]) %>% 
      filter(age_young) %>% select(all_of(genes_interest)) %>% `>`(0) %>% colSums(.)
    
    out_young$prop_mat[j, ] <- info_cells %>% filter(tissue_cell_name == tissue_cell_names[j]) %>% 
      filter(age_young) %>% select(all_of(genes_interest)) %>% `>`(0) %>% apply(., 2, mean)
    
    out_old$count_mat[j, ] <- info_cells %>% filter(tissue_cell_name == tissue_cell_names[j]) %>% 
      filter(!age_young) %>% select(all_of(genes_interest)) %>% `>`(0) %>% colSums(.)
    
    out_old$prop_mat[j, ] <- info_cells %>% filter(tissue_cell_name == tissue_cell_names[j]) %>% 
      filter(!age_young) %>% select(all_of(genes_interest)) %>% `>`(0) %>% apply(., 2, mean)
  }
  
  ## Assign column names: genes of interest
  colnames(out_all$count_mat) <- colnames(out_all$prop_mat) <- genes_interest
  colnames(out_young$count_mat) <- colnames(out_young$prop_mat) <- genes_interest
  colnames(out_old$count_mat) <- colnames(out_old$prop_mat) <- genes_interest
  
  ## Assign row names: tissue-cell-types
  rownames(out_all$count_mat) <- rownames(out_all$prop_mat) <- tissue_cell_names
  rownames(out_young$count_mat) <- rownames(out_young$prop_mat) <- tissue_cell_names
  rownames(out_old$count_mat) <- rownames(out_old$prop_mat) <- tissue_cell_names
  
  ## ---- LORs ----
  ## initialise output matrices
  out_odds <- matrix(NA, length(tissue_cell_names), length(genes_interest))
  p_odds <- matrix(NA, length(tissue_cell_names), length(genes_interest))
  
  ## loop over all tissue-cell-types
  for (j in seq_along(tissue_cell_names)) {
    ## find indices of current tissue-cell-type
    curr_idx <- which(info_cells$tissue_cell_name == tissue_cell_names[j])
    
    ## compute odds ratios for current tissue-cell-type
    tmp <- info_cells[curr_idx, ] %>% select(all_of(genes_interest)) %>% 
      my_getLogOR(., y = as.numeric(info_cells$Cdkn2a[curr_idx] > 0))
    out_odds[j, ] <- tmp$logor
    p_odds[j, ] <- tmp$pval
    rm(tmp)

  }
  
  ## assign column and row names
  colnames(out_odds)  <- genes_interest
  colnames(p_odds)  <- genes_interest
  rownames(out_odds) <- tissue_cell_names
  rownames(p_odds)  <- tissue_cell_names
  
  results_out <- list(
    'count_all' = out_all$count_mat,
    'count_young' = out_young$count_mat,
    'count_old' = out_old$count_mat,
    'prop_all' = out_all$prop_mat,
    'prop_young' = out_young$prop_mat,
    'prop_old' = out_old$prop_mat, 
    'logodds_ratios' = out_odds[, 1:4],       
    'logodds_pval' = p_odds[, 1:4],
  )
  
  # 保存结果到文件
  saveRDS(results_out, file = "results_out.rds")
  
} else {
  stop("Please ensure that the metadata contains an 'animal_type' column.")
}


### ---------- funcs.r ---------
### compute log-odds ratio
my_compLogOR <- function(x, y, labels = c('A', 'B')){
  
  mat <- table(x, y, dnn = labels)
  or <- prod(diag(mat)) / prod(diag(my_Rotate90(mat)))
  
  return(log(or))
}


### get log-odds ratios for matrix xx vs a vector y
my_getLogOR <- function(xx, y, xx_names = colnames(xx), y_name = 'CDKN2A'){
  
  ## initialise output
  xx <- as.matrix(xx)
  num_cols <- ncol(xx)
  logor <- pval <- numeric(num_cols)
  
  ## check for constant columns
  idx_const <- which(apply(xx, 2, sd) == 0)
  if(length(idx_const) > 0){ warning('constant columns detected!') }
  
  ## loop over columns in xx and compute odds ratio
  for(i in seq(num_cols)){
    logor[i] <- my_compLogOR(y > 0, xx[, i] > 0, 
                             labels = c(y_name, xx_names[i]))
    if(length(unique(xx[, i])) < 2 | length(unique(y)) < 2){
      pval[i] <- NA
    }else{
      pval[i] <- fisher.test(y > 0, xx[, i] > 0)$p.value
    }
  }
  
  ## set infinities to NA (resulting from a zero in the contingency table)
  logor[abs(logor) == Inf] = NA
  if(length(idx_const) > 0){ logor[idx_const] <- NA }
  
  ## output
  out <- list('logor' = logor, 'pval' = pval)
  
  return(out)
}


### rotate a matrix clock-wise by 90 degrees
my_Rotate90 <- function(xx){
  #' Rotate a matrix clock-wise by 90 degrees
  #'
  #' @param xx matrix
  #' @return \item{rr}{ rotated matrix }
  #' @export
  
  rr <- t(apply(xx, 2, rev))
  return(rr)
}

  