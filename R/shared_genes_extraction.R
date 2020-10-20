library(data.table)
library(dplyr)
library(readr)

vertices <- read.delim("../data/network_vertices_with_communities.tsv")
gene_corr <- fread("../data/gene_dependency_ranks_above_0-5.csv", sep = ",")
genes <- unlist(gene_corr[, 1])
names(genes) <- NULL
genes <- unlist(lapply(strsplit(genes, split = " "), "[[", 1))

gene_corr <- gene_corr[, -1]

cs <- unique(vertices$commun)
n <- length(genes)
nt <- n - n*0.05

### Genes in the top 5 per cent
all_genes <- lapply(cs, FUN = function(co){
  
  cell_lines <- vertices %>% filter(commun == co)
  cell_lines_m <- gene_corr %>% select(cell_lines$depmapId)
  cell_lines_m <- data.matrix(cell_lines_m)
  
  gene_ids <- lapply(1:nrow(cell_lines_m), FUN= function(i){
    if(all(!is.na(cell_lines_m[i,]))) {
      if(all(cell_lines_m[i,] > nt)) {
        return(i)
      }
    }
  })
  
  gene_ids <- unlist(gene_ids)
  cat("Community ", co, " has ", length(gene_ids), " common genes.\n")
  gene_ids = data.frame(gene = genes[gene_ids])
  gene_ids$community <- co
  return(gene_ids)
})

all_genes <- bind_rows(all_genes) 
in_all <- all_genes %>% group_by(gene) %>% tally(sort = T) %>% filter(n == 19)
in_one <- all_genes %>% group_by(gene) %>% tally(sort = T) %>% filter(n == 1)

in_all %>% select(gene) %>% write_tsv("../data/genes_in_all_communities.tsv")

in_one <- all_genes %>% semi_join(in_one, by = "gene") 

community_colors <- read_tsv("../data/community_colors.tsv")

in_one %>% group_by(community) %>% tally() %>% left_join(community_colors, by = 'community') %>% 
  write_tsv("../data/communities_with_exclusive_genes.tsv")
