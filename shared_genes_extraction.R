library(data.table)
library(dplyr)

vertices <- read.delim("./data/network_vertices_with_communities.tsv")
gene_corr <- fread("./data/gene_dependency_ranks_above_0-5.csv", sep = ",")
genes <- unlist(gene_corr[, 1])
names(genes) <- NULL
genes <- unlist(lapply(strsplit(genes, split = " "), "[[", 1))

gene_corr <- gene_corr[, -1]

cs <- unique(vertices$commun)
n <- length(genes)
nt <- n - n*0.1

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
  return(list(commun = co, genes = paste(genes[gene_ids], collapse = ",")))
})

all_genes <- bind_rows(all_genes)
write.table(all_genes, "data/network_top_genes_per_community.tsv", 
            row.names = F, quote = F, sep = "\t")
write(paste(genes, collapse = "\n"), "data/gene_list.tsv")
