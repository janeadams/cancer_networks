library(topGO)
library(dplyr)

gene_sets <- read.delim("data/network_top_genes_per_community_5_per.tsv", sep = "\t", 
                        stringsAsFactors = F)
all_genes <- read.delim("data/gene_list.tsv", stringsAsFactors = F, header = F)
all_genes <- as.character(all_genes[, 1])

all_enrichments <- lapply(1:nrow(gene_sets), function(n){
  if(nchar(gene_sets[n, 2]) > 0) {
    geneList <- unlist(strsplit(gene_sets[n, 2], split = ","))
    geneList <- setNames(factor(as.integer(all_genes %in% geneList)), all_genes)
    
    GOdata <- new("topGOdata",
                  description = "GO analysis of genes in community",
                  ontology = "BP",
                  allGenes = geneList,
                  nodeSize = 10,
                  annot = annFUN.org,
                  ID = "alias", 
                  mapping = "org.Hs.eg")
    
    resultFisher <- runTest(GOdata, "elim","fisher")
    ss <- score(resultFisher)
    cat("Community ", gene_sets[n, 1], " has ", sum(ss < 0.001), " significant processes")
    resultsTable <- GenTable(GOdata, classicFisher = resultFisher, numChar = 200,
                             topNodes = sum(ss < 0.001))
    resultsTable$commun <- gene_sets[n, 1]
    return(resultsTable)
  }
})
all_enrichments <- bind_rows(all_enrichments)
bps_in_all_comm <- all_enrichments %>% 
  group_by(GO.ID, Term) %>% tally() %>% filter(n == 23)

bps_not_in_all_comm <- all_enrichments %>% 
  group_by(GO.ID, Term) %>% tally() %>% filter(n != 23)

write.table(all_enrichments, file = "data/all_enrichment_scores_bp_above_10_genes_5_per.tsv",
            quote = F, row.names = F, sep = "\t")
write.table(bps_in_all_comm, file = "data/bps_in_all_comm_bp_above_10_genes_5_per.tsv",
            quote = F, row.names = F, sep = "\t")
write.table(bps_not_in_all_comm, file = "data/bps_not_in_all_comm_bp_above_10_genes_5_per.tsv",
            quote = F, row.names = F, sep = "\t")
