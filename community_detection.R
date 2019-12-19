library(igraph)
library(data.table)

setwd("/home/diana/Workspace/cnww/cancer_networks/")
links <- fread("data/edges_p3.txt")
colnames(links) <- c("from", "to", "weight")
vertices <- fread("data/cell_info.txt")
colnames(vertices)
# [1] "DepMap_ID"                    "CCLE_Name"                   
# [3] "Aliases"                      "COSMIC_ID"                   
# [5] "Sanger ID"                    "Primary Disease"             
# [7] "Subtype Disease"              "Source"                      
# [9] "Name"                         "Pathology"                   
# [11] "Site_Primary"                 "Site_Subtype1"               
# [13] "Site_Subtype2"                "Site_Subtype3"               
# [15] "Histology"                    "Hist_Subtype1"               
# [17] "Hist_Subtype2"                "Hist_Subtype3"               
# [19] "Life_Stage"                   "Age"                         
# [21] "Race"                         "Geo_Loc"                     
# [23] "inferred_ethnicity"           "Site_Of_Finding"             
# [25] "Disease"                      "Annotation_Source"           
# [27] "Characteristics"              "Growth.Medium"               
# [29] "Supplements"                  "Freezing.Medium"             
# [31] "Doubling.Time.from.Vendor"    "Doubling.Time.Calculated.hrs"
# [33] "type"                         "type_refined"                
# [35] "PATHOLOGIST_ANNOTATION"       "mutRate"                     
# [37] "tcga_code"                    "Gender
vertices <- vertices[, c("DepMap_ID", "Primary Disease", "Subtype Disease", "Name", 
                         "Pathology", "Site_Primary", "Site_Subtype1", "Site_Subtype2",
                         "Site_Subtype3", "Histology", "Hist_Subtype1", "Hist_Subtype2",
                         "Hist_Subtype3")]
from_to <- c(links$from, links$to)
from_to <- from_to[!duplicated(from_to)]
colnames(vertices) <- tolower(colnames(vertices))
colnames(vertices)[4] <- "cell_line_name"
colnames(vertices)[1] <- "name"
vertices <- vertices[vertices$name %in% from_to, ]
noinfo <- c("ACH-001366","ACH-001647","ACH-001838")
links <- links[!links$from %in% noinfo, ]
links <- links[!links$to %in% noinfo, ]
g <- graph_from_data_frame(links, vertices = vertices, directed = F)
commun <- cluster_louvain(g)
V(g)$louvain <- membership(commun)
g_df <- as_data_frame(g, what = "vertices")
write.table(g_df[, c("name", "louvain")], file = "data/louvain_clusters_3.txt",
            row.names = F, quote = F, sep = "\t")
