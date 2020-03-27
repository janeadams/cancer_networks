library(igraph)
library(viridis)
links <- read.delim("./data/cell_lines_corr-above0-5-0-323.net",
                    col.names = c("from", "to", "weight"))
vertices <- read.delim("./data/cell_lines_corr-above0-5-0-323.nodes", 
                       col.names = c("id", "depmapId", "NULL", "NULL"),
                       colClasses = c("numeric", 'character', "NULL", "NULL"))

vertices_info <- read.delim("data/cell_info.txt")
colnames(vertices_info)
# [1] "DepMap_ID"                    "CCLE_Name"                    "Aliases"                     
# [4] "COSMIC_ID"                    "Sanger.ID"                    "Primary.Disease"             
# [7] "Subtype.Disease"              "Source"                       "Name"                        
# [10] "Pathology"                    "Site_Primary"                 "Site_Subtype1"               
# [13] "Site_Subtype2"                "Site_Subtype3"                "Histology"                   
# [16] "Hist_Subtype1"                "Hist_Subtype2"                "Hist_Subtype3"               
# [19] "Life_Stage"                   "Age"                          "Race"                        
# [22] "Geo_Loc"                      "inferred_ethnicity"           "Site_Of_Finding"             
# [25] "Disease"                      "Annotation_Source"            "Characteristics"             
# [28] "Growth.Medium"                "Supplements"                  "Freezing.Medium"             
# [31] "Doubling.Time.from.Vendor"    "Doubling.Time.Calculated.hrs" "type"                        
# [34] "type_refined"                 "PATHOLOGIST_ANNOTATION"       "mutRate"                     
# [37] "tcga_code"                    "Gender"   
vertices_info <- vertices_info[, c("DepMap_ID", "Primary.Disease", "Subtype.Disease", "Name", 
                         "Pathology", "Site_Primary", "Site_Subtype1", 
                         "Histology", "Hist_Subtype1", "Gender")]

colnames(vertices_info) <- c("depmapId", "primary_disease", "subtype_disease", "cell_name", 
                             "pathology", "primary_site", "subtype_syte", 
                             "histology", "subtype_histology", "gender")
vertices <- merge(vertices, vertices_info, by = "depmapId")
vertices <- vertices[, c("id", colnames(vertices)[-2])]

from_to <- c(links$from, links$to)
from_to <- from_to[!duplicated(from_to)]
no_info <- vertices[!vertices$id %in% from_to, ]

g <- graph_from_data_frame(links, vertices = vertices, directed = F)

node_colors <- viridis(length(table(V(g)$primary_disease)))
V(g)$color <- node_colors[as.numeric(as.factor(V(g)$primary_disease))]

l = layout_with_fr(g)
png("figures/network_cell_type.png", 1200, 1200)
plot(g, vertex.label=NA, vertex.color=V(g)$color, 
     vertex.size=degree(g),  edge.width=E(g)$weight, 
     edge.color="grey50", layout = l)
dev.off()

commun <- cluster_infomap(g, nb.trials = 10)
V(g)$commun <- commun$membership

png("figures/network_cell_type_communities.png", 1200, 1200)
plot(g, vertex.label=NA, vertex.color=V(g)$color, 
        vertex.size=degree(g),  edge.width=E(g)$weight, 
        edge.color="grey50", layout = l, mark.groups= commun, mark.border=NA)
dev.off()

df_vertices <- as_data_frame(g, what = "vertices")
write.table(df_vertices, file = "data/network_vertices_with_communities.tsv",
            row.names = F, quote = F, sep = "\t")
