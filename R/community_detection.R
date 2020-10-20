library(igraph)
library(viridis)
library(dplyr)
library(readr)

links <- read_tsv("../data/cell_lines_corr-above0-5-0-323.net",
                    col_names =  c("from", "to", "weight"), skip = 1)
vertices <- read_tsv("../data/cell_lines_corr-above0-5-0-323-subtypes.nodes",
                     col_names = c("id", "depmapId", "primary_disease", "subtype_disease"), skip = 1)

from_to <- c(links$from, links$to)
from_to <- from_to[!duplicated(from_to)]
no_info <- vertices[!vertices$id %in% from_to, ]

g <- graph_from_data_frame(links, vertices = vertices, directed = F)

### the algorithm is deterministic
commun <- cluster_louvain(g)
save(commun, file = "../data/network_communities.RData")
#load("../data/network_communities.RData")
V(g)$commun <- commun$membership
df_vertices <- igraph::as_data_frame(g, what = "vertices")
write_tsv(df_vertices, "../data/network_vertices_with_communities.tsv")

### Network by Primary Disease
cancer_types <- names(table(V(g)$primary_disease))
node_colors <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
                 '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                 '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', 
                 '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080')
node_colors <- sample(node_colors)
names(node_colors) <- cancer_types

V(g)$color <- node_colors[V(g)$primary_disease]

l = layout_with_fr(g)
png("../figures/network_cell_type.png", 1800, 1800)
plot(g, vertex.label=NA, vertex.color=V(g)$color, 
     vertex.size= degree(g)/2 + 3,  edge.width=E(g)$weight, 
     edge.color="grey50", layout = l)
dev.off()

### Legend
png("../figures/network_legend.png", 1000, 2200)
plot.new()
legend("bottomleft", bty = "n",
       legend = names(node_colors),
       fill = node_colors, border=NA, cex = 5)

dev.off()

### Community highlighter
comm_colors = rainbow(length(commun), alpha = 0.3)
comm_colors = sample(comm_colors)
png("../figures/network_cell_type_communities_big.png", 3200, 3200)
plot(g, vertex.label=NA, vertex.color=V(g)$color, 
        vertex.size=degree(g),  edge.width=E(g)$weight, 
        edge.color="grey50", layout = l, mark.groups = groups(commun), 
        mark.border=NA, mark.col = comm_colors)
dev.off()

data.frame(color = comm_colors, community = names(groups(commun))) %>% 
        write_tsv("../data/community_colors.tsv")

### by subtype 
library(RColorBrewer)
cancer_subtype <- names(table(V(g)$subtype_disease))

node_colors <- colorRampPalette(brewer.pal(8, "Set1"))(length(cancer_subtype))
names(node_colors) <- cancer_subtype
V(g)$color <- node_colors[V(g)$subtype_disease]

png("../figures/network_cell_subtype.png", 1800, 1800)
plot(g, vertex.label=NA, vertex.color=V(g)$color, 
     vertex.size= degree(g)/2 + 3,  edge.width=E(g)$weight, 
     edge.color="grey50", layout = l)
dev.off()

png("../figures/network_subtype_legend.png", 1000, 2200)
plot.new()
legend("bottomleft", bty = "n",
       legend = names(node_colors),
       fill = node_colors, border=NA, cex = 5)
dev.off()
