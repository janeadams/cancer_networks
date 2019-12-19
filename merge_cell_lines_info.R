library(data.table)
setwd("/home/diana/Workspace/cnww/cancer_networks/")

d1 <- fread("data/depmap-2019q1-celllines_v2.csv")
d2 <- fread("data/cell_lines_annotations_20181226.txt")

colnames(d1)
colnames(d2)

merge <- merge(d1, d2, by.x = "DepMap_ID", by.y = "depMapID", all=T)
#### Manually changing the colnames to merge
merge <- merge[!is.na(merge$DepMap_ID), ]
merge$CCLE_Name <- ifelse(is.na(merge$CCLE_Name), merge$CCLE_ID, merge$CCLE_Name)
merge$Gender <- ifelse(is.na(merge$Gender.x), merge$Gender.y, merge$Gender.x)
merge$Source <- ifelse(is.na(merge$Source), merge$Original.Source.of.Cell.Line,
                       merge$Source)
merge <- merge[, -c("Gender.x", "Gender.y", "CCLE_ID", 
                    "Original.Source.of.Cell.Line")]
write.table(merge, file = "data/cell_info.txt", sep = "\t", row.names = F)
