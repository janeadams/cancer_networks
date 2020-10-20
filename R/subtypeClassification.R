library(dplyr)
library(readr)

vertices <- read_tsv("./data/cell_lines_corr-above0-5-0-323.nodes")
vertices <- vertices %>% select(Node_ID, DepMap_ID)
colnames(vertices) <- c("id", "depmap_id")

vertices_info <- read.delim("data/cell_info.txt", stringsAsFactors = F)

vertices_info <- vertices_info[, c("DepMap_ID", "Primary.Disease", "Subtype.Disease", "Name", 
                                   "Pathology", "Site_Primary", "Site_Subtype1", 
                                   "Histology", "Hist_Subtype1", "Gender")]

colnames(vertices_info) <- c("depmap_id", "primary_disease", "subtype_disease", "cell_name", 
                             "pathology", "primary_site", "subtype_syte", 
                             "histology", "subtype_histology", "gender")
vertices <- merge(vertices, vertices_info, by = "depmap_id")

vertices_subt <- vertices %>% group_by(primary_disease, subtype_disease) %>% tally(sort = T) 

#### Lung Cancer -- 
## NSCLC-Adenocarcinoma
## NSCLC-Large Cell 
## SCLC

vertices %>% filter(primary_disease == "Lung Cancer") %>% 
  group_by(primary_disease, subtype_disease) %>% tally(sort = T) 

vertices[vertices$primary_disease == "Lung Cancer" &
           vertices$subtype_disease == "Non-Small Cell Lung Cancer (NSCLC), Adenocarcinoma",
         "subtype_disease"] <- "NSCLC-Adenocarcinoma"
vertices[vertices$primary_disease == "Lung Cancer" &
           vertices$subtype_disease == "Non-Small Cell Lung Cancer (NSCLC), Bronchoalveolar Carcinoma",
         "subtype_disease"] <- "NSCLC-Adenocarcinoma"
vertices[vertices$primary_disease == "Lung Cancer" &
           vertices$subtype_disease == "Non-Small Cell Lung Cancer (NSCLC), Large Cell Carcinoma",
         "subtype_disease"] <- "NSCLC-Large Cell"
vertices[vertices$primary_disease == "Lung Cancer" &
           vertices$subtype_disease == "Small Cell Lung Cancer (SCLC)",
         "subtype_disease"] <- "SCLC"

## Pancreatic Cancer
## Exocrine - Endocrine. We only have exocrine.
## Exocrine: Adenocarcinoma, acinar cell carcinoma of the pancreas, Pancreatoblastoma, 
## adenosquamous carcinomas, etc
vertices %>% filter(primary_disease == "Pancreatic Cancer") %>% 
  group_by(primary_disease, subtype_disease) %>% tally(sort = T) 

vertices[vertices$primary_disease == "Pancreatic Cancer" &
           vertices$subtype_disease == "Ductal Adenocarcinoma, exocrine",
         "subtype_disease"] <- "Adenocarcinoma"

vertices[vertices$primary_disease == "Pancreatic Cancer" &
           vertices$subtype_disease == "Ductal Adenosquamous Carcinoma, exocrine",
         "subtype_disease"] <- "Adenosquamous Carcinoma"

## Neuroblastoma
## No subtypes identified
vertices[vertices$primary_disease == "Neuroblastoma", "subtype_disease"] <- "Neuroblastoma"

## Leukemia 
## AML
## ALL
vertices %>% filter(primary_disease == "Leukemia") %>% 
  group_by(primary_disease, subtype_disease) %>% tally(sort = T) 

vertices[vertices$primary_disease == "Leukemia" &
           grepl(x =vertices$subtype_disease, pattern = "AML"),
         "subtype_disease"] <- "AML"
vertices[vertices$primary_disease == "Leukemia" &
           grepl(x=vertices$subtype_disease, pattern = "ALL"),
         "subtype_disease"] <- "ALL"
vertices[vertices$primary_disease == "Leukemia" &
           vertices$subtype_disease == "Acute Myelogenous Leukemia (AML), M5 (Monocytic)",
         "subtype_disease"] <- "AML"

vertices[vertices$primary_disease == "Leukemia" &
           vertices$subtype_disease == "Unknown",]
#ACH-001735 615        Leukemia         Unknown      <NA>      <NA>         <NA>         <NA>      <NA>
#ACH-001736 616        Leukemia         Unknown      <NA>      <NA>         <NA>         <NA>      <NA>

# https://depmap.org/portal/cell_line/ACH-001735
# ACH-001735 SEMK2 SEMK2_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE  
# Acute Lymphoblastic Leukemia (ALL), B-cell

# https://depmap.org/portal/cell_line/ACH-001736
# ACH-001736 HB1119
# Acute Lymphoblastic Leukemia (ALL), B-cell

vertices[vertices$primary_disease == "Leukemia" &
           vertices$subtype_disease == "Unknown", "subtype_disease"] <- "ALL"

vertices[vertices$primary_disease == "Leukemia" &
           vertices$subtype_disease == "T-cell, unspecified",]
# https://depmap.org/portal/cell_line/ACH-000953
# SUP-T1 Acute Lymphoblastic Leukemia (ALL), T-cell
vertices[vertices$primary_disease == "Leukemia" &
           vertices$subtype_disease == "T-cell, unspecified", "subtype_disease"] <- "ALL"

## Brain Cancer
vertices[vertices$primary_disease == "Brain Cancer" &
           vertices$subtype_disease == "Unknown",]
# https://depmap.org/portal/cell_line/ACH-000655
# ACH-000655 Glioma
vertices[vertices$primary_disease == "Brain Cancer" &
           vertices$subtype_disease == "Unknown", "subtype_disease"] <- "Glioma"

## Endometrial/Uterine Cancer
vertices[vertices$primary_disease == "Endometrial/Uterine Cancer",]

# https://depmap.org/portal/cell_line/ACH-000831
# ACH-000831 HEC-50B Endometrial Adenocarcinoma
vertices[vertices$primary_disease == "Endometrial/Uterine Cancer" & 
           vertices$depmap_id == "ACH-000831", "subtype_disease"] <- "Endometrial Adenocarcinoma"

# Head and Neck Cancer
vertices[vertices$primary_disease == "Head and Neck Cancer", ] 

vertices[vertices$primary_disease == "Head and Neck Cancer" & 
           grepl(pattern = "Squamous", vertices$subtype_disease), "subtype_disease"] <- "Squamous Cell Carcinoma"

# ACH-000836 Squamous Cell Carcinoma, tongue
# ACH-000832 Squamous Cell Carcinoma
# ACH-001346 Squamous Cell Carcinoma, buccal mucosa


## Lymphoma
vertices[vertices$primary_disease == "Lymphoma" &
           vertices$subtype_disease == "B-cell, Hodgkins", "subtype_disease"] <- "Hodgkins"
vertices[vertices$primary_disease == "Lymphoma" &
           vertices$subtype_disease == "B-cell, Non-Hodgkins", "subtype_disease" ] <- "Non-Hodgkins, B-cell"
vertices[vertices$primary_disease == "Lymphoma" &
           vertices$subtype_disease == "Diffuse Large B-cell Lymphoma (DLBCL)", 
         "subtype_disease" ] <- "Non-Hodgkins, B-cell"
vertices[vertices$primary_disease == "Lymphoma" &
           vertices$subtype_disease == "T-cell, Non-Hodgkins, Anaplastic Large Cell (ALCL)", 
         "subtype_disease" ] <- "Non-Hodgkins, T-cell"


### Kidney

vertices[vertices$primary_disease == "Kidney Cancer" &
           vertices$depmap_id == "ACH-000649", 
         "subtype_disease" ] <- "Renal Carcinoma, clear cell"
# https://depmap.org/portal/cell_line/ACH-000649
# Renal Carcinoma, clear cell

### Ovarian
## 	ACH-000278 Adenocarcinoma, high grade serous
vertices[vertices$primary_disease == "Ovarian Cancer" &
           vertices$depmap_id == "ACH-000278", 
         "subtype_disease" ] <- "Adenocarcinoma, high grade serous"
# 	ACH-000906 Adenocarcinoma, clear cell
vertices[vertices$primary_disease == "Ovarian Cancer" &
           vertices$depmap_id == "ACH-000906", 
         "subtype_disease" ] <- "Adenocarcinoma, clear cell"

vertices[vertices$primary_disease == "Ovarian Cancer" &
           grepl("Adenocarcinoma", vertices$subtype_disease) , 
         "subtype_disease" ] <- "Adenocarcinoma"

## Skin Cancer
vertices[vertices$primary_disease == "Skin Cancer" &
           vertices$depmap_id == "ACH-000799", 
         "subtype_disease" ] <- "Melanoma"

## Breast Cancer
vertices[vertices$primary_disease == "Breast Cancer" &
           vertices$subtype_disease == "Breast Ductal Carcinoma",
         "subtype_disease"] <- "Carcinoma"

vertices[vertices$primary_disease == "Colon/Colorectal Cancer",
         "subtype_disease"] <- "Adenocarcinoma"

vertices[vertices$primary_disease == "Endometrial/Uterine Cancer" &
         vertices$subtype_disease == "Endometrial Adenocarcinoma",
         "subtype_disease"] <- "Adenocarcinoma"

## Sarcoma
vertices[vertices$primary_disease == "Sarcoma" &
           grepl(pattern = "Rhabdomyosarcoma" ,vertices$subtype_disease),
         "subtype_disease"] <- "Rhabdomyosarcoma"

vertices_subt <- vertices %>% group_by(primary_disease, subtype_disease) %>% tally()
vertices %>% select(id, depmap_id, primary_disease, subtype_disease) %>%
  rename(Node_ID = id, DepMap_ID = depmap_id, 
         "Primary Disease" = primary_disease, "Subtype Disease" = subtype_disease) %>%
  write_tsv("data/cell_lines_corr-above0-5-0-323-subtypes.nodes")
