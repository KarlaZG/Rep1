library(readr)
library(dplyr)
library(phyloseq)
library(ggplot2)
library(microbiome)
library(vegan)
library(hillR)
library(factoextra)
library(FactoMineR)
library(ALDEx2)
library(tibble)
library(readxl)
library(iNEXT)
library(xlsx)
library(tidyr)

#Importar la otu table, taxonomía y metadata

setwd("C:/Users/karla/OneDrive/Documents/4_PARCELAS/barplot")
otu_data <- read_tsv("merged-table.tsv")
#otu_data <- read_tsv("merged-filtered-table.tsv")
tax_data <- read_tsv("merged_taxonomy_trained.tsv")
tax_data<- tax_data %>% 
  #separamos columna taxonomia en los diferentes niveles
  separate(Taxon, c("Kingdom", "Phylum", "Class","Order", "Family", "Genus","Specie"),  sep = ";")
tax_data<- tax_data[,1:8] 
sample_data <- read_tsv("metadata_trained.tsv",
                        col_types = "ccccccccc")

otu_data <- otu_data %>%
  tibble::column_to_rownames("OTU ID") 
tax_data <- tax_data %>% 
  tibble::column_to_rownames("OTU ID")
sample_data <- sample_data %>% 
  tibble::column_to_rownames("Sample ID")

#transformar en matrices 

otu_data <- as.matrix(otu_data)
tax_data <- as.matrix(tax_data)

#transformar en objetos phyloseq 
OTU = otu_table(otu_data, taxa_are_rows = TRUE)
TAX = tax_table(tax_data)
samples = sample_data(sample_data)

#unir en un solo objeto
D0_D100_D100R<- phyloseq(OTU, TAX, samples)

#filtrar ASVs que no están asignadas a nivel de género
filterPhyla = c(NA, "p__")
D0_D100_D100R<- subset_taxa(D0_D100_D100R, !Phylum %in% filterPhyla)

#Convertir la otu table a un data frame
D0_D100_D100R_OTUTABLE<- otu_table(D0_D100_D100R)
D0_D100_D100R_OTU <- data.frame(D0_D100_D100R_OTUTABLE)
summary(sample_sums(D0_D100_D100R))

#Colapsar a nivel de género
TODO_GENERO <- tax_glom(D0_D100_D100R, "Genus", NArm = FALSE)
TODO_GENUS <- otu_table(TODO_GENERO)
Gen <- data.frame(TODO_GENUS)

#Obtener la taxonomía solo de género
TAXA <- tax_table(TODO_GENERO)
taxonomia_genus <- data.frame(TAXA[,6])
taxonomia_genus<- rownames_to_column(taxonomia_genus, "OTU ID")


GENERO_TODO <- data.frame(taxonomia_genus, Gen)
GENERO_TODO <- GENERO_TODO[,2:49]

#Asignar cada muestra a una variable 
data<- GENERO_TODO[2:48]

TAG0<- data[,1:3]
TAG0<-rownames_to_column(TAG0, "OTU ID")
TC0<- data[,4:6]
TC0<-rownames_to_column(TC0, "OTU ID")
TFE0<- data[,7:9]
TFE0<-rownames_to_column(TFE0, "OTU ID")
TSAG0<- data[,10:12]
TSAG0<-rownames_to_column(TSAG0, "OTU ID")
TSFE0<- data[,13:15]
TSFE0<-rownames_to_column(TSFE0, "OTU ID")

TAG100<- data[,16:19]
TAG100<-rownames_to_column(TAG100, "OTU ID")
TC100<- data[,20:23]
TC100<-rownames_to_column(TC100, "OTU ID")
TFE100<- data[,24:27]
TFE100<-rownames_to_column(TFE100, "OTU ID")
TSAG100<- data[,28:31]
TSAG100<-rownames_to_column(TSAG100, "OTU ID")
TSFE100<- data[,32:35]
TSFE100<-rownames_to_column(TSFE100, "OTU ID")

TAG100R<- data[,36:39]
TAG100R<-rownames_to_column(TAG100R, "OTU ID")
TC100R<- data[,40:43]
TC100R<-rownames_to_column(TC100R, "OTU ID")
TFE100R<- data[,44:47]
TFE100R<-rownames_to_column(TFE100R, "OTU ID")

setwd("C:/Users/karla/OneDrive/Documents/4_PARCELAS/barplot/subgrupos_counts")

#Hacer los subgrupos y exportarlos como archivo de excel
TAG100_TC100<- merge(TAG100,TC100)
write.xlsx(x=TAG100_TC100,file ="ag_c_bulk_counts.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE )

TAG100R_TC100R<- merge(TAG100R,TC100R)
write.xlsx(x=TAG100R_TC100R,file ="ag_c_riz_counts.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE )

TFE100_TC100<- merge(TC100,TFE100)
write.xlsx(x=TFE100_TC100,file ="fe_c_bulk_counts.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE )

TFE100R_TC100R<- merge(TC100R,TFE100R)
write.xlsx(x=TFE100R_TC100R,file ="fe_c_riz_counts.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE )

AG_C<- merge(TAG100,TC100)
NORIZ<- merge(AG_C,TFE100)
write.xlsx(x=NORIZ,file ="ag_c_FE_BULK_COUNTS.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE )

AGR_CR<- merge(TAG100R,TC100R)
RIZ<- merge(AGR_CR,TFE100R)
write.xlsx(x=RIZ,file ="ag_c_FE_RIZ_COUNTS.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE )

TSAG_TSFE<- merge(TSAG100,TSFE100)
write.xlsx(x=TSAG_TSFE,file ="ag_FE_NOCUL_COUNTS.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE )

T100R<- merge (TAG100R,TFE100R)
TS100<- merge (TSAG100,TSFE100)
INTERACCION<- merge(T100R,TS100)
write.xlsx(x=INTERACCION,file ="interaccion_RIZ_NOCULT_COUNTS.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE )

TSAG100_TAG100R <- merge(TSAG100,TAG100R)
write.xlsx(x=TSAG100_TAG100R,file ="TSAG_TAGR_COUNTS.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE )

TSFE100_TFE100R <- merge(TSFE100,TFE100R)
write.xlsx(x=TSFE100_TFE100R,file ="TSFE_TFER_counts.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE )

tiempo_control<- merge(TC0,TC100)
tiempo_control<-merge(tiempo_control,TC100R)
write.xlsx(x=tiempo_control,file ="TC0_TC100_counts.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE )

tiempo_plata<- merge(TAG0,TSAG0)
tiempo_plata<-merge(tiempo_plata,TAG100)
tiempo_plata<-merge(tiempo_plata,TSAG100)
tiempo_plata<-merge(tiempo_plata,TAG100R)
write.xlsx(x=tiempo_plata,file ="TAG0_TAG100_counts.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE )

tiempo_fierro<- merge(TFE0,TSFE0)
tiempo_fierro<-merge(tiempo_fierro,TFE100)
tiempo_fierro<-merge(tiempo_fierro,TSFE100)
tiempo_fierro<-merge(tiempo_fierro,TFE100R)
write.xlsx(x=tiempo_fierro,file ="TFE0_TFE100_counts.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE )
