suppressMessages(library(dplyr))
library(data.table)
library(stringr)
#library(scales)

epitope_h1 <- scan("IEDB_epitope_full_v3.csv", nlines = 1, what = character(), sep = ",")
epitope_h2 <- scan("IEDB_epitope_full_v3.csv", skip = 1, nlines = 1, what = character(), sep = ",")
epitope_data <- fread("IEDB_epitope_full_v3.csv", skip = 2, header = FALSE, sep = ",")
names(epitope_data) <- paste(gsub("\\s", "", epitope_h1), gsub("\\s", "", epitope_h2), sep="_")

mhc_bind_h1 <- scan("IEDB_mhc_ligand_full.csv", nlines = 1, what = character(), sep = ",")
mhc_bind_h2 <- scan("IEDB_mhc_ligand_full.csv", skip = 1, nlines = 1, what = character(), sep = ",")
mhc_bind_data <- fread("IEDB_mhc_ligand_full.csv", skip = 2, header = FALSE, sep = ",")
names(mhc_bind_data) <- paste(gsub("\\s", "", mhc_bind_h1), gsub("\\s", "", mhc_bind_h2), sep="_")

tcr_h1 <- scan("IEDB_tcell_full_v3.csv", nlines = 1, what = character(), sep = ",")
tcr_h2 <- scan("IEDB_tcell_full_v3.csv", skip = 1, nlines = 1, what = character(), sep = ",")
tcr_data <- fread("IEDB_tcell_full_v3.csv", skip = 2, header = FALSE, sep = ",")
names(tcr_data) <- paste(gsub("\\s", "", tcr_h1), gsub("\\s", "", tcr_h2), sep="_")

bcr_h1 <- scan("IEDB_bcell_full_v3.csv", nlines = 1, what = character(), sep = ",")
bcr_h2 <- scan("IEDB_bcell_full_v3.csv", skip = 1, nlines = 1, what = character(), sep = ",")
bcr_data <- fread("IEDB_bcell_full_v3.csv", skip = 2, header = FALSE, sep = ",")
names(bcr_data) <- paste(gsub("\\s", "", bcr_h1), gsub("\\s", "", bcr_h2), sep="_")

neoEpitopes <- epitope_data[epitope_data$RelatedObject_EpitopeRelationship=="neo-epitope" & grepl("sapien",epitope_data$RelatedObject_OrganismName) & !grepl("[^A-Z]",epitope_data$Epitope_Description) & !is.na(epitope_data$RelatedObject_EndingPosition) & epitope_data$RelatedObject_ParentProteinIRI!="", c("Epitope_Description","RelatedObject_StartingPosition","RelatedObject_EndingPosition","RelatedObject_ParentProteinIRI","RelatedObject_Description")]
neoEpitopes$PID <- str_split_fixed(neoEpitopes$RelatedObject_ParentProteinIRI, "/",5)[,5]

#convert qualitative Measurements to number 0-5 and then take the median for multiple measurements per peptide
#Negative: 598478; Positive-Low: 53538; Positive-High: 5361; Positive: 3877268; Positive-Intermediate: 45681
neoEpitopes_mhc_bind <- left_join(neoEpitopes,mhc_bind_data,by="Epitope_Description")
mhc_pos_assays <- group_by(neoEpitopes_mhc_bind[neoEpitopes_mhc_bind$Assay_QualitativeMeasure!="Negative",],`Assay_Method/Technique`,Assay_AssayGroup,Assay_AssayTypeIRI) %>% summarize(freq=n())
mhc_neg_assays <- group_by(neoEpitopes_mhc_bind[neoEpitopes_mhc_bind$Assay_QualitativeMeasure=="Negative",],`Assay_Method/Technique`,Assay_AssayGroup,Assay_AssayTypeIRI) %>% summarize(freq=n())
mhc_assays_merged <- full_join(mhc_pos_assays,mhc_neg_assays, by=c("Assay_Method/Technique", "Assay_AssayGroup", "Assay_AssayTypeIRI"), suffix = c(".positive",".negative"))

mhc_bind_data$MHC_aff[mhc_bind_data$Assay_QualitativeMeasure=="Positive-High"] <- 3
mhc_bind_data$MHC_aff[mhc_bind_data$Assay_QualitativeMeasure %in% c("Positive-Intermediate","Positive")] <- 2
mhc_bind_data$MHC_aff[mhc_bind_data$Assay_QualitativeMeasure=="Positive-Low"] <- 1
mhc_bind_data$MHC_aff[mhc_bind_data$Assay_QualitativeMeasure=="Negative"] <- 0
mhc_bind_uniq <- group_by(mhc_bind_data[grepl("sapien",mhc_bind_data$Host_Name) & !grepl("[^A-Z]",mhc_bind_data$Epitope_Description)],Epitope_Description) %>% summarize(MHC_aff=median(MHC_aff),MHC_AlleleName=first(MHC_AlleleName), MHC_MHCalleleclass=first(MHC_MHCalleleclass))
#AvgScore <- data17 %>% mutate_at(vars(avg_score.y), ~ if_else(is.na(.), avg_score.x, .))
neoEpitopes_aff1 <- left_join(neoEpitopes,mhc_bind_uniq,by="Epitope_Description")
neoEpitopes_aff2 <- left_join(neoEpitopes_aff1,mhc_bind_uniq,by=c("RelatedObject_Description"="Epitope_Description"))

#---------------------------- Summary of tcell,mhc,bcell assays -----------------------------------
#Negative: 13986; Positive-Low: 270; Positive-High: 93; Positive: 2112; Positive-Intermediate: 61
neoEpitopes_tcell_response <- left_join(neoEpitopes,tcr_data,by="Epitope_Description")
tcell_pos_assays <- group_by(neoEpitopes_tcell_response[neoEpitopes_tcell_response$Assay_QualitativeMeasure!="Negative",],`Assay_Method/Technique`,Assay_AssayGroup,Assay_AssayTypeIRI) %>% summarize(freq=n())
tcell_neg_assays <- group_by(neoEpitopes_tcell_response[neoEpitopes_tcell_response$Assay_QualitativeMeasure=="Negative",],`Assay_Method/Technique`,Assay_AssayGroup,Assay_AssayTypeIRI) %>% summarize(freq=n())
tcell_assays_merged <- full_join(tcell_pos_assays,tcell_neg_assays, by=c("Assay_Method/Technique", "Assay_AssayGroup", "Assay_AssayTypeIRI"), suffix = c(".positive",".negative"))

#Negative: 28; Positive-Low: 4; Positive-High: 2; Positive: 70;
neoEpitopes_bcell_response <- left_join(neoEpitopes,bcr_data,by="Epitope_Description")
bcell_pos_assays <- group_by(neoEpitopes_bcell_response[neoEpitopes_bcell_response$Assay_QualitativeMeasure!="Negative",],`Assay_Method/Technique`,Assay_AssayGroup,Assay_AssayTypeIRI) %>% summarize(freq=n())
bcell_neg_assays <- group_by(neoEpitopes_bcell_response[neoEpitopes_bcell_response$Assay_QualitativeMeasure=="Negative",],`Assay_Method/Technique`,Assay_AssayGroup,Assay_AssayTypeIRI) %>% summarize(freq=n())
bcell_assays_merged <- full_join(bcell_pos_assays,bcell_neg_assays, by=c("Assay_Method/Technique", "Assay_AssayGroup", "Assay_AssayTypeIRI"), suffix = c(".positive",".negative"))


#----------------------------tcell,mhc,bcell assays per neoepitope-----------------------------------
mhc_pos_epi <- group_by(neoEpitopes_mhc_bind[neoEpitopes_mhc_bind$Assay_QualitativeMeasure!="Negative",],Epitope_Description) %>% summarize(freq=n())
mhc_neg_epi <- group_by(neoEpitopes_mhc_bind[neoEpitopes_mhc_bind$Assay_QualitativeMeasure=="Negative",],Epitope_Description) %>% summarize(freq=n())
mhc_epi_merged <- full_join(mhc_pos_epi,mhc_neg_epi, by="Epitope_Description", suffix = c(".mhc.positive",".mhc.negative"))

tcell_pos_epi <- group_by(neoEpitopes_tcell_response[neoEpitopes_tcell_response$Assay_QualitativeMeasure!="Negative",],Epitope_Description) %>% summarize(freq=n())
tcell_neg_epi <- group_by(neoEpitopes_tcell_response[neoEpitopes_tcell_response$Assay_QualitativeMeasure=="Negative",],Epitope_Description) %>% summarize(freq=n())
tcell_epi_merged <- full_join(tcell_pos_epi,tcell_neg_epi, by="Epitope_Description", suffix = c(".tcell.positive",".tcell.negative"))

bcell_pos_epi <- group_by(neoEpitopes_bcell_response[neoEpitopes_bcell_response$Assay_QualitativeMeasure!="Negative",],Epitope_Description) %>% summarize(freq=n())
bcell_neg_epi <- group_by(neoEpitopes_bcell_response[neoEpitopes_bcell_response$Assay_QualitativeMeasure=="Negative",],Epitope_Description) %>% summarize(freq=n())
bcell_epi_merged <- full_join(bcell_pos_epi,bcell_neg_epi, by="Epitope_Description", suffix = c(".bcell.positive",".bcell.negative"))

assays_per_epitope <- full_join(tcell_epi_merged,bcell_epi_merged,by="Epitope_Description") %>% full_join(mhc_epi_merged,by="Epitope_Description")
assays_per_epitope[is.na(assays_per_epitope)]<-0
#Create separate files for neoepitopes with tcell,bcell and mhc data
assays_per_epitope[assays_per_epitope$freq.tcell.positive>0,] #1404
assays_per_epitope[assays_per_epitope$freq.tcell.positive==0 & assays_per_epitope$freq.bcell.positive>0,] #17
assays_per_epitope[assays_per_epitope$freq.tcell.positive==0 & assays_per_epitope$freq.bcell.positive==0 & assays_per_epitope$freq.tcell.negative==0 & assays_per_epitope$freq.mhc.positive>0,] #316
neoEpitopes_tcell_bcell <- inner_join(neoEpitopes,assays_per_epitope[assays_per_epitope$freq.tcell.positive>0 | assays_per_epitope$freq.bcell.positive>0,],by="Epitope_Description")


#----------------------------aggregate by various combinations of assays and count the number of neoepitopes for each combination------------------------
#----------------------------don't forget to count the the neoepitopes with no assays; could use bitflags like samtools-------------------------------

fwrite(data.table::data.table(neoEpitopes), file = 'IEDB_neoepitopes_sapien.tab', row.names = F, col.names = F, quote = F, sep = '\t')
fwrite(data.table::data.table(neoEpitopes_tcell_bcell), file = 'IEDB_tcellbcell_neoepitopes_sapien.tab', row.names = F, col.names = F, quote = F, sep = '\t')
epitopes_mhc <- unique(mhc_bind_data$Epitope_Description[grepl("sapien",mhc_bind_data$Host_Name) & !grepl("[^A-Z]",mhc_bind_data$Epitope_Description) & !(mhc_bind_data$Assay_QualitativeMeasure %in% c("Negative"))])
fwrite(data.table::data.table(epitopes_mhc), file = 'IEDB_allPositive_sapien.tab', row.names = F, col.names = F, quote = F, sep = '\t')
