library(ggplot2)
library(ggbeeswarm)
library(data.table)
library(dplyr)
library(ggsignif)
library(stringr)
library(pheatmap)
library(dendextend)
library(mRMRe)
library(proxy)
library(ggfortify)
unlink(".RData")
rm(list=ls())
gc(full=TRUE)

GRCh38_allbases <- 3209286105

ds="HCC2n5tissue" # 5tissue HCC2 HCC2n5tissue
db="" #"/multi_list"
setwd(paste0("/data/hemberg/nullomers/IEDB/",ds,db))
data<-fread("all_epitopeDB_neoepitopes_readCounts.tsv")
data2 <- group_by(data[data$"#reads">=10 & data$"#nullomers">2,c(1,2)],sample_id) %>% summarise(neoepitope_count=n())
metadata1 <- fread("HCC2_GEO_metadata.csv",select=c("Run","STAGE","Bases","AvgSpotLen"))
metadata1$stat[metadata1$STAGE==""] = "3"
metadata1$stat[metadata1$STAGE!=""] = "4"
metadata1$STAGE[metadata1$STAGE==""] = "Healthy2"
metadata1$STAGE[metadata1$STAGE!=""] = "HCC2"

metadata2 <- fread("5tissue_GEO_metadata.csv")[,c("Run","disease","Bases","AvgSpotLen")]
metadata2$STAGE <- factor(sub(" .*", "", metadata2$disease))
metadata2$stat <- as.integer(as.factor(metadata2$STAGE))

metadata <- rbind(metadata1[,c("Run","STAGE","stat","Bases","AvgSpotLen")],metadata2[,c("Run","STAGE","stat","Bases","AvgSpotLen")])
metadata$readCov = metadata$Bases / GRCh38_allbases

metadata_in <- metadata[match(data2$sample_id,metadata$Run),]
data_statS <- cbind(metadata_in[,c("stat","readCov")],data2)
data$epitope_ID <- paste0(data$ENST_id , "_", data$neoepitopes )
data3 <- data.frame(dcast(data[data$"#reads">=10 & data$"#nullomers">2,], sample_id ~ epitope_ID, value.var = "#reads", fill="0"))
rownames(data3) <- (data3$sample_id)
metadata_in <- metadata[match(data3$sample_id,metadata$Run),] 
#metadata_in$status<- factor(metadata_in$STAGE, labels = c("grey30","lightskyblue","salmon","violetred","palegreen2"))
data4 <- log10(data3[,-1]/metadata_in$readCov)
data4[data4<1] <- 0

data4$sample_id <- rownames(data4)
data_stat <- cbind(metadata_in[,"stat"],data4)

#plotting

data_statS$stat = as.character(data_statS$stat)

plot1 = ggplot(data_statS, aes(x=stat,y=log10(neoepitope_count/readCov))) + theme_classic() + theme(plot.title = element_text(size = 12)) +
geom_quasirandom() + 
geom_signif(comparisons=list(c("3","4")), map_signif_level = function(p) sprintf("p-value = %.2g", p)) +
stat_summary(fun = mean, geom = 'point', shape = 95, size = 10) +
stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', width = 0, size = 0.8) +
labs(title="Neoepitopes per patient (>=3 reads)") + 
scale_x_discrete(labels=c("Colorectal", "Esophagus", "Healthy", "Liver", "Lung", "Stomach")) +
xlab("Patients")+ ylab("#geneFusions / #reads")

plot2 = ggplot(data_statS, aes(x=stat,y=log10(neoepitope_count/readCov))) + theme_classic() + theme(plot.title = element_text(size = 12)) +
geom_quasirandom() + labs(title="Neoepitopes per patient (>=30 reads)") + 
scale_x_continuous(breaks=c(1,2,3,4,5,6), labels=c("Colorectal", "Esophagus", "Healthy", "Liver", "Lung", "Stomach")) +
xlab("Patients")+ ylab("#geneFusions / #reads")

            
plots_ph=NULL
plots_all=NULL
#select cancer in in "5tissues" dataset 
dataset <- c("CRC", "Esophagus", "Healthy", "HCC", "LUNG", "Stomach")

#outliers[1] <- c(12,52,65,90)   # -ouliers[x]
for (x in c(1,2,4,5,6)){ # start of function
	print(x)
	data_merged <- data_stat[data_stat$stat %in% c(x,3),]
	data_merged$stat[data_merged$stat==x] <- 1
	data_merged$stat[data_merged$stat==3] <- 0 
#	max_neo <- data_merged %>% group_by(stat) %>% summarize(c=n()) %>% select("c") %>% min()*0.8 #filter features in more than 80% patients 
	neo_all <- colSums(data_merged[,-c("stat","sample_id")])
#	patient_all <- rowSums(data_merged[,-c("stat","sample_id")])
	neo_inc <- neo_all[neo_all >0]
#	patient_inc <- which(patient_all < 3000)
	data_in <- data_merged[] %>% select(c("sample_id","stat",names(neo_inc)))
	metadata_in <- metadata[match(data_in$sample_id,metadata$Run),]
#	metadata_in$STAGE[metadata_in$stat==] <- factor(sub(" .*", "", metadata_in$disease))
#	metadata_in$stat <- as.integer(as.factor(metadata_in$STAGE))
	data_in <- data_in[,-c("sample_id")]
	data_in[] <- lapply(data_in, as.numeric)
	
	cancerCount <- length(data_in$stat[data_in$stat==1])
	healthyCount <- length(data_in$stat[data_in$stat==0])
	neoepitopes <- data.frame(colSums(data_in[data_in$stat==1]), colSums(data_in[data_in$stat==0]))
	colnames(neoepitopes) <- c("cancerCount","healthyCount")
	neoepitopes <- neoepitopes[-1,]
	neoepitopes$cancerMeans <- neoepitopes$cancerCount/cancerCount
	neoepitopes$healthyMeans <- neoepitopes$healthyCount/healthyCount
	neoepitopesF <- neoepitopes[(neoepitopes$cancerMeans/neoepitopes$healthyMeans<0.25 | neoepitopes$cancerMeans/neoepitopes$healthyMeans>4) & (neoepitopes$cancerMeans>0.05 | neoepitopes$healthyMeans>0.05),]
	data5 <- data.frame(data_in)[,c(rownames(neoepitopesF))]
	names_tmp <- str_split_fixed(colnames(data5), "[..]", n=3)
	colnames(data5)<-names_tmp[,1]
	ensts <- str_split_fixed(c(names_tmp[,1]), "_", n=2)[,1]
	enst_ann<-data.frame(unique(data[data$ENST_id %in% ensts, c("ENST_id","neoepitopes","annotation","db_name")]))
	enst_ann$neoepitopes <- gsub("[->;]", ".", enst_ann$neoepitopes)
	rownames(enst_ann)<- paste0(enst_ann$ENST_id,"_",enst_ann$neoepitopes)
	enst_ann$neoepitopes <- str_split_fixed(c(enst_ann$neoepitopes), "[..]", n=2)[,1]
	top_neoepitopes <- merge(neoepitopesF,enst_ann,by="row.names",all.x=TRUE)
	ann_row <- data.frame(data_in[,1])
	ann_row$stat = as.character(ann_row$stat)
	colnames(ann_row) <- c("cancer status")
	plots_ph[[x]] <- pheatmap(as.matrix(data5,rownames.force = TRUE), cluster_rows = F, clustering_distance_cols="euclidean",clustring_method="ward.D2", annotation_row = ann_row, annotation_colors = list(stat = c("0"="#00BFC4","1"="salmon")))
	fwrite(data.table::data.table(top_neoepitopes[,-1]), file = paste0(ds,"_topAntigens_logNormalizedReads_",dataset[x],".tsv"), row.names = F, col.names = T, quote = F, sep = '\t')

	#feature selection mRMR 0.005 or 0.009
	data_in$stat = as.numeric(data_in$stat)
	numFeatures <- ncol(data_in)-1 #ncol(data_in)-1 #75
	mRMR_o <- mRMR.data(data = data.frame(data_in[])) %>% 
	mRMR.classic("mRMRe.Filter", data = ., target_indices = 1, feature_count = numFeatures, method="exhaustive")
	col_sel <- scores( mRMR_o)[[1]]
	keep_cols <- c(1,solutions(mRMR_o)[[1]][col_sel > 0.005 & !is.na(col_sel)])
	# #wilcox t-test
	# data_in$stat = as.factor(data_in$stat)
	# wilcox_o <- lapply(colnames(data_in[,-1]), function(x){
	# t_result <- wilcox.test(as.formula(paste0(x, " ~ stat")), data=data_in)
	# return(c(x,t_result$p.value))
	# })
	# wilcox_df <- data.frame(matrix(unlist(wilcox_o), nrow=length(wilcox_o), byrow=TRUE),stringsAsFactors=FALSE)
	# keep_cols <- c("stat",wilcox_df$X1[as.numeric(wilcox_df$X2)<=0.1])

	pc <- prcomp(data_in[,-c("stat")], center = TRUE)
	plots_all[[paste0(x,"_1")]] <- cbind(data_in, metadata_in[,"stat"]) %>% autoplot(pc, data= ., colour = 'stat') + theme_classic() + 
	labs(title=paste0("PCA of all features in ",dataset[x])) 
	
	data_sub <- data_in[, ..keep_cols]

	pc_sub <- prcomp(data_sub[,-c("stat")], center = TRUE)
	plots_all[[paste0(x,"_2")]] <- cbind(data_sub, metadata_in[,"STAGE"]) %>% autoplot(pc_sub, x=1, y=2, data= ., colour = 'STAGE')+ theme_classic() +
	labs(title=paste0("PCA of ",length(keep_cols)-1," selected features in ",dataset[x]))
	
	#Hierarchical Clustering
	metadata_in$status<- factor(metadata_in$stat, labels = c("salmon","#00BFC4"))
	distance_mat <- dist(data_sub[,-c("stat")], method = 'eJaccard')
	set.seed(240)  # Setting seed
	Hierar_cl <- hclust(distance_mat, method = "ward.D2") #, ward.D2
	dend <- as.dendrogram(Hierar_cl)  %>% set("labels_cex",0.3)
	metadataHC <- metadata_in[Hierar_cl$order,]
	dend2 <- color_labels(dend,col=as.character(metadataHC$status))
	pdf(NULL)
	dev.control(displaylist="enable")
	plot(dend2, cex.axis = 0.4,cex.main=0.6, main=paste0("Heirarchical clustering using eJaccard distance of ",length(keep_cols)-1," selected features"))
	legend("topright", legend=levels(metadataHC$STAGE),fill=levels(metadataHC$status), cex = 0.4) #, title="Condition"
	plots_all[[paste0(x,"_3")]] <- recordPlot()
	invisible(dev.off())
}

pdf(paste0(ds,"_logNormalizedReads_per_fusion_v3.pdf"), onefile=TRUE, width=5.5, height=4)
plot1
plot2
plots_all
dev.off()

pdf(paste0(ds,"_logNormalized_epitopeDB_heatmap.pdf"), onefile=TRUE, width=10, height=10)
for (x in c(1,2,4,5,6)){
	grid::grid.newpage()
	grid::grid.draw(plots_ph[[x]]$gtable)
}
dev.off()

