library(ggplot2)
library(ggbeeswarm)
library(data.table)
library(class)
suppressMessages(library(dplyr))
library(caret)
library(glmnet)
library(caTools)
library(randomForest)
library(e1071)
library(ggfortify)
#library(mRMRe)
library(neuralnet)
#library(xgboost)
library(stringr)
library(pROC)
unlink(".RData")
rm(list=ls())
gc(full=TRUE)

#genome size for normalizing expression
GRCh38_allbases <- 3209286105
# generic neoepitoes selected using 2-fold expression and occurance based tests 
generic_neoEpi <- c("CALHM6_AQISAAAAL","CAVIN2_NDQEEESFAEGHAEASLASALVEGE","CXCR2_IYAFIGQKFC","MT_CO1_HIVTYYPGK", "MT_CO3_FTSKHHFSF","MT_CYB_GILALLLSI", "MT_CYB_ITNLLFAIPY","MT_ND4_HNTPGSLNI", "MT_ND4_SSLNILLLTL","ITGB3_ITIHDQKEF","PIP4K2A_IYIDDNSKKVF")

# Read the biomarker data in form of TEs, circRNAs, and neoepitopes. Normalize and transform the expression of each modality as specified in the publication
setwd("cfRNA/backup")
ds=""
dataTE<-fread(paste0(ds,"all_cancer_SalmonTE_quant.tsv"))

#dataset selection
ds="HCC2n5tissue" # 5tissue HCC2 HCC2n5tissue pancreatic

#setwd(paste0("/data/hemberg/nullomers/IEDB/",ds))
dataCirc<-fread(paste0("all_",ds,"_circ_quant.tsv"))
data1<-fread(paste0("all_",ds,"_epitopeDB_neoepitopes_readCounts6.tsv"))
data1$rCount <- data1$"#reads"
data1$firstSNeo <- sapply(data1$neoepitopes, function(i) first(sort(str_split(i, ";", simplify=T))))
data1$epitope_ID2 <- gsub("-", "_",paste0(data1$HGNC_symbol,"_",str_split(data1$firstSNeo, "->", simplify=T)[,2]),fixed = TRUE)
data1$epitope_ID2 <- as.factor(data1$epitope_ID2)
data <- data1[!(data1$epitope_ID2 %in% generic_neoEpi)]

#HCC2 metadata
metadata1 <- fread("metadata_HCC2_cfRNA.csv",select=c("Run","STAGE","Bases","AvgSpotLen"))
metadata1$stat[metadata1$STAGE==""] = "3"
metadata1$stat[metadata1$STAGE!=""] = "4"
metadata1$STAGE[metadata1$STAGE!=""] = "Liver2"
metadata1$STAGE[metadata1$STAGE==""] = "Healthy2"

#HCCnMM metadata
metadata3 <- fread("metadata_HCCnMM_cfRNA.csv")[,c("Run","source_name","Bases","AvgSpotLen")]
metadata3$stat[grepl("non-cancer",metadata3$source_name)] = 3
metadata3$stat[grepl("liver cancer",metadata3$source_name)] = 4
metadata3$stat[grepl("multiple",metadata3$source_name)] = 7
metadata3$STAGE[grepl("non-cancer",metadata3$source_name)] = "Healthy3"
metadata3$STAGE[grepl("liver cancer",metadata3$source_name)] = "Liver3"
metadata3$STAGE[grepl("multiple",metadata3$source_name)] = "MMyeloma"

#5tissue metadata
metadata2 <- fread("metadata_HCCn5tissue_cfRNA.csv")[,c("Run","disease","Bases","AvgSpotLen")]
metadata2$STAGE <- factor(sub(" .*", "", metadata2$disease))
metadata2$stat <- as.integer(as.factor(metadata2$STAGE))

#Calculate coverage for normalization and flag the smaples with no TEs or Neoepitopes
metadata <- rbind(metadata1[,c("Run","STAGE","stat","Bases","AvgSpotLen")], metadata2[,c("Run","STAGE","stat","Bases","AvgSpotLen")], metadata3[,c("Run","STAGE","stat","Bases","AvgSpotLen")])
metadata$readCov = metadata$Bases / 10/ GRCh38_allbases
metadata$readCovTE = metadata$Bases/metadata$AvgSpotLen/GRCh38_allbases/10
missingN <- setdiff(metadata$Run,unique(data$sample_id))
missingT <- setdiff(metadata$Run,unique(dataTE$sample_id))

#Normalize, logscale and filter neoepitopes
data2 <- left_join(data,metadata[,c("Run","stat","readCov")], by= c("sample_id"="Run"))
data2$nCov <- data2$rCount/data2$readCov
data2$nCov[data2$nCov < 1] <- 1
data2$lnCov <- as.numeric(log10(data2$nCov))
data_stat <- dcast(data2[data2$rCount>=3 & data2$"#nullomers">=1,], sample_id + stat ~ epitope_ID2, fun=max, value.var = "lnCov", fill="0")

#Normalize, logscale and filter TEs
data2TE <- left_join(dataTE,metadata[,c("Run","stat","readCovTE")], by= c("sample_id"="Run"))
data2TE$RKP <- (data2TE$NumReads/data2TE$readCovTE)/data2TE$EffectiveLength
data2TE$RKP[data2TE$RKP < 1] <- 1
data2TE$logRKP <- as.numeric(log10(data2TE$RKP))
data2TE$Name <- gsub("-", "_", data2TE$Name,fixed = TRUE)
data_statTE <- dcast(data2TE, sample_id + stat ~ Name, value.var = "logRKP", fill="0") #fun=mean

#Normalize and logscale circRNAs
data2Circ <- left_join(dataCirc[dataCirc$bsj>=2 & dataCirc$circ_type == "exon",],metadata[,c("Run","stat")], by= c("sample_id"="Run"))
data2Circ$CPM <- data2Circ$CPM * 10
data2Circ$CPM[data2Circ$CPM < 1] <- 1
data2Circ$logRKP <- as.numeric(log10(data2Circ$CPM))
data_statCirc <- dcast(data2Circ, sample_id + stat ~ Name, fun=max, value.var = "logRKP", fill="0") #fun=mean

dataset <- c("Colorectal", "Esophagus", "Healthy", "HCC", "Lung", "Stomach", "M.Myeloma", "all")
gROCplots<-NULL

#Loop across the cancer-types in the "dataset" array
for (x1 in c(1,4)){ # start of function
	print(x1)
	#add samples with no detectable neoepitopes, TEs, or circRNAs to simplify the merger later
	missingN <- setdiff(metadata$Run,unique(data_stat$sample_id))
	missingT <- setdiff(metadata$Run,unique(data_statTE$sample_id))	
	missingC <- setdiff(metadata$Run,unique(data_statCirc$sample_id))	
	addRowT<-metadata[metadata$Run %in% missingT & metadata$stat %in% c(x1,3),c("Run","stat")]
	addDfT <- cbind(addRowT, data.frame(matrix(0, ncol = ncol(data_statTE)-2, nrow = nrow(addRowT))))
	colnames(addDfT)<-colnames(data_statTE)
	addRowC<-metadata[metadata$Run %in% missingC & metadata$stat %in% c(x1,3),c("Run","stat")]
	addDfC <- cbind(addRowC, data.frame(matrix(0, ncol = ncol(data_statCirc)-2, nrow = nrow(addRowC))))
	colnames(addDfC)<-colnames(data_statCirc)
	addRowN<-metadata[metadata$Run %in% missingN & metadata$stat %in% c(x1,3),c("Run","stat")]
	addDfN <- cbind(addRowN, data.frame(matrix(0, ncol = ncol(data_stat)-2, nrow = nrow(addRowN))))
	colnames(addDfN)<-colnames(data_stat)
#	do.call("rbind",lapply(x , function(ctype){})) 
	if(x1==8){
	        data_merged <- data_stat
	        data_merged$stat[data_merged$stat==3] <- 0
	        data_merged$stat[data_merged$stat>0] <- 1
	}else{
	        data_merged <- rbind(data_stat[data_stat$stat %in% c(x1,3),],addDfN)
	        data_merged$stat[data_merged$stat==x1] <- 1
	        data_merged$stat[data_merged$stat==3] <- 0
			data_mergedTE <- rbind(data_statTE[data_statTE$stat %in% c(x1,3),],addDfT)
			data_mergedTE$stat[data_mergedTE$stat==x1] <- 1
			data_mergedTE$stat[data_mergedTE$stat==3] <- 0		
			data_mergedCirc <- rbind(data_statCirc[data_statCirc$stat %in% c(x1,3),],addDfC)
			data_mergedCirc$stat[data_mergedCirc$stat==x1] <- 1
			data_mergedCirc$stat[data_mergedCirc$stat==3] <- 0		
	}

	#keep neoepitopes found in 1 or more patients after filtering
	neo_all <- colSums(data_merged[,-c("stat","sample_id")])
	neo_inc <- neo_all[neo_all >0]
	data_in <- data_merged[] %>% select(c("sample_id","stat",names(neo_inc))) %>% arrange(sample_id)
#	data_in <- data_in[,-c("sample_id")]
#	data_in[] <- lapply(data_in, as.numeric)

	#keep TEs found in 1 or more patients after filtering
	neo_allTE <- colSums(data_mergedTE[,-c("stat","sample_id")])
	neo_incTE <- neo_allTE[neo_allTE >0]
	data_inTE <- data_mergedTE[] %>% select(c("sample_id","stat",names(neo_incTE))) %>% arrange(sample_id)
#	data_inTE <- data_inTE[,-c("sample_id")]
#	data_inTE[] <- lapply(data_inTE, as.numeric)
	#keep circRNAs found in 1 or more patients after filtering
	neo_allCirc <- colSums(data_mergedCirc[,-c("stat","sample_id")])
	neo_incCirc <- neo_allCirc[neo_allCirc >0]
	data_inCirc <- data_mergedCirc[] %>% select(c("sample_id","stat",names(neo_incCirc))) %>% arrange(sample_id)

  # Specify the methods and parameters  
	methods <- c("Ridge","Lasso","SVM","RF") # "Ridge","Lasso","NN","SVM","RF"
	nRepeats <- 10 #iterations to run
	K <- 5  # K-fold cross validation
#	fsType <- "log2fold" #"PCA" #"mRMR" "LVQ" "wCox" "xgb" "log2fold"
	#softplus <- function(x) log(1 + exp(x))
	lambda_seq <- 10^seq(2, -2, by = -.1)
	depthTree <- 2
#	layersNN <- c(200,150,20,20)
	layersNN <- c(50,50,10)
	minXgbGain <- 1e-4
	wCoxPval <- 1e-7
	wCox_minN <- 5
	wCoxPvalC <- 1e-4
	wCox_minNC <- 5
	minDonors <- 3 #log2fold

	mat_sensiv <- data.frame(matrix(ncol = K+1))
	mat_specif <- data.frame(matrix(ncol = K+1))
	mat_accuracy <- data.frame(matrix(ncol = K+1))
	mat_auroc <- data.frame(matrix(ncol = K+1))
	mat_F1 <- data.frame(matrix(ncol = K+1))
	rocList <- list()
  
  # run 10 iterations of 5-fold CV
	for (y in 1:nRepeats) {
		print(paste0("Round ",y," :: ",K,"-fold"))
		#Create a fresh set of folds with randomized patients from each class
		#in creating the folds we specify the target feature (dependent variable) and # of folds
		folds <- createFolds(data_in$stat, k = K)
		fldc<-0
		# in cv we are going to applying multiple classifiers to our 'folds'
		cv <- lapply(folds, function(x) { # start of function
		  # in the next two lines we will separate the Training set into it's 10 pieces
	#	  x <- folds[[1]]
		  fldc<<-fldc+1
#		  print(x)
		  training_fold_full <- data_in[-x, ] # training fold is minus (-) it's sub test fold
		  test_fold_full <- data_in[x, ] # here we describe the test fold individually
		  training_fold_fullT <- data_inTE[-x, ] # training fold is minus (-) it's sub test fold
		  test_fold_fullT <- data_inTE[x, ] # here we describe the test fold individually
		  training_fold_fullC <- data_inCirc[-x, ] # training fold is minus (-) it's sub test fold
		  test_fold_fullC <- data_inCirc[x, ] # here we describe the test fold individually

# Script for excluding or including samples and choosing train/test sets by sample id START
#		  test_fold <- paste0("SRR",seq(from=10822540, to=10822604))
#		  test_fold <- paste0("SRR",seq(from=14506659, to=14506888))
#		  test_fold <- paste0("SRR",seq(from=15618984, to=15619034))
		  
#		  length(which(colSums(test_fold != 0) == 0))
#		  counttC <- length(test_foldT$stat[test_foldT$stat==1])
#		  counttH <- length(test_foldT$stat[test_foldT$stat==0])
#		  zeroes <- data.frame(cbind(colSums(test_foldT[test_foldT$stat==0,-c(1,2)]== 0),colSums(test_foldT[test_foldT$stat==1,-c(1,2)]== 0)))
#		  colnames(zeroes) <- c("zeroH","zeroC")
#		  rownames(zeroes) <- colnames(test_foldT[,-c(1,2)])
#		  keep_list0t <- rownames(zeroes[zeroes$zeroH<=(counttH * 0.2) & zeroes$zeroC<=(counttC * 0.2),])	
		  		  
#		  training_fold_full <- data_in[!(data_in$sample_id %in% test_fold), ] # training fold is minus (-) it's sub test fold
#		  test_fold_full <- data_in[data_in$sample_id %in% test_fold, ] # here we describe the test fold individually
#		  training_fold_fullT <- data_inTE[!(data_in$sample_id %in% test_fold), ] # training fold is minus (-) it's sub test fold
#		  test_fold_fullT <- data_inTE[data_in$sample_id %in% test_fold, ] # here we describe the test fold individually

#			remove all neoepitopes that are all zero in test fold
#			remove all TE present in less that 80% of the test fold
		  
#			run feature selection for the training data in each fold
#		  	(fsType == "xgb"){
#			xgboost based feature selection
#			xgb_train = xgb.DMatrix(data = as.matrix(training_fold_full[,-c("stat")]), label = as.matrix(training_fold_full$stat))
#	        xgb_test = xgb.DMatrix(data = as.matrix(test_fold_full[,-c("stat")]), label = as.matrix(test_fold_full$stat))
#			watchlist = list(train=xgb_train, test=xgb_test)
#			XGmodel = xgb.train(data = xgb_train, max.depth = 2, watchlist=watchlist, nrounds = 70, nthread = 10, importance=T)
#		    impFall <- xgb.importance( model = XGmodel)
#		    impF <- impFall[impFall$Gain>minXgbGain,]
#		    FT <- c("stat",impF$Feature)
#			print(length(FT))
#			training_fold <- training_fold_full[,..FT]
#			test_fold <- test_fold_full[,..FT]

# Script for excluding or including samples and choosing train/test sets by sample id END

# Feature selection for TEs, circRNAs, and neoepitopes done seperately for each fold after separating the test and training set
#		   (fsType == "wCox"){.    # TE
			training_fold_fullT$stat <- as.factor(training_fold_fullT$stat)
			test_fold_fullT$stat <- as.factor(test_fold_fullT$stat)
			countAllC <- length(training_fold_fullT$stat[training_fold_fullT$stat==1])
			countAllH <- length(training_fold_fullT$stat[training_fold_fullT$stat==0])
			zeroes <- data.frame(cbind(colSums(training_fold_fullT[training_fold_fullT$stat==0,-c(1,2)]== 0),colSums(training_fold_fullT[training_fold_fullT$stat==1,-c(1,2)]== 0)))
			colnames(zeroes) <- c("zeroH","zeroC")
			rownames(zeroes) <- colnames(training_fold_fullT[,-c(1,2)])
			keep_list0 <- rownames(zeroes[zeroes$zeroH<=(countAllH * 0.2) & zeroes$zeroC<=(countAllC * 0.2),])	
			wilcox_o <- lapply(keep_list0, function(x){
				t_result <- pairwise.wilcox.test(training_fold_fullT[[x]], training_fold_fullT$stat, p.adjust.method="bonf")
				return(c(x,t_result$p.value))
			})
			wilcox_df <- data.frame(matrix(unlist(wilcox_o), nrow=length(wilcox_o), byrow=TRUE),stringsAsFactors=FALSE)
	        wilcox_df$X2 <- as.numeric(wilcox_df$X2)
			wilcox_df <- wilcox_df[as.numeric(wilcox_df$X2) < wCoxPval,]
	        wilcox_minN <- slice_min(wilcox_df, n=as.numeric(wCox_minN), X2)
#			print(wilcox_minN)
			keep_cols <- c("sample_id","stat",wilcox_minN$X1)
			training_foldT <- training_fold_fullT[,..keep_cols]
			test_foldT <- test_fold_fullT[,..keep_cols]

#		   (fsType == "wCox"){.    # circRNA
			training_fold_fullC$stat <- as.factor(training_fold_fullC$stat)
			test_fold_fullC$stat <- as.factor(test_fold_fullC$stat)
			countAllCC <- length(training_fold_fullC$stat[training_fold_fullC$stat==1])
			countAllCH <- length(training_fold_fullC$stat[training_fold_fullC$stat==0])
			zeroesC <- data.frame(cbind(colSums(training_fold_fullC[training_fold_fullC$stat==0,-c(1,2)]== 0),colSums(training_fold_fullC[training_fold_fullC$stat==1,-c(1,2)]== 0)))
			colnames(zeroesC) <- c("zeroH","zeroC")
			rownames(zeroesC) <- colnames(training_fold_fullC[,-c(1,2)])
			keep_list0C <- rownames(zeroesC[zeroesC$zeroH<=(countAllCH * 0.6) & zeroesC$zeroC<=(countAllCC * 0.6),])	
			wilcoxC_o <- lapply(keep_list0C, function(x){
				t_result <- pairwise.wilcox.test(training_fold_fullC[[x]], training_fold_fullC$stat, p.adjust.method="bonf")
				return(c(x,t_result$p.value))
			})
			wilcoxC_df <- data.frame(matrix(unlist(wilcoxC_o), nrow=length(wilcoxC_o), byrow=TRUE),stringsAsFactors=FALSE)
	        wilcoxC_df$X2 <- as.numeric(wilcoxC_df$X2)
			wilcoxC_df <- wilcoxC_df[as.numeric(wilcoxC_df$X2) < wCoxPvalC,]
	        wilcoxC_minN <- slice_min(wilcoxC_df, n=as.numeric(wCox_minNC), X2)
#			print(wilcox_minN)
			keep_colsC <- c("sample_id","stat",wilcoxC_minN$X1)
			training_foldC <- training_fold_fullC[,..keep_colsC]
			test_foldC <- test_fold_fullC[,..keep_colsC]

#		  (fsType == "log2fold"){      #Neoepitopes
			training_fold_full1 <- training_fold_full[,-1]
			training_fold_full1[] <- lapply(training_fold_full1, as.numeric)
			training_fold_C <-training_fold_full1[training_fold_full1$stat==1]
			training_fold_H <-training_fold_full1[training_fold_full1$stat==0]
			cancerCount <- length(training_fold_full1$stat[training_fold_full1$stat==1])
			healthyCount <- length(training_fold_full1$stat[training_fold_full1$stat==0])
			neoepitopes <- data.frame(colSums(training_fold_full1[training_fold_full1$stat==1]), colSums(training_fold_full1[training_fold_full1$stat==0]), colSums(training_fold_C!=0), colSums(training_fold_H!=0))
			colnames(neoepitopes) <- c("cancerLogCount","healthyLogCount","cancerCount","healthyCount")
			neoepitopes <- neoepitopes[-1,]
			neoepitopes$cancerMeans <- neoepitopes$cancerLogCount/cancerCount
			neoepitopes$healthyMeans <- neoepitopes$healthyLogCount/healthyCount
			neoepitopesFs <- neoepitopes[(neoepitopes$cancerMeans/neoepitopes$healthyMeans<0.5 & neoepitopes$healthyCount>=minDonors) | ( neoepitopes$cancerMeans/neoepitopes$healthyMeans>=2 & neoepitopes$cancerCount>=minDonors),]
			keep_cols <- c("sample_id","stat",rownames(neoepitopesFs))
#			print(length(keep_cols))
			training_foldN <- training_fold_full[,..keep_cols]
			test_foldN <- test_fold_full[,..keep_cols]

			tr2 <- inner_join(training_foldT,training_foldN, by=c("sample_id","stat"))
			tr3 <- inner_join(tr2,training_foldC, by=c("sample_id","stat"))
#			training_fold <- setDT(tr2)[, .SD[!all(.SD[, -1, with = F] == 0)], by = sample_id]
#			training_fold <- tr2[rowSums(tr2[,-c(1,2)] == 0) != ncol(tr2)-2,]

# Finally merge the three biomarkers
            te2 <- inner_join(test_foldT,test_foldN, by=c("sample_id","stat"))
    		te3 <- inner_join(te2,test_foldC, by=c("sample_id","stat"))
            training_fold <- tr3[,-1]
			test_fold <- te3[,-1]
		  # now apply (training_fold) on all the classifers and test predictions on test_fold
		  accList <- vector()
		  aucList <- vector()
		  F1List <- vector()
		  sensivList <- vector()
		  specifList <- vector()
	# Apply each method seperately and create the result file and AUROC curve for each method
		  for (m in methods) {
		      print(m)
			  # SVM/RF classifier
			  if(m %in% c("SVM","RF")){
				#setting Y_label as factor
				training_fold$stat <- as.factor(training_fold$stat)
				test_fold$stat <- as.factor(test_fold$stat)
				if(m == "SVM"){
				  classifier <- svm(formula = stat ~ ., data = training_fold,  kernel = "linear", cost = 2, probability = TRUE, scale = FALSE)
  				  y_pred <- predict(classifier, newdata = test_fold[,-c("stat")])
  				  y_prob <- predict(classifier, newdata = (test_fold[,-c("stat")]), probability = T)
				} else {
				  classifier <- randomForest(formula = stat ~ ., data = training_fold, nodesize=2)
  				  y_pred <- predict(classifier, newdata = test_fold[,-c("stat")])
  				  y_prob1 <- predict(classifier, newdata = (test_fold[,-c("stat")]), type = "prob")
				  y_prob <- y_prob1[,1]
				}
				# next step in the loop, we calculate the predictions and cm and we equate the accuracy
				# note we are training on training_fold and testing its accuracy on the test_fold
				#print(results)
				#predictors(results)
			  } else if(m %in% c("Ridge","Lasso")){
				if(m == "Ridge"){
				  alphaVal <- 0
				  bLambda <- 0.1 #0.05 for TEs+wcox; 0.6 for neo+log2fold
				} else {
				  alphaVal <- 1
				  bLambda <- 0.03 #0.02 for TEs+wcox; 0.04 for neo+log2fold
				}
#				lm_cv <- cv.glmnet(as.matrix(training_fold[,-c("stat")]), as.matrix(training_fold[,"stat"]), alpha = alphaVal, lambda = lambda_seq, family="binomial")
#				bLambda <- lm_cv$lambda.min
#				print(bLambda)
				best_lm <- glmnet(as.matrix(training_fold[,-c("stat")]), as.matrix(training_fold[,"stat"]), alpha = alphaVal, lambda = bLambda, family="binomial")
				y_pred <- predict(best_lm, s=bLambda, newx = as.matrix(test_fold[,-c("stat")]), type="class")
				y_prob <- predict(best_lm, s=bLambda, newx = as.matrix(test_fold[,-c("stat")]), type = "response")
			  } else if(m == "NN"){
				classifier=neuralnet(stat == 1 ~.,data=training_fold, hidden=layersNN, err.fct = "ce", act.fct = "logistic", linear.output = FALSE)
				preds=compute(classifier,test_fold[,-c("stat")])
				y_prob <- preds$net.result
				# Converting probabilities into binary classes setting threshold level 0.5
				y_pred <- ifelse(y_prob>0.5, 1, 0)
			  }
			  # calculate prediction accuracy for each fold
			  cmCaret <- caret::confusionMatrix(as.factor(y_pred),as.factor(test_fold$stat))
			  sensivList[[m]] <- cmCaret$byClass['Sensitivity']
			  specifList[[m]] <- cmCaret$byClass['Specificity']
			  F1List[[m]] <- cmCaret$byClass['F1']
	#		  cm <- table(as.factor(y_pred),as.factor(test_fold$stat))
	#		  accuracy <- (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
			  accList[[m]] <- cmCaret$byClass['Balanced Accuracy']
			  rocTest <- roc(as.factor(test_fold$stat), as.numeric(y_prob))
			  if(m=="RF"){
				  pC<-((y-1)*5)+fldc
				  print(pC)
#				  rocTest2 <- roc(as.factor(test_fold$stat), as.numeric(y_prob), smooth=TRUE)
				  if(pC == 1){
				  	rocList <<-list(rocTest)
				  }else{
			  		  rocList <<-c(rocList,list(rocTest))
				  }
			  }
			  aucList[[m]] <- as.numeric(auc(rocTest))
		  }
		  return(list(sensivList, specifList, F1List, accList, aucList))
		})

		cv_dfse <- list()
		cv_dfsp <- list()
		cv_dfF1 <- list()
		cv_dfac <- list()
		cv_dfauc <- list()
	    nti=0
		for(fold in cv)
		{
			nti=nti+1
			#same order as return()
			cv_dfse[[paste0("Fold",nti)]] <- fold[[1]]
			cv_dfsp[[paste0("Fold",nti)]] <- fold[[2]]
			cv_dfF1[[paste0("Fold",nti)]] <- fold[[3]]
			cv_dfac[[paste0("Fold",nti)]] <- fold[[4]]
			cv_dfauc[[paste0("Fold",nti)]] <- fold[[5]]
		}
		colnames(mat_sensiv) <- c(names(cv_dfse),"Average")
		mat_sensiv <- data.frame(cv_dfse) %>% transform(Average=rowSums(.)/K) %>% rbind(.,mat_sensiv)
		colnames(mat_specif) <- c(names(cv_dfsp),"Average")
		mat_specif <- data.frame(cv_dfsp) %>% transform(Average=rowSums(.)/K) %>% rbind(.,mat_specif)
		colnames(mat_F1) <- c(names(cv_dfF1),"Average")
		mat_F1 <- data.frame(cv_dfF1) %>% transform(Average=rowSums(.)/K) %>% rbind(.,mat_F1)
		colnames(mat_accuracy) <- c(names(cv_dfac),"Average")
		mat_accuracy <- data.frame(cv_dfac) %>% transform(Average=rowSums(.)/K) %>% rbind(.,mat_accuracy)
		colnames(mat_auroc) <- c(names(cv_dfauc),"Average")
		mat_auroc <- data.frame(cv_dfauc) %>% transform(Average=rowSums(.)/K) %>% rbind(.,mat_auroc)
	}
	RocC<-nRepeats*K
#	names(rocList)<-c(1:RocC)
#	colStr <- rep("grey", RocC)	
#	ltStr <- rep("dotted",RocC)
	alStr <- rep(0.09,RocC)
	gROCplots[[x1]] <- ggroc(rocList, aes="alpha", color = "#111111", size=0.5) + theme_classic() + ggtitle(paste0(dataset[x1]," ROC curve")) + 
	scale_alpha_manual(values = alStr, guide="none") + geom_smooth(aes(alpha = NULL, group=NULL), se=FALSE ,color="#448877", lty=2, size=0.5) +
	geom_vline(xintercept = 0.8, size = 0.2)
	#+ geom_smooth(aes(alpha = NULL), fill="#44887711" ,color="#448877", lty=2, size=0.4)
	
	fwrite(data.table::data.table(mat_sensiv, keep.rownames=TRUE), file = paste0(ds,"merged3SML_sensitivity",dataset[x1],".tsv"), row.names = F, col.names = T, quote = F, sep = '\t')
	fwrite(data.table::data.table(mat_specif, keep.rownames=TRUE), file = paste0(ds,"merged3SML_specificity",dataset[x1],".tsv"), row.names = F, col.names = T, quote = F, sep = '\t')
	fwrite(data.table::data.table(mat_accuracy, keep.rownames=TRUE), file = paste0(ds,"merged3SML_accuracy",dataset[x1],".tsv"), row.names = F, col.names = T, quote = F, sep = '\t')
	fwrite(data.table::data.table(mat_F1, keep.rownames=TRUE), file = paste0(ds,"merged3SML_F1",dataset[x1],".tsv"), row.names = F, col.names = T, quote = F, sep = '\t')
	fwrite(data.table::data.table(mat_auroc, keep.rownames=TRUE), file = paste0(ds,"merged3SML_AUROC",dataset[x1],".tsv"), row.names = F, col.names = T, quote = F, sep = '\t')
}
pdf(paste0(ds,"merged3_SML_ROCcurves.pdf"), onefile=TRUE, width=4, height=3.5)
	gROCplots
dev.off()
