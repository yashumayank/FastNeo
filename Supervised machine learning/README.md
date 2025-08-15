Fully automated script to run N iterations of K-fold cross-validation of Random forest, Linear regression (Lasso and Ridge), SVM, and neural networks for classifying cancer patients. The classifiers are trained using expression of neoepitopes, circRNAs and TEs. 

As a result it produces files with F1, Balanced accuracy, Specificity, Sensitivity, and Specificty values for each fold, including averaged values per iteration. It also plot an overlpped ROC curve for all the K fold in N iterations, eg. https://academic.oup.com/bioinformatics/article/41/5/btaf138/8124074#518357705
