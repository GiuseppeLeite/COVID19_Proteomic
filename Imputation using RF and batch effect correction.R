# UNIFESP, SAO PAULO, BRAZIL
# Giuseppe G. F. Leite, PhD
# Code used to for imputation and correction of batch effect in proteomics data
# Ref.: A proteomic network approach resolves stage-specific molecular phenotypes in chronic traumatic encephalopathy

#library

library(sva)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(randomForest)
library(tidyverse)

# Read data log2FC, normalized and filtered data
COVID19_Proteomic_Data <- read.delim("COVID19_Proteomic_Data.txt")

# Annotation for batch correction
Annotation_To_Batch <- read.delim("Annotation_To_Batch.txt")



################################################################################
#           Imputation with the "RandomForest" package                         # 
##                 1000 trees in 10 iterations                                ##
################################################################################

# Imputations using patients divided by day and severity

COVID19_Proteomic_Data_IMPUTATION <- rfImpute(as.factor(Clinical_parameter_Severity)~., 
                                                 ntree=1000, 
                                                 iter = 10, 
                                                 data = COVID19_Proteomic_Data)


# Saving the data that has now been imputed

write.table(COVID19_Proteomic_Data_IMPUTATION, 
            file = "Data_After_Imputation.txt", 
            sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE)

################################################################################
#                  Batch effect correction using COMBAT                        #
################################################################################

# Read data
Data_Imputated_RF_Transpor <- read.delim("Data_Imputated_RF_Transpor.txt", 
                                         row.names=1)

# Creating a model based on day and severityand correcting the batch effect

mod = model.matrix(~Clinical_parameter_Severity, data = Annotation_To_Batch)

Batch <- data.frame(annotation, prcomp(Data_Imputated_RF_Transpor, 
                                      scale = TRUE, center = TRUE)$rotation)

#ComBat considering sthe day of sample preparation as a batch effect

Data_combat <- ComBat(Data_Imputated_RF_Transpor, 
                       batch = Batch$Batch, 
                       mod = mod, par.prior = TRUE)


#Saving the data that has now been corrected

write.table(Data_combat, 
            file = "COVID_Imputed_and_Corrected.txt", 
            sep = "\t")


# analisar a normalização dos dados

################################################################################
#                  Principal component analysis (PCA)                          #
################################################################################

# Before the correction

data_pca <- data.frame(Annotation_To_Batch,prcomp(Data_Imputated_RF_Transpor, scale = TRUE, center = TRUE)$rotation)
ggplot(data_pca, aes(x = PC1, y = PC2, shape = Clinical_parameter_Day, color = Batch)) +
  geom_point(size= 7)

# After correction

data_pca <- data.frame(Annotation_To_Batch,prcomp(Data_combat, scale = TRUE, center = TRUE)$rotation)
ggplot(data_pca, aes(x = PC1, y = PC2, shape = Clinical_parameter_Day, color = Batch)) +
  geom_point(size= 7)


