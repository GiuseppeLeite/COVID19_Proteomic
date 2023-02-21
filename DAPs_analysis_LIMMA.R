# UNIFESP, SAO PAULO, BRAZIL
# Giuseppe G. F. Leite, PhD
# Code used to Differential protein abundance analysis
# Ref.: https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf


#library

library("limma")

# Loading Data 

data <- read.delim("Data_To_DAPs.txt", header=TRUE, row.names=1)

# Loading the data with group information 

targets <- read.delim("Target_Global.txt", header=TRUE)

# Matrix design 

f <- factor(targets$Target, levels = unique(targets$Target))
design <- model.matrix(~0 + f)
colnames(design) <- levels(f)

# Apply the intensity values to lmFit

fit <- lmFit(data, design)

# Save your data
write.table(fit, file="fit.txt", sep="\t", quote=FALSE)

# The appropriate contrast matrix can be created by
contrast.matrix <- makeContrasts("D0-HV", "D7-HV", "D30-HV", levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Note that the DAPs of each comparison can now be obtained as follows:

# D0-Control

D0_Control <- topTreat(fit2, coef=1, number=Inf, adjust.method="BH")
write.table(D0_Control, file="D0_Control.txt", sep="\t",
            quote=FALSE)

# D7-Control
D7_Control <- topTreat(fit2, coef=2, number=Inf, adjust.method="BH")
write.table(D7_Control, file="D7_Control.txt", sep="\t",
            quote=FALSE)

# D30_Control
D30_Control <- topTreat(fit2, coef=3, number=Inf, adjust.method="BH")
write.table(D30_Control, file="D30_Control.txt", sep="\t",
            quote=FALSE)

