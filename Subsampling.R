#### Subsampling Analysis
library(utils)
library(tidyverse)
library(limma)
library(progress)
## import spectronaut report file
data <- read.csv(file="eBayes_52Sup_directDIA_9rep.csv",header = TRUE)
## set temperature
temp <- "52"
## rename colnames
colnames(data) <- c("Protein.IDs","Gene.names",paste0("NorIntensity","_",temp,c(rep("_Staurosporine_",9),rep("_Vehicle_",9)),c(paste0("Rep",1:9))),"Log2FC_9rep")
data <- data %>% arrange(Protein.IDs)

## annotate kinases
data_Kinase.Uniprot <- read.table(file="UniprotID_family-isKinase_organism_isHuman_statis-isReviewed.txt",header = T) 
data$is.Kinase.Uniprot <- "no"
data[which(data$Protein.IDs %in% data_Kinase.Uniprot$Protein.IDs),]$is.Kinase.Uniprot <- "yes"
data_Kinase.KinHub <- read.csv(file="KinHub_Kinase.csv",header = T) 
data$is.Kinase.KinHub <- "no"
data[which(data$Protein.IDs %in% data_Kinase.KinHub$Protein.IDs),]$is.Kinase.KinHub <- "yes"
## combine kinase from two databases
myfunc <- function(x){
  ifelse(x[1]=="no" & x[2]=="no","no","yes")
}
kinase_col <- grep("is.Kinase",colnames(data))
data[,paste0("is.Kinase")] <- apply(data[,kinase_col], 1, myfunc)
## annotate ATP-binding proteins
data_ATPbinding <- read.table(file="UniprotID_ATPbinding_organism_isHuman_statis-isReviewed.txt",header = T)
data$is.ATPbinding <- "no"
data[which(data$Protein.IDs %in% data_ATPbinding$Protein.IDs),]$is.ATPbinding <- "yes"
## define staurosporine- and vehicle-treated groups
VehicleCols <- paste0("NorIntensity","_",temp,c(rep("_Vehicle_",9)),c(paste0("Rep",1:9)))
DrugCols <- paste0("NorIntensity","_",temp,c(rep("_Staurosporine_",9)),c(paste0("Rep",1:9)))

### 1.  -------------------------------------------------------------
##Sub-sample sets of "r" measurements from the groups of 9 (This yields choose(9,r)*choose(9,r) situations)
r <- 3
## creat a column in the dataframe to store the -Log10P.Value
data[,paste0("Log10P.Value",r)] = rep(0,nrow(data))
conds <- choose(9,r)
tconds <- conds^2
## using for loops
pb <- progress_bar$new(total = conds)
for (i in 1:conds) {
  pb$tick()
  SubVcols <- combn(VehicleCols,r)[,i]
  for (j in 1:conds) {
    SubDcols <- combn(DrugCols,r)[,j] ## a certain subdata was ready for downstream eBayes analysis
    subsample <- dplyr::select(data,Protein.IDs,Gene.names,append(SubVcols,SubDcols))
    # performing eBayes analysis
    ### Use limma for empirical bayes functions
    # set up model.matrix according to experimental design
    ExpDesign <- model.matrix(~ 0+factor(c(rep(1,r),rep(0,r)))) #0=Vehicle; 1=Drug
    colnames(ExpDesign) <- c("V","D") #0=Vehicle; 1=Drug
    fit <- lmFit(subsample[,append(SubDcols,SubVcols)], ExpDesign)
    fit <- eBayes(fit) ## Apply empirical Bayes smoothing to the standard errors.
    contrast <- makeContrasts("D-V", levels=ExpDesign)
    fit2 <- contrasts.fit(fit, contrast)
    fit2 <- eBayes(fit2, trend = TRUE)
    ptopf <- topTableF(fit2, adjust="BH",genelist = subsample[,"Protein.IDs"],number=Inf)
    metanms <- "Protein.IDs"
    # this table reports the moderated F-statistic
    # f.stat <- topTableF(fit2, adjust="BH",genelist = data[,metanms],number=Inf,p.value=1) 
    # View(f.stat)
    
    # this table reports the moderated t-statistic
    tab.ebayes <- topTable(fit2, adjust="BH",genelist = subsample[,metanms], coef=1,p.value=1, number=Inf)
    nms <- grep("logFC", colnames(tab.ebayes)) 
    colnames(tab.ebayes)[nms] <-  "Log2FC"	#change column names
    
    ## This reports the number that p-value<0.001 during the subsampling for each protein.
    # sig_list <- as.character(tab.ebayes$ID[tab.ebayes$P.Value <= 0.001])
    # data[,paste0("Rep",r)] = data[,paste0("Rep",r)] + (data[,"Protein.IDs"] %in% sig_list)
    
    ## This reports the sum of p-value during the subsampling for each protein
    tab.ebayes$Log10P.Value <- -log(tab.ebayes$P.Value,base = 10)
    tab.ebayes <- tab.ebayes %>% arrange(ID)
    data[,paste0("Log10P.Value",r)]=data[,paste0("Log10P.Value",r)] + tab.ebayes$Log10P.Value
  }
}
## Calculate average Log10P.Value
data[,paste0("Log10P.Value",r)]=data[,paste0("Log10P.Value",r)]/tconds

## output tabel
write.csv(file="Subsampling_final_data.csv",data,row.names = FALSE)


### 2. subsampling with an entire for loop -------------------------------------------------------------------------
## Sub-sample sets of "r" measurements from the groups of 9 (This yields choose(9,r)*choose(9,r) situations)
pb <- progress_bar$new(total = 7)
for (r in 3:9) {
## Create a column in the dataframe to store the frequency that a protein is identified as target of the drug.
  pb$tick()
  data[,paste0("Log10P.Value",r)] = rep(0,nrow(data))
  conds <- choose(9,r)
  tconds <- conds^2
## using for loops
 for (i in 1:conds) {
   SubVcols <- combn(VehicleCols,r)[,i]
   for (j in 1:conds) {
     SubDcols <- combn(DrugCols,r)[,j] ##a certain subdata was ready for downstream eBayes analysis
     subsample <- dplyr::select(data,Protein.IDs,Gene.names,append(SubVcols,SubDcols))
     # performing eBayes analysis
     ## Use limma for empirical bayes functions
     # set up model.matrix according to experimental design
     ExpDesign <- model.matrix(~ 0+factor(c(rep(1,r),rep(0,r)))) #0=Vehicle; 1=Drug
     colnames(ExpDesign) <- c("V","D") #0=Vehicle; 1=Drug
     fit <- lmFit(subsample[,append(SubDcols,SubVcols)], ExpDesign)
     fit <- eBayes(fit) ##Apply empirical Bayes smoothing to the standard errors.
     contrast <- makeContrasts("D-V", levels=ExpDesign)
     fit2 <- contrasts.fit(fit, contrast)
     fit2 <- eBayes(fit2, trend = TRUE)
     ptopf <- topTableF(fit2, adjust="BH",genelist = subsample[,"Protein.IDs"],number=Inf)
     ptopf$rank <- 1:length(ptopf$F)
     ptopf$numfp <- ptopf$rank * ptopf$adj.P.Val
     metanms <- "Protein.IDs"
     # this table reports the moderated F-statistic
     # f.stat <- topTableF(fit2, adjust="BH",genelist = data[,metanms],number=Inf,p.value=1) 
     # View(f.stat)
    
     # this table reports the moderated t-statistic
     tab.ebayes <- topTable(fit2, adjust="BH",genelist = subsample[,metanms], coef=1,p.value=1, number=Inf)
     
     ## This reports the sum of p-value during the subsampling for each protein
     tab.ebayes$Log10P.Value <- -log(tab.ebayes$P.Value,base=10)
     # arrange
     tab.ebayes <- arrange(tab.ebayes,ID)
     data[,paste0("Log10P.Value",r)]=data[,paste0("Log10P.Value",r)] + tab.ebayes$Log10P.Value
  }
 }
  ## Calculate average Log10P.Value
  data[,paste0("Log10P.Value",r)]=data[,paste0("Log10P.Value",r)]/tconds
}
## check
which(data$Protein.IDs!=tab.ebayes$ID)
### output
write.csv(file="subsampling_final_data.csv",data,row.names = F)

