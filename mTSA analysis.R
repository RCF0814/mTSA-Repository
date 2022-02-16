# 5x5 mTSA analysis
library(limma)
library(stringr)
library(tidyverse)
library(drc)
library(progress)
data <- read.csv(file="mTSA_heat52_Staurosporine_directDIA_Report.csv",header=TRUE)
## experiment design
# set temperature or solvent ratio ***
temp <- "52"
# set replication ***
n <- 5
# set the number of concentrations ***
d <- 5
## rename colnames (# set concentration )
colnames(data) <- c("Protein.IDs","Gene.names",paste0("NorIntensity_",temp,"_Staurosporine_",c(rep("000nM",n),rep("002nM",n),rep("020uM",n),rep("040nM",n),rep("800nM",n)),c(paste0("_Rep",1:n))))
data[data=="Filtered"]=NA
data_temp <- na.omit(data)
str(data_temp)
NorIntensity_col <- grep("NorIntensity", colnames(data_temp))
#change factor to numeric 
for (i in min(NorIntensity_col):max(NorIntensity_col)) {
  data_temp[,i] <- as.numeric(as.character(data_temp[,i]))
}
str(data_temp)
## choose the vehicle and highest drug concentration for ebayes analysis
data_sub <- data_temp %>% dplyr::select(Protein.IDs,Gene.names,contains(c("000nM","020uM")))
## set up model.matrix according to experimental design
ExpDesign <- model.matrix(~ 0+factor(c(rep(1,n),rep(0,n)))) #0=Vehicle; 1=Drug
colnames(ExpDesign) <- c("D","V") #0=Vehicle; 1=Drug
## rename colnames
colnames(data_sub) <- c("Protein.IDs","Gene.names",paste0("NorIntensity","_",temp,c(rep("_Vehicle_",n),rep("_Staurosporine_",n)),c(paste0("Rep",1:n))))
## transform
NorIntensity_col_sub <- grep("NorIntensity", colnames(data_sub))
data_sub[,NorIntensity_col_sub] <- log(data_sub[,NorIntensity_col_sub],base = 2)
## calculate CVs
Cal_CV=function(x){sd(x)/mean(x)}
veh_col=grep("Vehicle",colnames(data_sub))
sta_col=grep("Staurosporine",colnames(data_sub))
data_sub$CV_Vehicle=apply(data_sub[,veh_col],1,Cal_CV)
data_sub$CV_Staurosporine=apply(data_sub[,sta_col],1,Cal_CV)

## Use limma for empirical bayes functions
fit <- lmFit(data_sub[,NorIntensity_col_sub], ExpDesign)
fit <- eBayes(fit) ##Apply empirical Bayes smoothing to the standard errors.
contrast <- makeContrasts("D-V", levels=ExpDesign)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2, trend = TRUE)
ptopf <- topTableF(fit2, adjust="BH",genelist = data_sub[,"Protein.IDs"],number=Inf)
ptopf$rank <- 1:length(ptopf$F)
ptopf$numfp <- ptopf$rank * ptopf$adj.P.Val
metanms <- "Protein.IDs"

# this table reports the moderated t-statistic
t.stat <- topTable(fit2, adjust="BH",genelist = data_sub[,metanms], coef=1,p.value=1, number=Inf)
nms <- grep("logFC", colnames(t.stat)) 
colnames(t.stat)[nms] <-  "Log2FC"	#change column names
t.stat$Log10P.Value <- -log(t.stat$P.Value,base=10)
t.stat$Log10adj.P.Val <- -log(t.stat$adj.P.Val,base=10)
t.statNMS <-t.stat[c("ID","Log2FC","Log10P.Value","Log10adj.P.Val","P.Value","adj.P.Val","AveExpr","t","B") ]   
data_eBayes <- t.statNMS
# rename ID
ID_col <- grep("ID",colnames(data_eBayes))
colnames(data_eBayes)[ID_col] <- "Protein.IDs"
## combine original data with eBayes result data
data_sub <- data_sub %>% left_join(data_eBayes,by="Protein.IDs") 

## import database of kinase and ATP-binding proteins
data_Kinase.KinHub <- read.csv(file="KinHub_Kinase.csv",header = T)
data_Kinase.Uniprot <- read.table(file="UniprotID_family-isKinase_organism_isHuman_statis-isReviewed.txt",head=T)
data_ATPbinding <- read.table(file="UniprotID_ATPbinding_organism_isHuman_statis-isReviewed.txt",header = T)
data_sub$is.Kinase.KinHub <- "no"
data_sub$is.ATPbinding <- "no"
data_sub$is.Kinase.Uniprot <- "no"
data_sub[which(data_sub$Protein.IDs %in% data_Kinase.KinHub$Protein.IDs),]$is.Kinase.KinHub <- "yes"
table(data_sub$is.Kinase.KinHub)
data_sub[which(data_sub$Protein.IDs %in% data_Kinase.Uniprot$Protein.IDs),]$is.Kinase.Uniprot <- "yes"
table(data_sub$is.Kinase.Uniprot)
## combine KinHub database and Uniprot database
myfunc <- function(x){
  ifelse(x[1]=="no" & x[2]=="no","no","yes")
}
kinase_col <- grep("is.Kinase",colnames(data_sub))
data_sub[,paste0("is.Kinase")] <- apply(data_sub[,kinase_col], 1, myfunc)
table(data_sub$is.Kinase)
data_sub[which(data_sub$Protein.IDs%in%data_ATPbinding$Protein.IDs),]$is.ATPbinding <- "yes"
table(data_sub$is.ATPbinding)
## check
data_Plot <- data_sub %>% filter(CV_Vehicle < 0.055 & CV_Staurosporine < 0.055) %>% arrange(is.Kinase.KinHub)

p <- ggplot(data=data_Plot)+
  geom_point(aes(x=Log2FC,y=Log10P.Value,col=is.Kinase),shape=19,alpha=0.7,size=2.4)+
  geom_hline(yintercept = 2.3,linetype="dashed",size=0.5)+
  theme_classic()+
  scale_x_continuous(name = expression(Log["2 "]*staurosporine/vehicle),limits = c(-2.4,2.4))+
  scale_y_continuous(name = expression(-Log["10 "]*p-value),breaks = c(0,6,12),labels = c(0,6,12),limits = c(0,12))+
  scale_color_manual("",values = c("gray18","tomato"),labels=c("Others","Kinase"))+
  theme(legend.position = c(0.13,0.94),legend.background = element_blank(),legend.spacing.x = unit(0.01,"cm"),legend.key.size = unit(0.5,"cm"),legend.spacing.y = unit(0,"cm"),
        axis.title = element_text(family="sans"),axis.text = element_text(size=10))
p

ggsave("SIP8_20uM_0nM.tiff",p,width = 3.5,height = 3.77,dpi=600)
## combine with the original data
data_sub <- data_sub %>% 
  dplyr::select(Protein.IDs,Log2FC,Log10P.Value,is.Kinase.Uniprot,is.Kinase.KinHub,is.Kinase,is.ATPbinding) 
data <- data_temp %>% left_join(data_sub,by="Protein.IDs")
data <- rename(data,Log2FC_20uM_0nM=Log2FC,Log10P.Value_20uM_0nM=Log10P.Value)

## using mean intensity of 5 replicates to evaluate EC50----------------------------------------------------
for (i in c("000nM","002nM","020uM","040nM","800nM")) {
  grp=grep(i,colnames(data))
  data[,paste0("Average_Intensity_",i)]=apply(data[,grp], 1, mean)
}
## normalize to "average_0nM"
for (i in c("000nM","002nM","020uM","040nM","800nM")) {
  data[,paste0("Average_",i,"_000nM")]=data[,paste0("Average_Intensity_",i)]/data[,paste0("Average_Intensity_","000nM")]
}
## Transform dataset for ITDR analysis
TF_data1 <- data %>% filter(Log2FC_20uM_0nM >= 0)
for (x in c("000nM","002nM","020uM","040nM","800nM")) {
  TF_data1[,paste0("Transformed_Intensity_",temp,"_Staurosporine_Mean_",x)]=(TF_data1[,paste0("Average_",x,"_000nM")]-1)/(TF_data1[,paste0("Average_","020uM_000nM")]-1)  
}  
TF_data2 <- data %>% filter(Log2FC_20uM_0nM < 0)
for (x in c("000nM","002nM","020uM","040nM","800nM")) {
  TF_data2[,paste0("Transformed_Intensity_",temp,"_Staurosporine_Mean_",x)]=(TF_data2[,paste0("Average_",x,"_000nM")]-TF_data2[,paste0("Average_","020uM_000nM")])/(1-TF_data2[,paste0("Average_","020uM_000nM")])  
}
data_Final <- rbind(TF_data1,TF_data2)

data_Plot <- data_Final %>% 
  dplyr::select(Protein.IDs,Gene.names,contains(c("mean","Kinase"))) %>% 
  gather(Concentration,Intensity,paste0("Transformed_Intensity_",temp,"_Staurosporine_Mean_000nM"):paste0("Transformed_Intensity_",temp,"_Staurosporine_Mean_800nM")) 

## single plot
par(mfrow=c(1,1))
protein <- "PAK4"
data_protein <- data.frame(conc=c(log10(0),log10(2e-9),log10(2e-5),log10(4e-8),log10(8e-7)),intensity=data_Plot[data_Plot$Gene.names==protein,]$Intensity)
###library limma
protein_drm <- drm(intensity~conc,data=data_protein,fct = LL.4(),logDose = 10)
EC50 <- round(10^6*coef(protein_drm)[[4]],digits = 4)
EC50
pEC50 <- round(-log(coef(protein_drm)[[4]],base = 10),digits = 4)
pEC50
Correlations <- cor(data_protein$intensity, fitted(protein_drm), method = "pearson")
Correlations
CorPvalue <- cor.test(data_protein$intensity, fitted(protein_drm), method = "pearson")$p.value
CorPvalue
par(mar=c(3.5,3.5,2,1))
par(mgp=c(2,0.6,0))
plot(protein_drm,pch=17,broken = T,ylab = "apparent stability to vehicle",xlab= "staurosporine conc. (log M)",col="red",
     cex.lab=1.1,bty="n",lwd=2,family="sans",cex.axis = 0.9,font.lab=1,cex=1.25,xt=c(-10,-9,-8,-7,-6,-5,-4),
     xtlab = c("vehicle","-9","-8","-7","-6","-5","-4"),xttrim = FALSE,yt=c(0,0.2,0.4,0.6,0.8,1.0,1.4),ytlab = c("0","0.2","0.4","0.6","0.8","1","1.4"))
title(main=paste0(protein,"\n","EC50 = ",EC50," \u00b5","M"),adj=0.01,line=0,family="sans",font.main=1,cex.main=1)
## add text mannually 
text(x=-9.4,y=0.8,expression(R^2*"="*" 0.999"))

## extract EC50, Correlations, P-value of the whole dataset--------------------------------------------
data_Final$EC50 <- 0
data_Final$Correlations <- 0
data_Final$CorPavlue <- 0
EC50.col <- grep("EC50",colnames(data_Final))
Correlations.col <- grep("Correlations",colnames(data_Final))
CorPavlue.col <- grep("CorPavlue",colnames(data_Final))
#reset state counter
state <- c(0)
## add progress bar for monitoring
pb <- progress_bar$new(total = length(data_Final$Protein.IDs))

for (proteins in data_Final$Gene.names) {
  pb$tick()
  state <- state + 1
  data_proteins <- data.frame(conc=c(log10(0),log10(2e-9),log10(2e-5),log10(4e-8),log10(8e-7)),intensity=data_Plot[data_Plot$Gene.names==proteins,]$Intensity)
  protein_drm <- tryCatch(
    drm(intensity~conc,data=data_proteins,fct = LL.4(),logDose = 10),
    error=function(e){ print("Correlation LL.4 failed, L.4 backup performed");return(NA)}
  )
  
  # export EC50
  data_Final[state,EC50.col] <- tryCatch(
    round(10^6*coef(protein_drm)[[4]],digits = 4),
    error=function(e){return(NA)}
  )
  # export correlations
  data_Final[state,Correlations.col] <- tryCatch(
    round(cor(data_proteins$intensity, fitted(protein_drm), method = "pearson"),digits = 4),
    error=function(e){return(NA)}
  )
  # export p-value of correlations
  data_Final[state,CorPavlue.col] <- tryCatch(
    round(cor.test(data_proteins$intensity, fitted(protein_drm), method = "pearson")$p.value,digits = 5),
    error=function(e){return(NA)}
  )
}

## adding a cutoff line to select the targets-----------------------------------------
myfun3 <- function(x){
  10^x
}
correlations.col <- grep("Correlations",colnames(data_Final))
data_Final[,paste0("cur.x")] <- apply(as.data.frame(data_Final[,correlations.col]),1,myfun3)
# set criterions
c <- 5.8
x0 <- 6.3095
y0 <- 0.68

cur.f <- function(x){c/(x-x0)+y0}
cur.x.col <- grep("cur.x",colnames(data_Final))
data_Final[,paste0("cur.y")] <- apply(as.data.frame(data_Final[,cur.x.col]), 1, cur.f)
cur.y.col <- grep("cur.y",colnames(data_Final))

myfun4 <- function(x){
  ifelse(x[1]>x[2],"yes","no")
}

Pvalue.col <- grep("Log10P.Value_20uM_0nM",colnames(data_Final))
data_Final[which(data_Final$Log2FC_20uM_0nM > 0 & data_Final$cur.x >= x0+0.1),paste0("is.Target")] <- apply(as.data.frame(data_Final[which(data_Final$Log2FC_20uM_0nM>0 & data_Final$cur.x >= x0+0.1),c(Pvalue.col,cur.y.col)]), 1, myfun4)
data_Final[which(data_Final$Log2FC_20uM_0nM < 0 & data_Final$cur.x >= x0+0.1),paste0("is.Target")] <- apply(as.data.frame(data_Final[which(data_Final$Log2FC_20uM_0nM<0 & data_Final$cur.x >= x0+0.1),c(Pvalue.col,cur.y.col)]), 1, myfun4)
data_Final[data_Final$Log10P.Value_20uM_0nM<=2.33,]$is.Target <- "no"  #  P-value as another filter criteria in Staurosporine analysis
# annotate ADK, CMPK2, and FUK manually.
data_Final[data_Final$Gene.names%in%c("ADK","CMPK2","FUK"),]$is.Kinase.KinHub <- "yes"
table(data_Final[data_Final$is.Target=="yes",]$is.Kinase.KinHub)

# output data
write.csv(file=paste0("data_mean_",temp,"_final_20220216.csv"),data_Final,row.names = FALSE)