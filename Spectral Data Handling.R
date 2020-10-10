
rm(list=ls())
# =================== #
# Required libraries
# =================== #


# =================== #
# Required libraries
# =================== #
pckg <- c('readxl',         # Read Excel files
          'prospectr',      # Spectral data pre-processing 
          'corrplot',       # Visualization of a correlation matrix   
          'Boruta',         # Feature selection with the Boruta algorithm
          'magrittr',       # A forward-Pipe operator for R
          'aqp',            # Algorithms for quantitative pedology
          'tidyr',          # Tidy Messy Data
          'hddtools',       # Hydrological Data Discovery Tools
          'caret',          # Classification And REgression Training
          'doParallel',     # Foreach parallel adaptator for the 'parallel' package
          'factoextra',     # Extract and visualize the results of multivariate data analysis
          'FactoMineR',     # Multivariate exploratory data analysis and mining
          'dplyr',
          'abind',
          'psych'          # Procedures for psychological, psychometric, and personality research
      )   

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}

lapply(pckg,usePackage)


setwd("G:\\My Drive\\ESALQ_USP\\REMOTE_SENSING\\PROJECT_PERNAMBUCO\\Scripts")


#==================================
## Data loading
#==================================
VISNIR <- read_excel("G:\\My Drive\\ESALQ_USP\\REMOTE_SENSING\\PROJECT_PERNAMBUCO\\Data\\PROFILES_DATA.xlsx",
                     sheet="horizonsVISNIR") %>% data.frame
MIR <- read_excel("G:\\My Drive\\ESALQ_USP\\REMOTE_SENSING\\PROJECT_PERNAMBUCO\\Data\\PROFILES_DATA.xlsx",
                  sheet="horizonsMIR") %>% data.frame
MIR <- MIR[,c(4,47:dim(MIR)[2])]


#XRF <- read_excel("G:\\My Drive\\ESALQ_USP\\REMOTE_SENSING\\PROJECT_PERNAMBUCO\\PROFILES_DATA.xlsx",
#                  sheet="horizonsXRF") %>% data.frame


sites <- read_excel("G:\\My Drive\\ESALQ_USP\\REMOTE_SENSING\\PROJECT_PERNAMBUCO\\Data\\PROFILES_DATA.xlsx",
                    sheet="sites") %>% data.frame


#========================================
## Matrix embedding into data frames
#========================================
names(VISNIR)
data1 <- cbind(VISNIR[1:47], VISNIR = I(as.matrix(VISNIR[-c(1:47)])),MIR=I(as.matrix(MIR[-1])))
str(data1)
colnames(data1$VISNIR) <- gsub("X", "", colnames(data1$VISNIR))
colnames(data1$MIR) <- gsub("X", "", colnames(data1$MIR))
wavmir <- as.numeric(c(colnames(data1$MIR)))
wavvnir <- as.numeric(c(colnames(data1$VISNIR)))



# ========================================== #
# Soil spectral data pre-processing VIS-NIR
# ========================================== #
# 1) Continuum removal
data1$VISNIRcr <- continuumRemoval(data1$VISNIR,wavvnir,type='R')##Adding new columns to the previous data frame
matplot(wavvnir,t(data1$VISNIRcr[1:70,]),type='l',ylim=c(0,1),xlab='Wavenumber (cm-1)',ylab='Reflectance')
#matlines(as.numeric(colnames(data1$VISNIR)), t(data1$VISNIR[1:10, ]))


# 2) Kubelka-Munk function
data1$VISNIRkm <- ((1-data1$VISNIRcr)^2)/(2*data1$VISNIRcr)
matplot(wavvnir,t(data1$VISNIRkm[1:70,]),type='l',ylim=c(0,1),xlab='Wavenumber (cm-1)',ylab='Reflectance')


# 3) Computing absorbance
data1$VISNIRabs <- log10(1/data1$VISNIRkm)
matplot(wavvnir,t(data1$VISNIRabs[1:70,]),type='l',ylim=c(0,20),xlab='Wavenumber (cm-1)',ylab='Absorbance')
#matlines(as.numeric(colnames(data1$VISNIR)), t(data1$VISNIR[1:10, ]))

# ========================================== #
# Soil spectral data pre-processing MIR
# ========================================== #

# 1) Continuum removal
data1$MIRcr <- continuumRemoval(data1$MIR,wavmir,type='R')##Adding new columns to the previous data frame
matplot(wavmir,t(data1$MIRcr[1:70,]),type='l',ylim=c(0,100),xlab='Wavenumber (cm-1)',ylab='Reflectance')
#matlines(as.numeric(colnames(data1$VISNIR)), t(data1$VISNIR[1:10, ]))


# 2) Kubelka-Munk function
data1$MIRkm <- ((1-data1$MIRcr)^2)/(2*data1$MIRcr)
matplot(wavmir,t(data1$MIRkm[1:70,]),type='l',ylim=c(0,15),xlab='Wavenumber (cm-1)',ylab='Reflectance')


# 3) Computing absorbance
data1$MIRabs <- log10(1/data1$MIRkm)
matplot(wavmir,t(data1$MIRabs[1:70,]),type='l',ylim=c(0,15),xlab='Wavenumber (cm-1)',ylab='Absorbance')
#matlines(as.numeric(colnames(data1$VISNIR)), t(data1$VISNIR[1:10, ]))


# =============================== #
# Features selection
# =============================== #

# ---------- #
# PCA
# ---------- #
#PCA of Vis-Nir  with Principal Component Analysis#
PCAVISNIR <- prcomp(data1$VISNIRcr, scale. = F)  #select only Vis-Nir and MIR bands
##Plot PCAs
fviz_eig(PCAVISNIR)   #To select the number of optimal PCs
summary(PCAVISNIR)

PCAVISNIR$center       ##outputs the mean of variables
PCAVISNIR$scale        ##outputs the standard deviation of variables
PCAVISNIR$rotation     ##The rotation measure provides the principal component leading. Each column of rotation matrix contains the principal component loading vector. This is the most important measure we should be interested in.

PCA <- as.data.frame(PCAVISNIR$rotation[,1:3])
write.csv(PCA, "PCA_P1.csv")
fviz_pca_var(PCAVISNIR,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)    # Avoid text overlapping

(corvar <- PCAVISNIR$rotation %*% diag(PCAVISNIR$sdev))# loadings


# -------------------------------------------------------- #
# Data preparation for features selection (RFE and Boruta)
# -------------------------------------------------------- #
VISNIRcr <- matrix(data1$VISNIR,
       nrow=dim(data1$VISNIR)[1],
       ncol=dim(data1$VISNIR)[2]) %>% data.frame
colnames(VISNIRcr) <- wavvnir
data2 <- data.frame(ORDER=as.factor(VISNIR$ORDER),VISNIRcr)


MIRcr <- matrix(data1$MIR,
                nrow=dim(data1$MIR)[1],
                ncol=dim(data1$MIR)[2]) %>% data.frame
colnames(MIRcr) <- wavmir
data3 <- data.frame(ORDER=as.factor(VISNIR$ORDER),MIRcr)
dim(data3)
sapply(data3, function(x) sum(is.na(x)))
data3 <- data3[ , colSums(is.na(data3)) == 0]



# ---------- #
# RFE
# ---------- #
start <- Sys.time()
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)
control2 <- rfeControl(functions=rfFuncs, method="repeatedcv", number=5, repeats=5)
(rfmodel <- rfe(x=data2[,-1], y=data2[,1], sizes=c(1:10), rfeControl=control2))
plot(rfmodel, type=c("g", "o"))
predictors(rfmodel)[1:5]
print(Sys.time() - start)




# ---------- #
# Boruta
# ---------- #
start <- Sys.time()
(bor <- Boruta(ORDER ~ ., data = data3, doTrace = 0, ntree = 500,maxRuns=100,))
plot(bor, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(bor$ImpHistory),function(i)
  bor$ImpHistory[is.finite(bor$ImpHistory[,i]),i])
names(lz) <- colnames(bor$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(bor$ImpHistory), cex.axis = 0.7)
print(Sys.time() - start)

print(bor)
names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])

final.bor<-TentativeRoughFix(bor)
print(final.bor)

getConfirmedFormula(final.bor)
namesdef<-names(final.bor$finalDecision[final.bor$finalDecision %in% c("Confirmed")])
namesdef


data2 <- data2[,namesdef]
names(data2)
dim(data2)

data3 <- data3[,namesdef] 
names(data3)
dim(data3)

# =========================== #
# Outer product VISNIR x MIR
#============================ #

str(data2)
dim(data3)

samples <- 1:70

opa <- c()
for (i in 1:length(samples)) {
  z<- as.matrix(data2[i,],1)%o%as.matrix(data3[i,],1)
  temp<- matrix(z,nrow = ncol(data2),ncol=ncol(data3),byrow = T)
  opa <- abind(opa,temp,along=3) 
}
opa
dim(opa)
# class(opa)
# opa[,,70]

# =========================== #
# Unfold OPA matrix
#============================ #
opafinal <- c()
for (i in 1:length(samples)) {
  temp <- c(opa[,,i])
  opafinal <- rbind(opafinal,temp)
}
class(opafinal)
dim(opafinal)

# ====================================== #
# Hierarchical clustering
# ====================================== #


# ------------------------------------ #
# Option 1 - Hierarchical Clustering
# ------------------------------------ #
dist <- dist(df2[,-1])
hc <- hclust(dist,"ward.D")
plot(hc)

# ------------------------------------------ #
# Option 2 - Hierarchical clustering on PCA
# ------------------------------------------ #

#Computing median by profile ID and then perform HCPCA
df <- cbind(profileID=VISNIR$profileID,data.frame(opafinal)) %>% data.frame()
str(df)
names(df)
df <- as_tibble(df)
df2 <- df %>%
  group_by(profileID ) %>%
  summarise(across(where(is.numeric), median, na.rm= TRUE)) %>% data.frame
dim(df2)

View(df2)

# Compute PCA with ncp = 3
res.pca <- PCA(df2[,-1], ncp = 3, graph = FALSE)
# Compute hierarchical clustering on principal components
res.hcpc <- HCPC(res.pca,nb.clust = -1, graph = FALSE)
fviz_dend(res.hcpc, 
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.8      # Augment the room for labels
)
fviz_cluster(res.hcpc,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
)

 
