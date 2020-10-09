rm(list=ls())
# =================== #
# Required libraries
# =================== #
pckg <- c('readxl',         # Read cel files
          'prospectr',             # Classes and methods for spatial data 
          'corrplot',          # Bindings for the 'Geospatial' data
          'Boruta', # Assessment moel for agriculture soil conditions and crop suitability
          'magrittr',       # A forward-Pipe operator for R
          'aqp',     # Spatial analysis and modelling utilities
          'tidyr',      # Visualizaton methods for raster data
          'hddtools',      # Color schemes for dichromats
          'caret',
          'doParallel',
          'randomForest',
          'factoextra',
          'psych',
          'SuperLearner')   # ColorBrewer palettes

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
matplot(wavmir,t(data1$MIRcr[1:70,]),type='l',ylim=c(0,15),xlab='Wavenumber (cm-1)',ylab='Reflectance')
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
PCAVISNIR <- prcomp(data1$VISNIRkm, scale. = F)  #select only Vis-Nir and MIR bands
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


# ---------- #
# RFE
# ---------- #





# ---------- #
# Boruta
# ---------- #