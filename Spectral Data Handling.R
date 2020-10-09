# =================== #
# Required libraries
# =================== #
pckg <- c('readxl',         # Geographic data analysis and modelling
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
data1 <- cbind(VISNIR[1:47], VISNIR = I(as.matrix(VISNIR[-c(1:46)])),MIR=I(as.matrix(MIR[-1])))
str(data1)
colnames(data1$VISNIR) <- gsub("X", "", colnames(data1$VISNIR))
colnames(data1$MIR) <- gsub("X", "", colnames(data1$MIR))
wavmir <- as.numeric(c(colnames(data1$MIR)))
wavvnir <- as.numeric(c(colnames(data1$VISNIR)))
