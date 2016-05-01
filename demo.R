# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
library(iRefR)
library(GEOquery)
library(Biobase)
library(org.Hs.eg.db)
library(TCGAbiolinks)
samplesList=c("TCGA-BH-A0H7-11A-13R-A089-07","TCGA-BH-A0H7-01A-13R-A056-07")
query <- TCGAquery(samples = samplesList,platform = "AgilentG4502A_07_3",level = 3)
TCGAdownload(query, path = "/tmp/dataBrca",samples = samplesList)
expressionDataTCGA <- data.frame(TCGAprepare(query, dir = "/tmp/dataBrca",summarizedExperiment = F))
gds = getGEO('GDS3837',AnnotGPL = TRUE, getGPL = TRUE,destdir = "/tmp")
gpl = getGEO(Meta(gds)$platform,destdir = "/tmp")
eset = GDS2eSet(gds, GPL=gpl,do.log2=FALSE)
numPairs = dim(pData(eset))[1]/2
expressionDataGEO=prepareExpressionData(eset,gpl,probeCombinerMean)
iref_mitab=get_irefindex(tax_id="9606",data_folder="/tmp/",iref_version = "13.0")
irefPPIGraph = loadPPIGraphIREF(iref_mitab)
us=getUniProtToHGNCSymbolMapping();
reactome_mitab= downloadReactomeInteractions(data_folder="/tmp/")
reactomePPIGraph=createIGraphObject(subset(unique(mapgraph),A!=B))
reactome_pdata = read.csv(url('http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt'),sep='\t',header = FALSE,skip = 1,stringsAsFactors = FALSE)
reactomePathwayData=loadPathwayDataReactome(reactome_pdata,us)
experimentDataGEO=runTestOnData(expressionDataGEO,1:60,61:120,logPairedTTestFunction)
experimentDataTCGA=runTestOnData(tcgaData,1,2,foldChangeFunction,list(exponent=2))
