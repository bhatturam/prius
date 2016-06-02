library(GEOquery)
library(Biobase)
library(org.Hs.eg.db)
library(TCGAbiolinks)
library(devtools)
install_github("bhatturam/prius")
library(Prius)

#Temporary Folder for experiments
data_folder = tempdir();

#Reactome Interactions and Pathway Data
us=getUniProtToHGNCSymbolMapping(data_folder);
mitab <- downloadReactomeInteractionsMITAB(data_folder=data_folder)
geneInteraction <- getGeneInteractionsFromReactomeMITAB(mitab,us)
ppiGraph<-createIGraphObject(geneInteraction)
reactome_pdata = downloadReactomePathways(data_folder,'HomoSapiens')
reactomePathwayData=loadPathwayDataReactome(reactome_pdata,us)

#Individual Analysis
downloadPathTCGA = file.path(data_folder,"dataBrca");
samplesList=c("TCGA-BH-A0H7-11A-13R-A089-07","TCGA-BH-A0H7-01A-13R-A056-07")
query <- TCGAquery(samples = samplesList,platform = "AgilentG4502A_07_3",level = 3)
TCGAdownload(query, path = downloadPathTCGA,samples = samplesList)
expressionDataTCGA <- data.frame(TCGAprepare(query, dir = downloadPathTCGA,summarizedExperiment = F))
experimentDataTCGA=runTestOnData(expressionDataTCGA,1,2,foldChangeFunction,list(exponent=2))
personalizationVectorsTCGA=computePersonalizationVectors(experimentDataTCGA,ppiGraph,data.frame(foldChange=1),outdegreeNormalizedFCScoreFunction)
pageRanksTCGA=computePageRanks(personalizationVectorsTCGA,ppiGraph,0.7)
pathwayScoresTCGA=computePathwayScores(pageRanksTCGA,reactomePathwayData$pathwayMembership,meanAbsoluteDeviationDistanceFunction,list())
resultTCGA=cbind(reactomePathwayData$pathwayInfo,pathwayScoresTCGA)
resultTCGA=resultTCGA[with(resultTCGA,order(-meanAbsoluteDeviation)),]

#Cohort Analysis
gds = getGEO('GDS3837',AnnotGPL = TRUE, getGPL = TRUE,destdir = data_folder)
gpl = getGEO(Meta(gds)$platform,destdir = data_folder)
eset = GDS2eSet(gds, GPL=gpl,do.log2=FALSE)
numPairs = dim(pData(eset))[1]/2
expressionDataGEO=importExpressionDataGEO(eset,gpl)
experimentDataGEO=runTestOnData(expressionDataGEO,1:60,61:120,logPairedTTestFunction,list())
personalizationVectorsGEO=computePersonalizationVectors(experimentDataGEO,ppiGraph,data.frame(foldChange=1,pValue=1),outdegreeNormalizedFCOneMinusPScoreFunction)
pageRanksGEO=computePageRanks(personalizationVectorsGEO,ppiGraph,0.7)
pathwayScoresGEO=computePathwayScores(pageRanksGEO,reactomePathwayData$pathwayMembership,meanAbsoluteDeviationDistanceFunction,list())
resultGEO=cbind(reactomePathwayData$pathwayInfo,pathwayScoresGEO)
resultGEO=resultGEO[with(resultGEO,order(-meanAbsoluteDeviation)),]
