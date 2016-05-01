library(sqldf)
library(stringr)
library(igraph)
library(entropy)
library(tools)

CONSTANTS<-list(
    uniProt2HGNCDataURL = 'http://www.genenames.org/cgi-bin/download?col=gd_hgnc_id&col=gd_app_sym&col=md_prot_id&status=Approved&status_opt=2&where=&order_by=gd_hgnc_id&format=text&limit=&submit=submit',
    reactome=list(
        pathways=list(fileURL='http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt',
                      downloadedFilename = 'UniProt2Reactome_All_Levels.txt',
                      speciesMap=list(HomoSapiens='Homo sapiens')
        ),
        interactions=list(
            HomoSapiens=list(
                mitabURL='http://www.reactome.org/download/current/homo_sapiens.mitab.interactions.txt.gz',
                downloadedFilename = 'homo_sapiens.mitab.interactions.txt.gz'
            )
        )
    )
)

getUniProtToHGNCSymbolMapping <- function(){
    us <- read.csv(
        url(CONSTANTS$uniProt2HGNCDataURL)
        ,sep="\t",header=FALSE,skip=1,stringsAsFactors = FALSE);
    us<-subset(us[,c(2,3)],V3 != '')
    colnames(us) <- c('HGNCSymbol','uniProtID')
    return(us)
}

createIGraphObject <- function(geneInteractionList){
    g<-graph.edgelist(as.matrix(geneInteractionList))
    vmap<-data.frame(indegree=degree(g,mode = "in"),outdegree=degree(g,mode = "out"))
    row.names(vmap)=V(g)$name
    return(list(igraph_object=g,vertex_map=vmap))
}

downloadReactomePathways <- function(data_folder,species){
    dFilePath <- file.path(data_folder,CONSTANTS$reactome$pathways$downloadedFilename)
    download.file(CONSTANTS$reactome$pathways$fileURL,dFilePath);
    pdata=read.csv(dFilePath,sep='\t',header = FALSE,skip = 1,stringsAsFactors = FALSE)
    hpdata<-subset(pdata,pdata[,6]==CONSTANTS$reactome$pathways$speciesMap[species])[,c(1,2,4)]
    colnames(hpdata)<-c("uniProtID","reactomeID","reactomeDescription")
    return(hpdata)
}

downloadReactomeInteractions <- function(data_folder,species){
    dFilePath <- file.path(data_folder,CONSTANTS$reactome$interactions$`species`$downloadedFilename)
    download.file(CONSTANTS$reactome$interactions[species]$mitabURL,dFilePath);
    return(read.csv(gzfile(dFilePath),header=FALSE,sep="\t",skip=1,stringsAsFactors= FALSE));
}

loadPPIGraphReactome<-function(mitab,usmap){
    idata<-mitab[,c(1,2)]
    idata[,1]<-gsub('uniprotkb:','',gsub('-[0-9]','',idata[,1]));
    idata[,2]<-gsub('uniprotkb:','',gsub('-[0-9]','',idata[,2]));
    colnames(idata) <- c("aup","bup")
    idata<-subset(idata,aup != bup);
    mapgraphA = merge(x=idata,y=usmap,by.x="aup",by.y="uniProtID")[,c(3,2)]
    mapgraph = merge(x=mapgraphA,y=usmap,by.x="bup",by.y="uniProtID")[,c(2,3)]
    colnames(mapgraph)=c("A","B")
    return(createIGraphObject(subset(unique(mapgraph),A!=B)))
}

loadPPIGraphIREF<-function(mitab){
    factori <- sapply(mitab, is.factor)
    mitab[factori] <- lapply(mitab[factori], as.character)
    idata<-unique(
        subset(mitab[,c('aliasA','aliasB')],
               (str_count(aliasA,'hgnc') == 1 & str_count(aliasB,'hgnc') == 1)))
    idata$aliasA<-gsub('hgnc:','',
                       grep('hgnc:',
                            unlist(strsplit(idata$aliasA,'\\|')),value = TRUE))
    idata$aliasB<-gsub('hgnc:','',
                       grep('hgnc:',
                            unlist(strsplit(idata$aliasB,'\\|')),value = TRUE))
    idata<-unique(idata)
    colnames(idata)<-c("A","B")
    idataE<-unique(
        subset(mitab[,c('aliasA','aliasB')],
               (str_count(aliasA,'hgnc') > 1 & str_count(aliasB,'hgnc') > 1)))
    list_function <-
        function(cvec){return(gsub('hgnc:','',grep('hgnc:',cvec,value = TRUE)))}
    idataE$aliasA <-lapply(str_split(idataE$aliasA,'\\|'),list_function)
    idataE$aliasB <-lapply(str_split(idataE$aliasB,'\\|'),list_function)
    idataES<-do.call(rbind,apply(idataE,1,
                                 function(item){
                                     return(merge(x=unlist(item[1]),
                                                  y=unlist(item[2]),
                                                  all=TRUE,
                                                  stringsAsFactors=FALSE))}))
    colnames(idataES)<-c("A","B")
    return(createIGraphObject(subset(unique(rbind(idata,idataES)),A != B)))
}

loadPathwayDataReactome <- function(pdata,usmap){
    memberships=unique(merge(pdata,usmap,by = "uniProtID")[,c(4,2)])
    pathwayMembership=split(memberships,memberships$reactomeID)
    pathways<-unique(pdata[,c(2,3)])
    row.names(pathways)=pathways$reactomeID
    pathways=pathways[,c("reactomeDescription"),drop=FALSE]
    listFunction <-function(df){
        return (df[,c("HGNCSymbol")])
    }
    pathwayMembership=lapply(pathwayMembership, listFunction)
    pathways=pathways[names(pathwayMembership),,drop=FALSE]
    return(list(pathwayInfo=pathways,
                pathwayMembership=pathwayMembership))
}

probeCombinerMean <- function(pls,expdata){
    probes = unlist(strsplit(pls,','))
    if(length(probes) > 1){
        return(colMeans(expdata[probes,]))
    }else{
        return(expdata[probes,])
    }
}

prepareExpressionData <- function(eset,gpl,combinerFunction){
    expdata = exprs(eset)
    dataRowNames<-data.frame(probeName=row.names(expdata))
    hgnc2probe<-Table(gpl)[,c("ID","Gene Symbol")]
    names(hgnc2probe)<-c("probeName","HGNCSymbol")
    selectedProbes<-subset(hgnc2probe,
                           (grepl("[0-9]+_at", probeName) |
                                grepl("[0-9]+_a_at", probeName) ) & !grepl(" ",HGNCSymbol))
    selectedProbes<-merge(x=selectedProbes,y=dataRowNames,by="probeName")[,c(1,2)]
    hgncProbeList =  sqldf("select a.HGNCSymbol, group_concat(a.probeName) as
                         probes from hgnc2probe a inner join selectedProbes b
                         on b.probeName=a.probeName where a.HGNCSymbol != ''
                         and a.HGNCSymbol is not null and a.probeName is not null group by a.HGNCSymbol")
    data = t(sapply(hgncProbeList$probes,combinerFunction,expdata))
    row.names(data) = hgncProbeList$HGNCSymbol
    return(data)
}

logPairedTTestFunction <- function(nData,dData,hgncSymbol,extraArgs){
    if(is.null(nData) && is.null(dData) && is.null(hgncSymbol)){
        return(data.frame(foldChange=numeric(0),pValue=numeric(0),nMean=numeric(0),dMean=numeric(0)))
    }
    nDataLog = log(1+nData)
    dDataLog = log(1+dData)
    result = t.test(dDataLog,nDataLog,paired=TRUE)
    fc=exp(result$estimate)
    fc[fc<1]=1/fc[fc<1]
    return(data.frame(foldChange=fc,pValue=result$p.value,row.names = c(hgncSymbol)))
}

foldChangeFunction <- function(nData,dData,hgncSymbol,extraArgs){
    if(is.null(nData) && is.null(dData) && is.null(hgncSymbol)){
        return(data.frame(foldChange=numeric(0)))
    }
    if(is.null(extraArgs$exponent)){
        exponent=1
    }else{
        exponent=extraArgs$exponent
    }
    if(exponent==1){
        fc=dData/nData;
    }else{
        fc=dData-nData;
    }
    fc=exponent^fc
    fc[fc<1]=1/fc[fc<1]
    fc[is.na(fc)]=1
    return(data.frame(foldChange=fc,row.names = c(hgncSymbol)))
}

runTestOnData <- function (data,normalSampleIndexes,diseaseSampleIndexes,testFunction,testFunctionExtraArgs){
    hgncSymbols = row.names(data)
    normalData = t(data[,normalSampleIndexes])
    diseaseData = t(data[,diseaseSampleIndexes])
    result = testFunction(NULL,NULL,NULL);
    for(i in 1:length(hgncSymbols)){
        result=rbind(result,testFunction(normalData[,i],diseaseData[,i],hgncSymbols[i],testFunctionExtraArgs))
    }
    row.names(result)=hgncSymbols
    return(result)
}

outdegreeNormalizedFCScoreFunction <- function(experimentalData,PPIGraph,extraArgs){
    npv=PPIGraph$vertex_map$outdegree/sum(PPIGraph$vertex_map$outdegree)
    dpv=PPIGraph$vertex_map$outdegree*experimentalData$foldChange
    dpv=dpv/sum(dpv)
    result = data.frame(normal=npv,disease=dpv)
    row.names(result) = row.names(experimentalData)
    return(result)
}

outdegreeNormalizedFCOneMinusPScoreFunction <- function(experimentalData,PPIGraph,extraArgs){
    npv=PPIGraph$vertex_map$outdegree/sum(PPIGraph$vertex_map$outdegree)
    dpv=PPIGraph$vertex_map$outdegree*experimentalData$foldChange*(1-experimentalData$pValue)
    dpv=dpv/sum(dpv)
    result = data.frame(normal=npv,disease=dpv)
    row.names(result) = row.names(experimentalData)
    return(result)
}

computePersonalizationVectors <- function (experimentalData,PPIGraph,defaults,scoreFunction,scoreFunctionExtraArgs){
    common=intersect(row.names(experimentalData),row.names(PPIGraph$vertex_map))
    absent=setdiff(row.names(PPIGraph$vertex_map),row.names(experimentalData))
    if(ncol(experimentalData)==1){
        experimentalDataSubset=experimentalData[common,1,drop=FALSE]
        experimentalDataSubset[absent,1]=defaults
    }else{
        experimentalDataSubset=experimentalData[common,]
        experimentalDataSubset[absent,]=defaults
    }
    experimentalDataSubset=experimentalDataSubset[row.names(PPIGraph$vertex_map),,drop=FALSE]
    result=scoreFunction(experimentalDataSubset,PPIGraph,scoreFunctionExtraArgs)
    row.names(result)=row.names(experimentalDataSubset)
    return(result)
}

computePageRanks <- function(personalizationVectors,PPIGraph,alpha){
    pvs=personalizationVectors[row.names(PPIGraph$vertex_map),]
    npr=page.rank(PPIGraph$igraph_object,personalized = pvs$normal,damping = alpha)
    ppr=page.rank(PPIGraph$igraph_object,personalized = pvs$disease,damping = alpha)
    result=data.frame(normalPageRank=npr$vector,diseasePageRank=ppr$vector)
    return(result)
}

meanAbsoluteDeviationDistanceFunction <- function(normalPageRanks,diseasePageRanks,extraArgs){
    return(data.frame(meanAbsoluteDeviation=sum(abs(diseasePageRanks-normalPageRanks))/length(normalPageRanks)))
}

KLDivergenceDistanceFunction<- function(normalPageRanks,diseasePageRanks,extraArgs){
    return(data.frame(KLDivergence=KL.plugin(normalPageRanks,diseasePageRanks)))
}

computePathwayScores <- function(pageRanks,pathwayMembership,distanceFunction,distanceFunctionExtraArgs){
    listFunction <- function(item,pageranks,analysis,analysisparam){
        return(data.frame(analysis(pageranks[item,"normalPageRank"],pageranks[item,"diseasePageRank"],analysisparam),stringsAsFactors = FALSE))
    }
    return(do.call("rbind",lapply(pathwayMembership,listFunction,pageranks=pageRanks,analysis=distanceFunction,analysisparam=distanceFunctionExtraArgs)))
}
