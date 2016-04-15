library(sqldf)
library(stringr)
library(igraph)
library(entropy)
library(GEOquery)
library(Biobase)
uniProt2HGNCDataURL = 'http://www.genenames.org/cgi-bin/download?col=gd_hgnc_id&col=gd_app_sym&col=md_prot_id&status=Approved&status_opt=2&where=&order_by=gd_hgnc_id&format=text&limit=&submit=submit'
getUniProtToHGNCSymbolMapping <- function(){
    us <- read.csv(
        url(uniProt2HGNCDataURL)
        ,sep="\t",header=FALSE,skip=1,stringsAsFactors = FALSE);
    us<-subset(us[,c(2,3)],V3 != '')
    colnames(us) <- c('HGNCSymbol','uniProtID')
    return(us)
}

createIGraphObject <- function(geneInteractionList){
    g<-graph.edgelist(as.matrix(geneInteractionList))
    vmap<-as.data.frame(cbind(hgnc=V(g)$name))
    vmap$vertex_id<-as.numeric(row.names(vmap))
    vmap$indegree<-degree(g,mode = "in")
    vmap$outdegree<-degree(g,mode = "out")
    return(list(igraph_object=g,vertex_map=vmap))
}

loadPPIGraphReactome<-function(mitab,usmap=NULL){
    idata<-mitab[,c(1,2)]
    idata[,1]<-gsub('uniprotkb:','',gsub('-[0-9]','',idata[,1]));
    idata[,2]<-gsub('uniprotkb:','',gsub('-[0-9]','',idata[,2]));
    colnames(idata) <- c("aup","bup")
    idata<-subset(idata,aup != bup);
    if(is.null(usmap)){
        usmap<-getUniProtToHGNCSymbolMapping()
    }
    mapgraphA = merge(x=idata,y=usmap,by.x="aup",by.y="uniProtID")[,c(3,2)]
    mapgraph = merge(x=mangraphA,y=usmap,by.x="bup",by.y="uniProtID")[,c(1,3)]
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
    hpdata<-subset(pdata,pdata[,6]=='Homo sapiens')[,c(1,2,4)]
    colnames(hpdata)<-c("uniProtID","reactomeID","reactomeDescription")
    if(is.null(usmap)){
        usmap<-getUniProtToHGNCSymbolMapping()
    }
    memberships=unique(merge(hpdata,usmap,by = "uniProtID")[,c(4,2)])
    pathways<-unique(hpdata[,c(2,3)])
    return(list(pathwayInfo=pathways,
                pathwayMembership=split(memberships,memberships$reactomeID)))
}

loadExpressionDataGEO <- function(eset,gpl,normalIdx,diseaseIdx){
    dataRowNames<-data.frame(probeName=featureNames(featureData(eset)))
    hgnc2probe<-Table(gpl)[,c("ID","Gene Symbol")]
    names(hgnc2probe)<-c("probeName","HGNCSymbol")
    selectedProbes<-subset(hgnc2probe,
                           (grepl("[0-9]+_at", probeName) |
                                grepl("[0-9]+_a_at", probeName) ) & !grepl(" ",HGNCSymbol))
    selectedProbes<-merge(x=selectedProbes,y=dataRowNames,by="probeName")[,c(1,2)]
    expdata = exprs(eset[selectedProbes$probe,])
    hgncProbeList =  sqldf("select a.hgnc, group_concat(a.probe) as
                         probes from hgnc2probe a inner join selectedProbes b
                         on b.probe=a.probe where hgnc != ''
                         and hgnc is not null group by hgnc ")
    listFunction <- function(pls,expdata){
        probes = unlist(strsplit(pls,','))
        if(length(probes) > 1){
            return(colMeans(expdata[probes,]))
        }else{
            return(expdata[probes,])
        }
    }
    data = t(sapply(hgncProbeList$probes,listFunction,expdata))
    row.names(data) = hgncProbeList$hgnc
    normalData = as.data.frame(data[,normalIdx])
    normalData$hgnc = row.names(data)
    diseaseData = as.data.frame(data[,diseaseIdx])
    diseaseData$hgnc = row.names(data)
    return(list(normal=normalData,disease=diseaseData))
}
