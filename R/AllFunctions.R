CONSTANTS <- list(
    uniProt2HGNCDataURL = 'http://www.genenames.org/cgi-bin/download?col=gd_hgnc_id&col=gd_app_sym&col=md_prot_id&status=Approved&status_opt=2&where=&order_by=gd_hgnc_id&format=text&limit=&submit=submit',
    reactome = list(
        pathways = list(
            fileURL = 'http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt',
            downloadedFilename = 'UniProt2Reactome_All_Levels.txt',
            speciesMap = list(HomoSapiens = 'Homo sapiens')
        ),
        interactions = list(
            HomoSapiens = list(mitabURL = 'http://www.reactome.org/download/current/homo_sapiens.mitab.interactions.txt.gz',
                               downloadedFilename = 'homo_sapiens.mitab.interactions.txt.gz')
        )
    )
)

#'Get a mapping table between UniProt ID's and HGNC Gene Symbols from
#'genenames.org
#'
#'This function downloads data from genenames.org to create a mapping table
#'between Universal Protein Resource (UniProt) Identifiers to HUGO Gene
#'Nomenclature Committee; gene symbols - which are the default ID's used in this
#'package. This is particularly useful in the case of reactome data, in which
#'interactions as well as pathway membership use UniProt ID's.
#'
#'
#'
#'@export
#' @references \url{http://www.genenames.org/cgi-bin/download}
#' @examples
#' usmap=getUniProtToHGNCSymbolMapping()
getUniProtToHGNCSymbolMapping <- function() {
    us <- utils::read.csv(
        url(CONSTANTS$uniProt2HGNCDataURL)
        ,sep = "\t",header = FALSE,skip = 1,stringsAsFactors = FALSE
    );
    us <- subset(us[,c(2,3)],us$V3 != '')
    colnames(us) <- c('HGNCSymbol','uniProtID')
    return(us)
}

#' Create graph object from gene interaction list
#'
#' Convert a data frame containing the gene interaction list into
#' a form usable by the other functions in the library. This is just
#' a wrapper around the igraph graph.edgelist function, which creates
#' and returns an igraph object with a data frame containing the in and
#' outdegrees.
#'
#'
#'
#' @param geneInteractionList - a data frame that contains the gene interaction
#'                              graph in two columns, the first being
#'                              the source HGNCSymbol and the second being the
#'                              targer HGNCSymbol of an. interaction
#'
#' @return A list containing igraph_object - the igraph object and
#'                           vertex_map - a frame containing node attributes
#'                           currently indegree and outdegree
#' @export
#'
#' @references \url{http://igraph.org/r/doc/graph_from_edgelist.html}
#' @seealso \code{\link{getGeneInteractionListFromReactomeMITAB}}
#'
#' @examples
#' \dontrun{
#'      mitab <- downloadReactomeInteractions(data_folder="/tmp/")
#'      geneInteractionList <- getGeneInteractionListFromReactomeMITAB(mitab)
#'      ppiGraph<-createIGraphObject(geneInteractionList)
#' }
#'
createIGraphObject <- function(geneInteractionList) {
    g <- igraph::graph.edgelist(as.matrix(geneInteractionList))
    vmap <-
        data.frame(indegree = igraph::degree(g,mode = "in"),outdegree = igraph::degree(g,mode = "out"))
    row.names(vmap) <- igraph::V(g)$name
    return(list(igraph_object = g,vertex_map = vmap))
}

#' Title
#'
#' @param data_folder
#' @param species
#'
#' @return
#' @export
#'
#' @examples
downloadReactomePathways <- function(data_folder,species) {
    dFilePath <-
        file.path(data_folder,CONSTANTS$reactome$pathways$downloadedFilename)
    utils::download.file(CONSTANTS$reactome$pathways$fileURL,dFilePath);
    pdata <- utils::read.csv(
        dFilePath,sep = '\t',header = FALSE,skip = 1,stringsAsFactors = FALSE
    )
    hpdata <-
        subset(pdata,pdata[,6] == CONSTANTS$reactome$pathways$speciesMap[species])[,c(1,2,4)]
    colnames(hpdata) <-
        c("uniProtID","reactomeID","reactomeDescription")
    return(hpdata)
}

#' Title
#'
#' @param data_folder
#' @param species
#'
#' @return
#' @export
#'
#' @examples
downloadReactomeInteractionsMITAB <- function(data_folder,species) {
    dFilePath <-
        file.path(data_folder,CONSTANTS$reactome$interactions$`species`$downloadedFilename)
    utils::download.file(CONSTANTS$reactome$interactions[species]$mitabURL,dFilePath);
    return(utils::read.csv(
        gzfile(dFilePath),header = FALSE,sep = "\t",skip = 1,stringsAsFactors = FALSE
    ));
}

#' Title
#'
#' @param mitab
#' @param usmap
#'
#' @return
#' @export
#'
#' @examples
getGeneInteractionListFromReactomeMITAB <- function(mitab,usmap) {
    idata <- mitab[,c(1,2)]
    idata[,1] <- gsub('uniprotkb:','',gsub('-[0-9]','',idata[,1]));
    idata[,2] <- gsub('uniprotkb:','',gsub('-[0-9]','',idata[,2]));
    colnames(idata) <- c("aup","bup")
    idata <- subset(idata,idata$aup != idata$bup);
    mapgraphA <- merge(
        x = idata,y = usmap,by.x = "aup",by.y = "uniProtID"
    )[,c(3,2)]
    mapgraph <- merge(
        x = mapgraphA,y = usmap,by.x = "bup",by.y = "uniProtID"
    )[,c(2,3)]
    colnames(mapgraph) <- c("A","B")
    geneInteractionList <- subset(unique(mapgraph),mapgraph$A != mapgraph$B)
    return(geneInteractionList)
}

#' Title
#'
#' @param mitab
#'
#' @return
#' @export
#'
#' @examples
getGeneInteractionListFromIRefMITAB <- function(mitab) {
    factori <- sapply(mitab, is.factor)
    mitab[factori] <- lapply(mitab[factori], as.character)
    idata <- unique(subset(mitab[,c('aliasA','aliasB')],
                           (
                               stringr::str_count(mitab$aliasA,'hgnc') == 1 & stringr::str_count(mitab$aliasB,'hgnc') == 1
                           )))
    idata$aliasA <- gsub('hgnc:','',
                         grep('hgnc:',
                              unlist(strsplit(
                                  idata$aliasA,'\\|'
                              )),value = TRUE))
    idata$aliasB <- gsub('hgnc:','',
                         grep('hgnc:',
                              unlist(strsplit(
                                  idata$aliasB,'\\|'
                              )),value = TRUE))
    idata <- unique(idata)
    colnames(idata) <- c("A","B")
    idataE <- unique(subset(mitab[,c('aliasA','aliasB')],
                            (
                                stringr::str_count(mitab$aliasA,'hgnc') > 1 & stringr::str_count(mitab$aliasB,'hgnc') > 1
                            )))
    list_function <-
        function(cvec) {
            return(gsub('hgnc:','',grep('hgnc:',cvec,value = TRUE)))
        }
    idataE$aliasA <-
        lapply(stringr::str_split(idataE$aliasA,'\\|'),list_function)
    idataE$aliasB <-
        lapply(stringr::str_split(idataE$aliasB,'\\|'),list_function)
    idataES <- do.call(rbind,apply(idataE,1,
                                   function(item) {
                                       return(merge(
                                           x = unlist(item[1]),
                                           y = unlist(item[2]),
                                           all = TRUE,
                                           stringsAsFactors = FALSE
                                       ))
                                   }))
    colnames(idataES) <- c("A","B")
    geneInteractionList<-unique(rbind(idata,idataES))
    return(subset(geneInteractionList,geneInteractionList$A != geneInteractionList$B))
}

#' Title
#'
#' @param pdata
#' @param usmap
#'
#' @return
#' @export
#'
#' @examples
loadPathwayDataReactome <- function(pdata,usmap) {
    memberships = unique(merge(pdata,usmap,by = "uniProtID")[,c(4,2)])
    pathwayMembership = split(memberships,memberships$reactomeID)
    pathways <- unique(pdata[,c(2,3)])
    row.names(pathways) = pathways$reactomeID
    pathways = pathways[,c("reactomeDescription"),drop = FALSE]
    listFunction <- function(df) {
        return (df[,c("HGNCSymbol")])
    }
    pathwayMembership = lapply(pathwayMembership, listFunction)
    pathways = pathways[names(pathwayMembership),,drop = FALSE]
    return(list(pathwayInfo = pathways,
                pathwayMembership = pathwayMembership))
}

#' Title
#'
#' @param pls
#' @param expdata
#'
#' @return
#' @export
#'
#' @examples
probeCombinerMean <- function(pls,expdata) {
    probes = unlist(strsplit(pls,','))
    if (length(probes) > 1) {
        return(colMeans(expdata[probes,]))
    }else{
        return(expdata[probes,])
    }
}

#' Title
#'
#' @param probeName
#'
#' @return
#' @export
#'
#' @examples
defaultProbeSelector <- function(probeName){
    return( grepl("[0-9]+_at", probeName) |
                grepl("[0-9]+_a_at", probeName))
}

#' Title
#'
#' @param eset
#' @param gpl
#' @param selectorFunction
#' @param combinerFunction
#'
#' @return
#' @export
#'
#' @examples
prepareExpressionData <- function(eset,gpl,selectorFunction,combinerFunction) {
    expdata <- Biobase::exprs(eset)
    dataRowNames <- data.frame(probeName = row.names(expdata))
    hgnc2probe <- GEOquery::Table(gpl)[,c("ID","Gene Symbol")]
    names(hgnc2probe) <- c("probeName","HGNCSymbol")
    selectedProbes <- subset(hgnc2probe,
                             (
                                 selectorFunction(hgnc2probe$probeName)
                             ) & !grepl(" ",hgnc2probe$HGNCSymbol))
    selectedProbes <-
        merge(x = selectedProbes,y = dataRowNames,by = "probeName")[,c(1,2)]
    hgncProbeList <-  sqldf::sqldf(
        "select a.HGNCSymbol, group_concat(a.probeName) as
        probes from hgnc2probe a inner join selectedProbes b
        on b.probeName=a.probeName where a.HGNCSymbol != ''
        and a.HGNCSymbol is not null and a.probeName is not null group by a.HGNCSymbol"
    )
    data <- t(sapply(hgncProbeList$probes,combinerFunction,expdata))
    row.names(data) <- hgncProbeList$HGNCSymbol
    return(data)
}

#' Title
#'
#' @param nData
#' @param dData
#' @param hgncSymbol
#' @param extraArgs
#'
#' @return
#' @export
#'
#' @examples
logPairedTTestFunction <-
    function(nData,dData,hgncSymbol,extraArgs) {
        if (is.null(nData) && is.null(dData) && is.null(hgncSymbol)) {
            return(
                data.frame(
                    foldChange = numeric(0),pValue = numeric(0),nMean = numeric(0),dMean = numeric(0)
                )
            )
        }
        nDataLog = log(1 + nData)
        dDataLog = log(1 + dData)
        result = stats::t.test(dDataLog,nDataLog,paired = TRUE)
        fc = exp(result$estimate)
        fc[fc < 1] = 1 / fc[fc < 1]
        return(data.frame(
            foldChange = fc,pValue = result$p.value,row.names = c(hgncSymbol)
        ))
    }

#' Title
#'
#' @param nData
#' @param dData
#' @param hgncSymbol
#' @param extraArgs
#'
#' @return
#' @export
#'
#' @examples
foldChangeFunction <- function(nData,dData,hgncSymbol,extraArgs) {
    if (is.null(nData) && is.null(dData) && is.null(hgncSymbol)) {
        return(data.frame(foldChange = numeric(0)))
    }
    if (is.null(extraArgs$exponent)) {
        exponent = 1
    }else{
        exponent = extraArgs$exponent
    }
    if (exponent == 1) {
        fc = dData / nData;
    }else{
        fc = dData - nData;
    }
    fc = exponent ^ fc
    fc[fc < 1] = 1 / fc[fc < 1]
    fc[is.na(fc)] = 1
    return(data.frame(foldChange = fc,row.names = c(hgncSymbol)))
}

#' Title
#'
#' @param data
#' @param normalSampleIndexes
#' @param diseaseSampleIndexes
#' @param testFunction
#' @param testFunctionExtraArgs
#'
#' @return
#' @export
#'
#' @examples
runTestOnData <-
    function (data,normalSampleIndexes,diseaseSampleIndexes,testFunction,testFunctionExtraArgs) {
        hgncSymbols = row.names(data)
        normalData = t(data[,normalSampleIndexes])
        diseaseData = t(data[,diseaseSampleIndexes])
        result = testFunction(NULL,NULL,NULL);
        for (i in 1:length(hgncSymbols)) {
            result = rbind(
                result,testFunction(
                    normalData[,i],diseaseData[,i],hgncSymbols[i],testFunctionExtraArgs
                )
            )
        }
        row.names(result) = hgncSymbols
        return(result)
    }

#' Title
#'
#' @param experimentalData
#' @param PPIGraph
#' @param extraArgs
#'
#' @return
#' @export
#'
#' @examples
outdegreeNormalizedFCScoreFunction <-
    function(experimentalData,PPIGraph,extraArgs) {
        npv = PPIGraph$vertex_map$outdegree / sum(PPIGraph$vertex_map$outdegree)
        dpv = PPIGraph$vertex_map$outdegree * experimentalData$foldChange
        dpv = dpv / sum(dpv)
        result = data.frame(normal = npv,disease = dpv)
        row.names(result) = row.names(experimentalData)
        return(result)
    }

#' Title
#'
#' @param experimentalData
#' @param PPIGraph
#' @param extraArgs
#'
#' @return
#' @export
#'
#' @examples
outdegreeNormalizedFCOneMinusPScoreFunction <-
    function(experimentalData,PPIGraph,extraArgs) {
        npv = PPIGraph$vertex_map$outdegree / sum(PPIGraph$vertex_map$outdegree)
        dpv = PPIGraph$vertex_map$outdegree * experimentalData$foldChange *
            (1 - experimentalData$pValue)
        dpv = dpv / sum(dpv)
        result = data.frame(normal = npv,disease = dpv)
        row.names(result) = row.names(experimentalData)
        return(result)
    }

#' Title
#'
#' @param experimentalData
#' @param PPIGraph
#' @param defaults
#' @param scoreFunction
#' @param scoreFunctionExtraArgs
#'
#' @return
#' @export
#'
#' @examples
computePersonalizationVectors <-
    function (experimentalData,PPIGraph,defaults,scoreFunction,scoreFunctionExtraArgs) {
        common = intersect(row.names(experimentalData),row.names(PPIGraph$vertex_map))
        absent = setdiff(row.names(PPIGraph$vertex_map),row.names(experimentalData))
        if (ncol(experimentalData) == 1) {
            experimentalDataSubset = experimentalData[common,1,drop = FALSE]
            experimentalDataSubset[absent,1] = defaults
        }else{
            experimentalDataSubset = experimentalData[common,]
            experimentalDataSubset[absent,] = defaults
        }
        experimentalDataSubset = experimentalDataSubset[row.names(PPIGraph$vertex_map),,drop =
                                                            FALSE]
        result = scoreFunction(experimentalDataSubset,PPIGraph,scoreFunctionExtraArgs)
        row.names(result) = row.names(experimentalDataSubset)
        return(result)
    }

#' Title
#'
#' @param personalizationVectors
#' @param PPIGraph
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples
computePageRanks <-
    function(personalizationVectors,PPIGraph,alpha) {
        pvs = personalizationVectors[row.names(PPIGraph$vertex_map),]
        npr = igraph::page.rank(PPIGraph$igraph_object,personalized = pvs$normal,damping = alpha)
        ppr = igraph::page.rank(PPIGraph$igraph_object,personalized = pvs$disease,damping = alpha)
        result = data.frame(normalPageRank = npr$vector,diseasePageRank = ppr$vector)
        return(result)
    }

#' Title
#'
#' @param normalPageRanks
#' @param diseasePageRanks
#' @param extraArgs
#'
#' @return
#' @export
#'
#' @examples
meanAbsoluteDeviationDistanceFunction <-
    function(normalPageRanks,diseasePageRanks,extraArgs) {
        return(data.frame(meanAbsoluteDeviation = sum(
            abs(diseasePageRanks - normalPageRanks)
        ) / length(normalPageRanks)))
    }

#' Title
#'
#' @param normalPageRanks
#' @param diseasePageRanks
#' @param extraArgs
#'
#' @return
#' @export
#'
#' @examples
KLDivergenceDistanceFunction <-
    function(normalPageRanks,diseasePageRanks,extraArgs) {
        return(data.frame(KLDivergence = entropy::KL.plugin(normalPageRanks,diseasePageRanks)))
    }

computePathwayScores <-
    function(pageRanks,pathwayMembership,distanceFunction,distanceFunctionExtraArgs) {
#' Title
#'
#' @param item
#' @param pageranks
#' @param analysis
#' @param analysisparam
#'
#' @return
#' @export
#'
#' @examples
        listFunction <- function(item,pageranks,analysis,analysisparam) {
            return(data.frame(
                analysis(pageranks[item,"normalPageRank"],pageranks[item,"diseasePageRank"],analysisparam),stringsAsFactors = FALSE
            ))
        }
        return(do.call(
            "rbind",lapply(
                pathwayMembership,listFunction,pageranks = pageRanks,analysis = distanceFunction,analysisparam =
                    distanceFunctionExtraArgs
            )
        ))
    }
