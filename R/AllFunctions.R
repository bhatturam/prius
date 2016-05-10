CONSTANTS <- list(
    uniProt2HGNC = list(
                fileURL='http://www.genenames.org/cgi-bin/download?col=gd_hgnc_id&col=gd_app_sym&col=md_prot_id&status=Approved&status_opt=2&where=&order_by=gd_hgnc_id&format=text&limit=&submit=submit',
                downloadedFilename='uniProt2HGNC.csv'
    ),
    reactome = list(
        pathways = list(
            fileURL = 'http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt',
            downloadedFilename = 'UniProt2Reactome_All_Levels.txt',
            speciesMap = list(HomoSapiens = 'Homo sapiens')
        ),
        interactions = list(
            HomoSapiens = list(
                mitabURL = 'http://www.reactome.org/download/current/homo_sapiens.mitab.interactions.txt.gz',
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
#'interactions as well as pathway membership use UniProt ID's. It can reuse
#'downloaded files if the downloaded folder is passed as an argument
#'
#'@param data_folder - The folder where the data resides, or
#'                     where it should be downloaded if it doesnt exist
#'
#'@export
#' @references \url{http://www.genenames.org/cgi-bin/download}
#' @examples
#' \dontrun{
#' data_folder<-tempdir()
#' usmap<-getUniProtToHGNCSymbolMapping(data_folder)
#' }
getUniProtToHGNCSymbolMapping <- function(data_folder) {
    dFilePath <-
        file.path(data_folder,CONSTANTS$uniProt2HGNC$downloadedFilename)
    if(!file.exists(dFilePath)){
        utils::download.file(CONSTANTS$uniProt2HGNC$fileURL,dFilePath);
    }
    us <- utils::read.csv(
        dFilePath
        ,sep = "\t",header = FALSE,skip = 1,stringsAsFactors = FALSE
    );
    us <- subset(us[,c(2,3)],us$V3 != '')
    colnames(us) <- c('HGNCSymbol','uniProtID')
    return(us)
}

#' Create graph object from gene interaction list
#'
#' Convert a data frame containing the gene interaction list into a form usable
#' by the other functions in the library. This is just a wrapper around the
#' igraph graph.edgelist function, which creates and returns an igraph object
#' with a data frame containing the in and outdegrees.
#'
#'
#'
#' @param geneInteractions - a data frame that contains the gene interaction
#'   graph in two columns, the first being the source HGNCSymbol and the second
#'   being the target HGNCSymbol of an interaction. See
#'   \code{\link{getGeneInteractionsFromReactomeMITAB}},
#'   \code{\link{getGeneInteractionsFromIRefMITAB}} for more details
#'
#' @return A list containing igraph_object - the igraph object and vertex_map -
#'   a frame containing node attributes currently indegree and outdegree
#' @export
#'
#' @references \url{http://igraph.org/r/doc/graph_from_edgelist.html}
#' @seealso \code{\link{getGeneInteractionsromReactomeMITAB}},
#'   \code{\link{getGeneInteractionsFromIRefMITAB}}
#'
#' @examples
#' \dontrun{
#'      data_folder=tempdir()
#'      usmap<-getUniProtToHGNCSymbolMapping(data_folder)
#'      mitab <- downloadReactomeInteractionsMITAB(data_folder)
#'      geneInteraction <- getGeneInteractionListFromReactomeMITAB(mitab,usmap)
#'      ppiGraph<-createIGraphObject(geneInteraction)
#' }
#'
createIGraphObject <- function(geneInteractions) {
    g <- igraph::graph.edgelist(as.matrix(geneInteractions))
    vmap <-
        data.frame(indegree = igraph::degree(g,mode = "in"),
                   outdegree = igraph::degree(g,mode = "out"))
    row.names(vmap) <- igraph::V(g)$name
    return(list(igraph_object = g,vertex_map = vmap))
}

#' Download current Pathway Data from Reactome
#'
#' Downloads current pathway data (protein-to-pathway mapping) from reactome,
#' saves it in a local directory and returns a data frame usable by the library.
#' It can reuse downloaded files if the downloaded folder is passed as an
#' argument
#'
#' @param data_folder - The folder where the data resides, or where it should be
#'   downloaded if it doesnt exist
#' @param species - The species being analysed. Currently only supports the
#'   string 'HomoSapiens'
#'
#' @return a data frame containing three columns, uniProtID, reactomeID,
#'   reactomeDescription which stand for the protein identifier, pathway
#'   identifier and the pathway description
#' @export
#'
#' @references \url{http://www.reactome.org/pages/download-data/}
#' @examples
#' \dontrun{
#'      dataFolder <- tempdir()
#'      pathwayData<-downloadReactomePathways(dataFolder)
#' }
downloadReactomePathways <- function(data_folder,species='HomoSapiens') {
    dFilePath <-
        file.path(data_folder,CONSTANTS$reactome$pathways$downloadedFilename)
    if(!file.exists(dFilePath)){
        utils::download.file(CONSTANTS$reactome$pathways$fileURL,dFilePath);
    }
    pdata <- utils::read.csv(
        dFilePath,sep = '\t',header = FALSE,skip = 1,stringsAsFactors = FALSE
    )
    hpdata <-
        subset(pdata,
               pdata[,6] == CONSTANTS$reactome$pathways$speciesMap[species]
               )[,c(1,2,4)]
    colnames(hpdata) <-
        c("uniProtID","reactomeID","reactomeDescription")
    return(hpdata)
}

#' Download PPI Reactions from Reactome in MITAB format
#'
#' Downloads current PPI data (protein-to-protein interactions) in MITAB format
#' from reactome, saves it in a local directory and returns a data frame usable
#' by the library. It can reuse downloaded files if the downloaded folder is
#' passed as an argument
#'
#' @param data_folder - The folder where the data resides, or where it should be
#'   downloaded if it doesnt exist
#' @param species - The species being analysed. Currently only supports the
#'   string 'HomoSapiens'
#'
#' @return a data frame that results when a mitab file is read
#' @export
#' @references \url{http://www.reactome.org/pages/download-data/},
#'   \url{http://www.psidev.info/groups/molecular-interactions}
#' @examples
#' \dontrun{
#'      data_folder <- tempdir()
#'      reactomeMITAB<-downloadReactomeInteractionsMITAB(data_folder)
#' }
downloadReactomeInteractionsMITAB <- function(data_folder,species='HomoSapiens')
{
    dFilePath <-
        file.path(data_folder,
                  CONSTANTS$reactome$interactions$`species`$downloadedFilename)
    if(!file.exists(dFilePath)){
        utils::download.file(CONSTANTS$reactome$interactions[species]$mitabURL,
                             dFilePath);
    }
    return(utils::read.csv(
        gzfile(dFilePath),header = FALSE,sep = "\t",skip = 1,
        stringsAsFactors = FALSE
    ));
}

#' Get a gene interaction data frame from a reactome ppi mitab data frame
#'
#' @param mitab - The mitab data frame containing PPI's
#'                See  \code{\link{downloadReactomeInteractionsMITAB}}
#'
#' @param usmap - A mapping file between UniProt ID's and HGNC Symbols
#'                See \code{\link{getUniProtToHGNCSymbolMapping}}
#'
#' @return  a data frame that contains the gene interaction
#'   graph in two columns, the first being the source HGNCSymbol and the second
#'   being the target HGNCSymbol of an interaction
#' @export
#' @seealso \code{\link{downloadReactomeInteractionsMITAB}}
#'
#' @examples
#' \dontrun{
#'      data_folder=tempdir()
#'      usmap<-getUniProtToHGNCSymbolMapping(data_folder)
#'      mitab <- downloadReactomeInteractionsMITAB(data_folder)
#'      geneInteractions <- getGeneInteractionsFromReactomeMITAB(mitab,usmap)
#' }
getGeneInteractionsFromReactomeMITAB <- function(mitab,usmap) {
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
    geneInteractions <- subset(unique(mapgraph),mapgraph$A != mapgraph$B)
    return(geneInteractions)
}

#' Get a gene interaction data frame from a iRef ppi mitab data frame
#'
#' @param mitab @param mitab - The mitab data frame containing PPI's e.g
#'   Obtained using
#'   \url{http://www.inside-r.org/packages/cran/iRefR/docs/get_irefindex}
#'
#' @return a data frame that contains the gene interaction
#'   graph in two columns, the first being the source HGNCSymbol and the second
#'   being the target HGNCSymbol of an interaction
#' @export
#'
#' @examples
getGeneInteractionsFromIRefMITAB <- function(mitab) {
    factori <- sapply(mitab, is.factor)
    mitab[factori] <- lapply(mitab[factori], as.character)
    idata <- unique(subset(mitab[,c('aliasA','aliasB')],
                           (
                               stringr::str_count(mitab$aliasA,'hgnc') == 1 &
                                   stringr::str_count(mitab$aliasB,'hgnc') == 1
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
                                stringr::str_count(mitab$aliasA,'hgnc') > 1 &
                                    stringr::str_count(mitab$aliasB,'hgnc') > 1
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
    geneInteractions<-unique(rbind(idata,idataES))
    return(subset(geneInteractions,
                  geneInteractions$A != geneInteractionList$B))
}

#' Load pathway data from Reactome
#'
#' @param pdata a data frame containing three columns, uniProtID, reactomeID,
#'   reactomeDescription which stand for the protein identifier, pathway
#'   identifier and the pathway description. See
#'   \code{\link{downloadReactomePathways}},for more details
#' @param usmap A mapping file between UniProt ID's and HGNC Symbols
#'                See \code{\link{getUniProtToHGNCSymbolMapping}}
#'
#' @return A list containing two elements -
#'          pathways- a single column data frame with description
#'          string for each unique pathway id
#'          memberships - a list of character vectors of gene ids
#'          indexed by pathway id, containing the list of hgnc symbols
#'          of genes in each pathway.
#' @export
#'
#' @examples
#' \dontrun{
#'      dataFolder <- tempdir()
#'      rawPathwayData<-downloadReactomePathways(dataFolder)
#'      pdata<-loadPathwayDataReactome(pdata)
#' }
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

#' Combine probes by mean expression value
#'
#' Multiple probes may map to the same HGNC Symbol, and this function
#' that is to be used as input to the \code{\link{prepareExpressionData}} combines
#' these probes by taking the mean of the expression values
#' @param A string containing comma separated probe names
#' @param expdata a matrix of expression values returned from functions such as
#'        \url{'http://svitsrv25.epfl.ch/R-doc/library/Biobase/html/exprs.html'}
#'
#' @return a single value representing the gene expression
#' @export
#'
probeCombinerMean <- function(pls,expdata) {
    probes = unlist(strsplit(pls,','))
    if (length(probes) > 1) {
        return(colMeans(expdata[probes,]))
    }else{
        return(expdata[probes,])
    }
}

#' Select only _at and _a_at probe sets for analysis
#'
#' There are probe sets with different suffixes depending on probe behavior,(See
#' \url{http://www.affymetrix.com/support/help/faqs/hgu133/index.jsp} for more
#' information) and this function that is to be used as input to the
#' \code{\link{prepareExpressionData}} selects only the expression values
#' corresponding probe set suffixes _at and _a_at for analysis
#' @param probeName the name of the probe set
#'
#' @return a binary value indicating if probeName is to be selected
#' @export
#'
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
        and a.HGNCSymbol is not null and
        a.probeName is not null group by a.HGNCSymbol"
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
                    foldChange = numeric(0),
                    pValue = numeric(0)
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
    function (data,normalSampleIndexes,diseaseSampleIndexes,testFunction,
              testFunctionExtraArgs) {
        hgncSymbols = row.names(data)
        normalData = t(data[,normalSampleIndexes])
        diseaseData = t(data[,diseaseSampleIndexes])
        result = testFunction(NULL,NULL,NULL);
        for (i in 1:length(hgncSymbols)) {
            result = rbind(
                result,testFunction(
                    normalData[,i],diseaseData[,i],hgncSymbols[i],
                    testFunctionExtraArgs
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
    function (experimentalData,PPIGraph,defaults,
              scoreFunction,scoreFunctionExtraArgs) {
        common = intersect(row.names(experimentalData),
                           row.names(PPIGraph$vertex_map))
        absent = setdiff(row.names(PPIGraph$vertex_map),
                         row.names(experimentalData))
        if (ncol(experimentalData) == 1) {
            experimentalDataSubset = experimentalData[common,1,drop = FALSE]
            experimentalDataSubset[absent,1] = defaults
        }else{
            experimentalDataSubset = experimentalData[common,]
            experimentalDataSubset[absent,] = defaults
        }
        experimentalDataSubset =
            experimentalDataSubset[row.names(PPIGraph$vertex_map),,drop =
                                                            FALSE]
        result = scoreFunction(experimentalDataSubset,PPIGraph,
                               scoreFunctionExtraArgs)
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
        npr = igraph::page.rank(PPIGraph$igraph_object,
                                personalized = pvs$normal,damping = alpha)
        ppr = igraph::page.rank(PPIGraph$igraph_object,
                                personalized = pvs$disease,damping = alpha)
        result = data.frame(normalPageRank = npr$vector,
                            diseasePageRank = ppr$vector)
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
        return(
            data.frame(
                KLDivergence =
                    entropy::KL.plugin(normalPageRanks,diseasePageRanks)))
    }

#' Title
#'
#' @param pageRanks
#' @param pathwayMembership
#' @param distanceFunction
#' @param distanceFunctionExtraArgs
#'
#' @return
#' @export
#'
#' @examples
computePathwayScores <-
    function(pageRanks,pathwayMembership,distanceFunction,
             distanceFunctionExtraArgs) {
        listFunction <- function(item,pageranks,analysis,analysisparam) {
            return(data.frame(
                analysis(pageranks[item,"normalPageRank"],
                         pageranks[item,"diseasePageRank"],
                         analysisparam),stringsAsFactors = FALSE
            ))
        }
        return(do.call(
            "rbind",lapply(
                pathwayMembership,listFunction,
                pageranks = pageRanks,
                analysis = distanceFunction,
                analysisparam = distanceFunctionExtraArgs
            )
        ))
    }
