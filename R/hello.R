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

gds = getGEO('GDS3837',AnnotGPL = TRUE, getGPL = TRUE,destdir = ".")
gpl = getGEO(Meta(gds)$platform,destdir = ".")
eset = GDS2eSet(gds, GPL=gpl,do.log2=FALSE)
pheno = pData(eset)
data = data[,sqldf('select * from pheno order by "disease.state",individual')$sample]
numPairs = dim(data)[2]/2
normalIdex = 1:numPairs
diseaseIdx = (numPairs+1):(2*numPairs)
iref_mitab=get_irefindex(tax_id="9606",data_folder="/tmp/",iref_version = "13.0")
reactome_mitab=reacall <- read.csv(
    url('http://www.reactome.org/download/current/homo_sapiens.mitab.interactions.txt.gz'),
    header=FALSE,
    sep="\t",
    skip=1,
    stringsAsFactors= FALSE);
reactome_pdata = read.csv(url('http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt'),sep='\t',header = FALSE,skip = 1,stringsAsFactors = FALSE)
hello <- function() {
  print("Hello, world!")
}
