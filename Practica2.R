
########################
## Original data sets ##
########################

vantVDir <- system.file("extdata/vantVeer", package = "seventyGeneData")
dir(vantVDir)
vantDeVDir <- system.file("extdata/vanDeVijver", package = "seventyGeneData")
dir(vantDeVDir)

###Load the libraries
library(limma)
library(gdata)
library(Biobase)
library(breastCancerNKI)

#################################
### EXPRESSIONSET PREPARATION ###
#################################

##############################
## Feature data preparation ##
##############################

data(nki)
dim(nki)
#24481 genes, 337 pacientes

#resumen de formato
str(nki)
#resumen de datos
show(nki)
#datos de expresión
exprs(nki)[1:5,1:5]

#para cada paciente (sampleName): 21 variables (phenoData)
#datos fenotípicos
head(pData(nki))
#guardamos phenoData en una variable
nkiPheno<-pData(nki)

#para cada gen (featureName): 10 variables (featureData)
#Comprobar featureData
str(featureData(nki))
#nombres genes
featureNames(nki)[1:20]
#guardamos featureData en una variable
nkiAnn <- featureData(nki)
###Turn all annotation information into character
nkiAnn@data <- as.data.frame(apply(nkiAnn@data, 2, as.character),stringsAsFactors=FALSE)

###Read and manipulate the annotation contained in files from Nature and NKI

###Set working directory
setwd(vantVDir)
###Read GBACC information for van't Veer dataset
featAcc <- read.xls("ArrayNomenclature_contig_accession.xls",skip=0,header=TRUE,stringsAsFactors=FALSE)
###Read seq information for van't Veer dataset
featSeq <- read.xls("ProbeSeq.xls",skip=0,header=TRUE,stringsAsFactors=FALSE)
###Read 70-genes signature information for van't Veer dataset
gns231 <- read.xls("415530a-s9.xls",skip=0,header=TRUE,stringsAsFactors=FALSE)
###Remove special characters in the colums header,
###which are due to white spaces present in the Excel files
colnames(gns231) <- gsub("\\.\\.", "", colnames(gns231))
###Remove GO annotation
gns231 <- gns231[, -grep("sp_xref_keyword_list", colnames(gns231))]
###Reorder the genes in decreasing order by absolute correlation
gns231 <- gns231[order(abs(gns231$correlation), decreasing=TRUE),]
###Select the feature identifiers corresponding to the top 231 and 70 genes
gns231$genes231 <- TRUE
gns231$genes70 <- gns231$accession %in% gns231$accession[1:70]
###Merge all information (including 70-gene signature information)
###with the annotation obtained from the breastCancerNKI package
newAnn <- nkiAnn@data
newAnn <- merge(newAnn, featAcc, by.x=1, by.y=1, all=TRUE, sort=FALSE)
newAnn <- merge(newAnn, featSeq, by.x=1, by.y=1, all=TRUE, sort=FALSE)
newAnn <- merge(newAnn, gns231, by.x=1, by.y=1, all=TRUE, sort=FALSE)
###Check the structure of the new annotation data.frame
newAnn <- newAnn[order(newAnn[,1]),]
str(newAnn)

#########################################
## Assemble the vantVeer ExpressionSet ##
#########################################

###Create list of files to be read in
filesVtV <- dir(vantVDir, full.names=TRUE, pattern="^ArrayData")
filesVtV

phenoVtV <- read.xls("415530a-s8.xls",skip=0,header=TRUE,stringsAsFactors=FALSE)
str(phenoVtV)
###Remove the special characters in the colums headers
###due to white spaces present in the Excel file
colnames(phenoVtV) <- gsub("\\.$", "", gsub("\\.$", "", colnames(phenoVtV)))
####Remove columns that do not contain useful information
phenoVtV$SampleName <- paste("Sample", phenoVtV$Sample)
rownames(phenoVtV) <- phenoVtV$SampleName
###Read sample names from the 6 expression data tables
samplesVtV <- lapply(filesVtV,read.table,nrow=1,header=FALSE,sep="\t",stringsAsFactors=FALSE,fill=TRUE,strip.white=TRUE)
###Format the samples strings
samplesVtV <- lapply(samplesVtV, function(x) x[ grep("^Sample", x) ])
headerDesc <- samplesVtV
samplesVtV <- lapply(samplesVtV, function(x) gsub(",.+", "", x) )
###Check sample lables obtained from expression data files
str(samplesVtV)
###Combine the lables in one unique vector
allSamplesVtV <- do.call("c", samplesVtV)
###Compare order the order the samples between the expression data
###and phenotypic information data.frames
if (all(rownames(phenoVtV) %in% allSamplesVtV)) {
  print("All sample names match phenoData")
  if (all(rownames(phenoVtV) == allSamplesVtV)) {
    print("All sample names match phenoData")
  } else {
    print("Sample names from tables and phenoData need reordering")
    phenoVtV <- phenoVtV[order(phenoVtV$SampleName), ]
  }
} else {
  print("Sample names DO NOT match phenoData")
}
###Read expression data from the 4 converted TAB-delimited text files
dataVtV <- lapply(filesVtV, read.table, skip=1, sep="\t", quote="",header=TRUE,
                  row.names=NULL,stringsAsFactors=FALSE, fill=FALSE, strip.white=FALSE)
sapply(dataVtV, dim)
###Extract annotation: note that column headers are slightly different
sapply(dataVtV, function(x) head(colnames(x)) )
sapply(dataVtV, function(x) tail(colnames(x)) )
###Extract the associated annotation
annVtV <- lapply(dataVtV, function(x) x[,c("Systematic.name", "Gene.name")])
annVtV <- lapply(annVtV, function(x) {x[x==""] <- NA ; x })
annVtV <- do.call("cbind", annVtV)
###Check annotation order in all data files
if ( all(apply(annVtV[, seq(1, 8, by=2)], 1, function(x) length(unique(x)) == 1 )) ) {
  print("OK")
  annVtV <- annVtV[,1:2]
} else {
  print("Check annotation")
}
###Define the function
extractColumns <- function(x, pattern, ann) {
  sel <- grep(pattern, colnames(x), value=TRUE)
  x <- x[,sel ]
  rownames(x) <- ann
  x <- x[order(rownames(x)), ]
}
###Extract log ratio data from all the spreadsheets
logRat <- lapply(dataVtV, extractColumns, pattern="Log10\\.ratio", ann=annVtV[,1])
logRat <- do.call("cbind", logRat)
###Assign colnames and reorder the columns
colnames(logRat) <- allSamplesVtVAs
logRat <- logRat[, order(colnames(logRat)),]
###Check order
all(phenoVtV$SampleName == colnames(logRat))
###Extract p-values from all the spreadsheets
pVal <- lapply(dataVtV, extractColumns, pattern="value", ann=annVtV[,1])
pVal <- do.call("cbind", pVal)
###Assign colnames and reorder the columns
colnames(pVal) <- allSamplesVtV
pVal <- pVal[, order(colnames(pVal)),]
###Check order
all(phenoVtV$SampleName == colnames(pVal))
intensity <- lapply(dataVtV, extractColumns, pattern="Intensity", ann=annVtV[,1])
intensity <- do.call("cbind", intensity)
###Assign colnames and reorder the columns
colnames(intensity) <- allSamplesVtV
intensity <- intensity[, order(colnames(intensity)),]
###Check order
all(phenoVtV$SampleName == colnames(intensity))
###Merge annotation objects and check order
annVtV <- merge(annVtV, newAnn, by=1, all=TRUE, sort=TRUE)
rownames(annVtV) <- annVtV[,1]
all(rownames(annVtV) == rownames(logRat))
all(rownames(annVtV) == rownames(pVal))
all(rownames(annVtV) == rownames(intensity))
###Create the new assayData
myAssayData <- assayDataNew(exprs=logRat, pValue=pVal, intensity=intensity)
###Create the new phenoData
myPhenoData <- new("AnnotatedDataFrame", phenoVtV)
###Create the new featureData
myFeatureData <- new("AnnotatedDataFrame", annVtV)
###Create the new experimentData
myExperimentData <- new("MIAME", name = "Marc J Van De Vijver, Hongyue Dai, and Laura J van't Veer",
                        lab = "The Netherland Cancer Institute, Amsterdam, The Netherlands",
                        contact = "Luigi Marchionni <marchion@gmail.com>",
                        title = "Gene expression profiling predicts clinical outcome of breast cancer",
                        abstract = "Breast cancer patients with the same stage of disease can have 
                        markedlydifferent treatment responses and overall outcome. The strongest 
                        predictors for metastases (for example, lymph node status and 
                        histological grade) fail to calssify accurately breast tumours according to 
                        their clinical behaviour. Chemotherapy or hormonal therapy reduces the risk of 
                        distant metastases by approximately one-third; however, 70-80% of patients 
                        receiving this treatment would have survived without it. None of the signatures 
                        of breast cancer gene expression reported to date allow for patient-tailored 
                        therapy strategies. Here we used DNA microarray analysis on primary breast 
                        tumours of 117 young patients, and applied supervised classification to 
                        identify a gene expression signature strongly predictive of a short interval 
                        to distant metastases ('poor prognosis' signature) in patients without
                        tumor cells in local lymph nodes al diagnosis (lymph node negative).
                        In addition, we established a signature that identifies tumours of 
                        BRCA1 carriers. The poor prognosis signature consist of genes regulating 
                        cell cycle, invasion, metastasis and angiogenesis. This gene expression 
                        profile will outperform all currently used clinical parameters in predicting
                        disease outcome. Our findings provide a strategy to select patients who would
                        benefit from adjuvant therapy.",
                        url = "http://www.ncbi.nlm.nih.gov/pubmed/?term=11823860",
                        pubMedIds = "11823860" )
###Create the expression set
vantVeer <- new("ExpressionSet",
                assayData = myAssayData,
                phenoData = myPhenoData,
                featureData = myFeatureData,
                experimentData = myExperimentData)


############################################
## Assemble the vanDeVijver ExpressionSet ##
############################################

###Set working directory
setwd(vantDeVDir)

###Create list of files containing expression
filesVdV <- dir(vantDeVDir, full.names=TRUE, pattern="NKI")
filesVdV

###Read phenotypic information
myFile <- dir(vantDeVDir, full.names=TRUE, pattern="Table1_ClinicalData_Table.xls")
phenoVdV <- read.xls(myFile, skip=2, header=TRUE, stringsAsFactors=FALSE)
###Remove columns that do not contain useful information
phenoVdV <- phenoVdV[ , apply(phenoVdV, 2, function(x) length(unique(x)) > 1)]
phenoVdV$SampleName <- paste("Sample", phenoVdV$SampleID)
rownames(phenoVdV) <- phenoVdV$SampleName
###Read sample names form the expression data spreadsheets
samplesVdV <- lapply(filesVdV, scan, what="character", nlines=1, sep="\t", strip.white=FALSE)
samplesVdV <- lapply(samplesVdV, function(x) x[x!=""])
allSamplesVdV <- do.call("c", samplesVdV)
###Read all data contained in the expression data spreadsheets
dataVdV <- lapply(filesVdV, read.table, header=TRUE, skip=1, sep="\t", quote="",
                  stringsAsFactors=FALSE, fill=TRUE, strip.white=TRUE)
###Extract feature annotation
annVdV <- lapply(dataVdV, function(x) x[,c("Substance", "Gene")])
annVdV <- lapply(annVdV, function(x) {x[x==""] <- NA ; x })
annVdV <- do.call("cbind", annVdV)

###Check annotation order in all data files
if( all(apply(annVdV[, seq(1, 12, by=2)], 1, function(x) length(unique(x)) == 1 )) ) {
  print("OK")
  annVdV <- annVdV[,1:2]
} else {
  print("Check annotation")
}

###Define the function
extractColumns <- function(x, pattern, annVdV) {
  colnames(x) <- gsub("Log\\.Ratio\\.Error", "Error", colnames(x))
  sel <- grep(pattern, colnames(x), value=TRUE)
  x <- x[,sel ]
  rownames(x) <- annVdV
  x <- x[order(rownames(x)), ]
}

###Extract and assemble the log ratio values
logRat <- lapply(dataVdV, extractColumns, pattern="Log\\.Ratio", ann=annVdV[,1])
logRat <- do.call("cbind", logRat)
###Set the column names
colnames(logRat) <- allSamplesVdV
###Check order
all(phenoVdV$SampleName == colnames(logRat))

###Extract log ratio error
logRatError <- lapply(dataVdV, extractColumns, pattern="Error", ann=annVdV[,1])
logRatError <- do.call("cbind", logRatError)
###Set the column names
colnames(logRatError) <- allSamplesVdV
###Check order
all(phenoVdV$SampleName == colnames(logRatError))

###Extract P-value
pVal <- lapply(dataVdV, extractColumns, pattern="alue", ann=annVdV[,1])
pVal <- do.call("cbind", pVal)
###Set the column names
colnames(pVal) <- allSamplesVdV
###Check order
all(phenoVdV$SampleName == colnames(pVal))

###Extract Intensity
intensity <- lapply(dataVdV, extractColumns, pattern="Intensity", ann=annVdV[,1])
intensity <- do.call("cbind", intensity)
###Set the column names
colnames(intensity) <- allSamplesVdV
###Check order
all(phenoVdV$SampleName == colnames(intensity))


###Create an instance of class ExpressionSet for the van't Veer cohort:

###Merge and check order
annVdV <- merge(annVdV, newAnn, by=1, all=TRUE, sort=TRUE)
rownames(annVdV) <- annVdV[,1]
all(rownames(annVdV) == rownames(logRat))
all(rownames(annVdV) == rownames(logRatError))
all(rownames(annVdV) == rownames(pVal))
all(rownames(annVdV) == rownames(intensity))
###Create the new assayData
myAssayData <- assayDataNew(exprs=logRat, exprsError=logRatError,
                            pValue=pVal, intensity=intensity)
###Create the new phenoData
myPhenoData <- new("AnnotatedDataFrame", phenoVdV)
###Create the new featureData
myFeatureData <- new("AnnotatedDataFrame", annVdV)
###Create the new experimentData
myExperimentData <- new("MIAME", name = "Marc J Van De Vijver, Yudong D He, and Laura J van't Veer", 
                        lab = "The Netherland Cancer Institute, Amsterdam, The Netherlands",
                        contact = "Luigi Marchionni <marchion@gmail.com>",
                        title = "A gene-expression signature as a predictor of survival in breast cancer", 
                        abstract = "Background: A more accurate means of prognostication in breast cancer 
                        will improve the selection of patients for adjuvant systemic therapy. 
                        Methods: Using microarray analysis to evaluate our previously established 70-gene 
                        prognosis profile, we classified a series of 295 consecutive patients with primary 
                        breast carcinomas as having a gene-expression signature associated with either a 
                        poor prognosis or a good prognosis. All patients had stage I or II breast cancer 
                        and were younger than 53 years old; 151 had lymph-node-negative disease, and 144 
                        had lymph-node-positive disease. We evaluated the predictive power of the prognosis 
                        profile using univariable and multivariable statistical analyses. 
                        Results: Among the 295 patients, 180 had a poor-prognosis signature and 115 had a 
                        good-prognosis signature, and the mean (±SE) overall 10-year survival rates were 
                        54.6±4.4 percent and 94.5±2.6 percent, respectively. At 10 years, the probability 
                        of remaining free of distant metastases was 50.6±4.5 percent in the group with a 
                        poor-prognosis signature and 85.2±4.3 percent in the group with a good-prognosis 
                        signature. The estimated hazard ratio for distant metastases in the group with a 
                        poor-prognosis signature, as compared with the group with the good-prognosis 
                        signature, was 5.1 (95 percent confidence interval, 2.9 to 9.0; P<0.001). This 
                        ratio remained significant when the groups were analyzed according to lymph-node 
                        status. Multivariable Cox regression analysis showed that the prognosis profile 
                        was a strong independent factor in predicting disease outcome. 
                        Conclusions: The gene-expression profile we studied is a more powerful predictor 
                        of the outcome of disease in young patients with breast cancer than standard  
                        systems based on clinical and histologic criteria.",
                        url = "http://www.ncbi.nlm.nih.gov/pubmed/?term=12490681",
                        pubMedIds = "12490681")

###Create the expression set
vanDeVijver <- new("ExpressionSet",
                   assayData = myAssayData,
                   phenoData = myPhenoData,
                   featureData = myFeatureData,
                   experimentData = myExperimentData)


#########################################
### PHENOTYPIC INFORMATION PROCESSING ###
#########################################

########################################################
## Patients' phenotypes in the vantVeer ExpressionSet ##
########################################################

###Define the data set type from file of origin
type <- gsub("..txt", "", gsub(".+ArrayData_", "", filesVtV))
dataSetType <- mapply(x=samplesVtV, y=type, FUN=function(x,y) {
  rep(y, length(x))
  })
###Combine with sample information
dataSetType <- do.call("c", dataSetType)
names(dataSetType) <- allSamplesVtV
###Reorder
dataSetType <- dataSetType[order(names(dataSetType))]
###Add the information to pData(vantVeer)
if (all(rownames(pData(vantVeer)) == names(dataSetType) )) {
  pData(vantVeer)$DataSetType <- dataSetType
  print("Adding information about data set type to pData")
} else {
  print("Check order pData and data set type information")
}

###Process time metastases (TTM)
pData(vantVeer)$TTM <- pData(vantVeer)$followup.time.yr
####Process TTM event
pData(vantVeer)$TTMevent <- pData(vantVeer)$metastases
####Create binary TTM at 5 years groups
pData(vantVeer)$FiveYearMetastasis <- pData(vantVeer)$TTM < 5 & pData(vantVeer)$TTMevent == 1
###Show structure of updated phenotypes
str(pData(vantVeer))
###Save the final ExpressionSet object
dataDirLoc <- system.file("data", package = "seventyGeneData")
save(vantVeer, file=paste(dataDirLoc, "/vantVeer.rda", sep=""))


############################################################
## Patients' phenotypes in the vanDeVijver ExpressionSet ##
############################################################

###Select new cases not included in the van't Veer study
pVDV <- pData(vanDeVijver)
###Rename columns
selNames <- c("TIMEmeta", "EVENTmeta", "TIMEsurvival", "EVENTdeath", "TIMErecurrence")
newNames <- c("TTM", "TTMevent", "OS", "OSevent", "RFS")
colnames(pVDV)[ sapply(selNames, grep, colnames(pVDV)) ] <- newNames
###Process time metastases (TTM)
pVDV$TTM[is.nan(pVDV$TTM)] <- pVDV$OS[is.nan(pVDV$TTM)]
###Process recurrence free survival (RFS) adding RFSevent
pVDV$RFSevent <- pVDV$RFS < pVDV$OS
###Create binary TTM at 5 years groups selecting:
###1) the cases with metastases as first event within 5 years
badCases <- which(
  pVDV$TTM <= pVDV$RFS ###Met is 1st recurrence
  & pVDV$TTMevent == 1 ### Metastases occurred
  & pVDV$TTM < 5 ### Recurrence within 5 years
)
###2) the cases disease free for at least 5 years
goodCases <- which(
  pVDV$TTM > 5 ### No metastasis before 5 years
  & pVDV$RFS > 5 ###No recurrence before 5 years
  & pVDV$TTMevent == 0 ### Metastases did notoccurred
)
###Check if there are duplicated cased present in both prognostic groups
all (!goodCases %in% badCases)

###Create groups by setting all cases to NA and then identifying bad cases
pVDV$FiveYearMetastasis <- NA
pVDV$FiveYearMetastasis[badCases] <- TRUE
###And then excluding patients with a relapse before a metastasis within 5 years
pVDV$FiveYearMetastasis[goodCases] <- FALSE
###Assign updated phenotypic data
pData(vanDeVijver) <- pVDV
###Show structure of updated phenotypes
str(pData(vanDeVijver))
###Save the final ExpressionSet object
dataDirLoc <- system.file("data", package = "seventyGeneData")
save(vanDeVijver, file=paste(dataDirLoc, "/vanDeVijver.rda", sep=""))



