library(seventyGeneData)
data(vantVeer)

show(vantVeer)

###Extract the 78 selected tumors from vantVeer data
met5y <- pData(vantVeer)$DataSetType
less5year <- which(met5y=="less_than_5y")
greater5year <- which(met5y=="greater_than_5y")
tumors78 <- sort(c(less5year,greater5year))

###ExpressionSet whith only 78 selected tumors
vantVeer78tumors <- vantVeer[,tumors78]
summary(pData(vantVeer78tumors))

###Save pData in a variable
pData78 <- pData(vantVeer78tumors)

###Check dataSetType for the 78 tumors & change in to factor type
dataType <- pData78$DataSetType
dataType <- as.factor(dataType)
summary(dataType)
#reordenamos el vector dataType con respecto al número de muestra
dataType <- dataType[order(pData78$Sample)]
dataType

####CLASE A ESTUDIAR: variable dataType! 

###Reorder expression data, samples in rows (78), genes are column names
expData78 <- t(exprs(vantVeer78tumors))
samples <- rownames(expData78)
samples.num <- as.numeric(sub(".* ", "", samples))
rownames(expData78) <- samples.num
expData78 <- expData78[order(as.numeric(rownames(expData78))), ]
row.names(expData78) <- seq(nrow(expData78))

###Transform expression data into numeric with NA's
expData78 <- as.data.frame(expData78)
expData78[expData78=="   NaN"] <- NA
expData78 <- as.data.frame(lapply(expData78, function(x) sub(",", ".", x)))
expData78 <- as.data.frame(lapply(expData78, as.character), stringsAsFactors=FALSE)
expData78 <- as.data.frame(lapply(expData78, as.numeric))


###Función para contar los NA en cada fila
find.na <- function(data){
  nas <- vector()
  for(i in 1:78){ 
    #print(length(which(is.na(data[i,]))))
    nas <- c(nas, length(which(is.na(data[i,]))))
  }
  return(nas)
}

row.nas <- find.na(expData78)
row.nas

#nos damos cuenta de que en la mayoría de las muestras el número de na's es 293 o muy cercano
#vamos a comprobar si se trata de los mismos genes en todas ellas:

which((is.na(expData78[1,]))) %in% which((is.na(expData78[5,])))

#parece que esos 293 genes se repiten como valores perdidos en todas (o casi todas) las muestras
#por lo tanto los eliminaremos:

omitedGenes <- which((is.na(expData78[1,])))
expData78clean <- expData78[ ,-omitedGenes]

##mostramos los NA por muestra del nuevo conjunto:
row.nas.cleaned <- find.na(expData78clean)
row.nas.cleaned

##el resto de NAs los eliminamos por imputación, utilizaremos el valor medio de expresión de cada gen
for(i in 1:24188){
  genAux <- expData78clean[,i]
  meanExp <- summary(genAux)[4]
  expData78clean[is.na(expData78clean[,i]),i] <- meanExp
}

##volvemos a mostrar los NAs por muestra del conjunto, tras la imputación:
row.nas.imp <- find.na(expData78clean)
row.nas.imp
##podemos observar que hemos eliminado todos los datos perdidos.


###Data Normalization 

expData78.norm <- sapply(expData78clean[,-24189], Normalize <- function(dataSet){
  datanorm <-( (dataSet-min(dataSet)) / (max(dataSet)-min(dataSet)) )
})

expData78.norm <- as.data.frame(expData78.norm)


###Data Filtration

#Conversión de la clase a binario ("greater_than_5y" = 0, "less_than_5y" = 1)
dataType.binary <- (dataType =="less_than_5y")*1

#cálculo de la correlación
correlation <- sapply(expData78.filt, function(i){
  cor(i,dataType.binary)
})

### Subconjunto de 1000 genes:

#Número de genes que queremos seleccionar
gene.number <- 1000
#Ordenación por correlación (en valor absoluto) y selección de los 200 primeros
best.correlation <- sort(abs(correlation), decreasing = TRUE)[1:gene.number] 
#Nombres de los genes seleccionados
best.correlation.names <- names(best.correlation)
#Extracción de la matriz de datos de los genes seleccionados
expData78.best1000.corr <- expData78.norm[,best.correlation.names]


### Subconjunto de 70 genes:

#Número de genes que queremos seleccionar
gene.number <- 90
#Ordenación por correlación (en valor absoluto) y selección de los 200 primeros
best.correlation <- sort(abs(correlation), decreasing = TRUE)[1:gene.number] 
#Nombres de los genes seleccionados
best.correlation.names <- names(best.correlation)
#Extracción de la matriz de datos de los genes seleccionados
expData78.best90.corr <- expData78.filt[,best.correlation.names]


## Visualización la expresión de los genes seleccionados en un heatmap
library(gplots)
color.map <- function(class) { if (class=="greater_than_5y") "#0000FF" else "#FF8000" }
classcolors <- unlist(lapply(dataType, color.map))
heatmap(t(as.matrix(expData78.best90.corr)), Colv=NA, col=redgreen(75), ColSideColors=classcolors)



library(randomForest)
# Añadimos la clase al data frame
expData78.best1000.corr$Class <- dataType.binary
# prepare training scheme
control.train <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the model
model <- train(Class~., data=expData78.best1000.corr, method="knn", trControl=control.train)
# estimate variable importance
importance <- varImp(model, scale=FALSE)
# plot importance
plot(importance, top=70)


##########################
### NEW GENE SELECTION ###
##########################

library(genefilter)


###Specific filter selection
###Genes differentially expressed in the groups defined by DataSetType, with p-value=0.01

filter <- ttest(dataType.binary, p=0.1)
selection <- genefilter(t(expData78.norm), filterfun(filter))
sum(selection)

expData78.filt <- expData78.norm[,selection]


###Non-specific filter selection
###Genes with expression measure above 0.6 in at least 3 tumors
nsFilter <- kOverA(5, 0.8)
nsSelection <- genefilter(t(expData78.norm), filterfun(nsFilter))
sum(nsSelection)

###Combine the two filters
combinedFilter <- filterfun(filter, nsFilter)
combinedSelection <- genefilter(t(expData78.norm), combinedFilter)
sum(combinedSelection)
which(combinedSelection)

###Check for coincidences
###Extract the 70 selected genes from vantVeer data
gene70 <- which(featureData(vantVeer)@data$genes70)
expData.vantveer70 <- expData78[,gene70]
###
selection1 <- unname(which(nsSelection))
gene70 %in% selection1
length(which(gene70 %in% selection1))


knnCV <- function(EXPR, selectfun, cov, Agg, pselect = 0.01, Scale=FALSE) {
  nc <- ncol(EXPR)
  outvals <- rep(NA, nc)
  for(i in 1:nc) {
    v1 <- EXPR[,i]
    expr <- EXPR[,-i]
    glist <- selectfun(expr, cov[-i], p=pselect)
    expr <- expr[glist,]
    if( Scale ) {
      expr <- scale(expr)
      v1 <- as.vector(scale(v1[glist]))
    }
    else
      v1 <- v1[glist]
    out <- paste("iter ",i, " num genes= ", sum(glist), sep="")
    print(out)
    Aggregate(row.names(expr), Agg)
    if( length(v1) == 1)
      outvals[i] <- knn(expr, v1, cov[-i], k=5)
    else
      outvals[i] <- knn(t(expr), v1, cov[-i], k=5)
  }
  return(outvals)
}

gfun <- function(expr, cov, p=0.01) {
  f2 <- ttest(cov, p=p)
  ffun <- filterfun(f2)
  which <- genefilter(expr, ffun)
}

library(class)
geneData <- genescale(exprs(vanDeVijver), 1)
Agg <- new ("aggregator")
testcase <- knnCV(geneData, gfun, pData(vanDeVijver)$FiveYearMetastasis, Agg, pselect=0.01)

library(CMA)


expData78.best90.corr$Class <- dataType

library(pROC)
library(nnet)
nn.fit <- nnet(Class~., data=expData78.best90.corr, size=2, maxit=20)
nn.pred <- predict(nn.fit, expData78.best90.corr, type="raw")
nn.roc <- roc(expData78.best90.corr$Class, nn.pred, smooth=FALSE, auc=TRUE)
plot(nn.roc, main="Redes Neuronales - Curva ROC")
auc(nn.roc)


gene70 <- which(featureData(vantVeer)@data$genes70)
expData.vantveer70 <- expData78[,gene70]

for(i in 1:70){
  genAux <- expData.vantveer70[,i]
  meanExp <- summary(genAux)[4]
  expData.vantveer70[is.na(expData.vantveer70[,i]),i] <- meanExp
}

expData.vv70.norm <- sapply(expData.vantveer70, Normalize <- function(dataSet){
  datanorm <-( (dataSet-min(dataSet)) / (max(dataSet)-min(dataSet)) )
})

expData.vv70.norm <- as.data.frame(expData.vv70.norm)
expData.vv70.norm$Class <- dataType

nn70.fit <- nnet(Class~., data=expData.vv70.norm, size=2, maxit=20)
nn70.pred <- predict(nn70.fit, expData.vv70.norm, type="raw")
nn70.roc <- roc(expData78.best90.corr$Class, nn70.pred, smooth=FALSE, auc=TRUE)
plot(nn.roc, main="Redes Neuronales - Curva ROC")
auc(nn.roc)
