################################################################
##        Main script of the R PACKAGE : RTools4TB
##
## Date : 8 / 10 / 2008
## Authors : BERGON A. 
##  with the collaboration of LOPEZ F., TEXTORIS J., GRANJEAUD S. and PUTHIER D.
## Laboratory : TAGC / INSERM U928
##              163, Avenue de Luminy, case 928
##              13 288 Marseille cedex 09 - FRANCE
##
## R package require : Biobase, tools
##
## R CMD SHLIB dbf.c -o dbf
##
## Package version : 1.1.3
## Last modification : 30 / 06 / 2009
## News :
##   * 1.1.0 to 1.1.1
##      o Add spm and spgm distance methods
##      o Add 'fac' argument into heatmapFromCDT and plotGeneExpProfiles
##   * 1.1.1 to 1.1.2
##      o Extract (limma's) normalizeQuantiles function from the library code
##   * 1.1.2 to 1.1.3
##	    o Modification of the name of 2 methods:
##		    cluster ==> clusterEisen
##		    getInfo ==> getTBInfo
##   * 1.1.3 to 1.1.4
##      o Reduce time to process vignette
##      o Correction of a bug in 'DBFMCLresult-class.Rd'
#################################################################
## verification que la longueur de la ligne est inferieur a 80 caracteres
## perl -ne 'chomp; print length $_,"\t",$_,"\n" if(length $_ > 80) ;'  new.R

#################################################################
##    DEFINITION OF A SPECIFIC CLASS OBJECT : DBFMCLresult
#################################################################

setClass("DBFMCLresult",

   representation=list(
       name="character",
       data="matrix",
       cluster="vector",
       size="vector",
       center="matrix",
       parameters="list"
## qui contiendra:
#                normalizationMethod="character",
#                distanceMethod="character",
#                inflation="numeric",
#                FDR="numeric",
#                knn="numeric",
#                random="numeric",
#                setSeed="numeric"
   ),
   prototype=list(
       name=character(),
       data=matrix(nr=0, nc=0),
       cluster=numeric(),
       size=numeric(),
       center=matrix(nc=0, nr=0),
       parameters=list()
## qui contiendra:
#              normalizationMethod=character(),
#              distanceMethod=character(),
#              inflation=numeric(),
#              FDR=numeric(),
#              knn=numeric(),
#              random=numeric(),
#              setSeed=numeric()
   )
) 
#################################################################
##    REDEFINE SHOW METHOD FOR CLASS OBJECT : DBFMCLresult
#################################################################

setMethod("show", signature("DBFMCLresult"),

    function(object){

        cat("\t\tAn object of class DBFMCLresult\n")
        cat("\t\tName:", slot(object, "name"), "\n")
        cat("\t\tMemory used: ", object.size(object), "\n")
        cat("\t\tNumber of samples: ", ncol(slot(object, "data")), "\n")
        cat("\t\tNumber of informative genes: ", 
            nrow(slot(object, "data")),"\n")
        cat("\t\tNumber of clusters: ", length(slot(object, "size")), "\n")
        cat("\t\tThis object contains the following informations:\n")

        for (i in slotNames(object)){
            cat("\t\t\t - ", i, "\n")
        }
        if(length(slot(object, "parameters"))>0){
            for(i in 1:length(slot(object, "parameters"))){
                cat("\t\t\t\t * ", names(slot(object, "parameters"))[[i]], 
                    " = ", slot(object, "parameters")[[i]], "\n")
            }
        }
    }
)

#################################################################
##    PRINCIPAL METHOD WICH ALLOWS TO RUN DBF-MCL
##   take an ExpressionSet, a data.frame, a file or a matrix 
##   (convertion into matrix format)
#################################################################

createSignatures4TB <- function(data=NULL, filename=NULL, path=".", 
            name=NULL, 
            normalization=c("rank", "gaussian", "quantiles", "none"), 
            distance.method=c("spearman", "pearson", "euclidean", "spm", "spgm"), 
            silent=FALSE, verbose=TRUE, k=150, random=3, memory.used=1024, 
            fdr=10, inflation=2.0, median.center=FALSE, set.seed=123, 
            returnRank=FALSE){

    ## testing the system
    if(.Platform$OS.type != "windows"){
        ## changing working directory
        dirAtStart <- getwd()
        setwd(path)

        ## getting parameters
        normalization <- match.arg(normalization)
        distance.method <- match.arg(distance.method)
        if(is.null(inflation)) inflation <- 2.0

        ## getting data
        if(verbose) cat("Reading data...\n\n")
        sourceData <- getData4DBFMCL(data=data, filename=filename, path=".")
        ## data (matrix format)
        m <- sourceData$data
        ## experiment name 
        ## By default "exprs"
        if(is.null(name)) name <- sourceData$name
        if(is.null(name)) name <- "exprs"

        ## writting all parameters 
        if(verbose){
            cat("The following parameters will be used :",
            "\n\tWorking directory: ", getwd(),
            "\n\tName: ", name,
            "\n\tNormalization method: ", normalization,
            "\n\tDistance method: ", distance.method,
            "\n\tNumber of neighbors: ", k,
            "\n\tNumber of randomizations: ", random,
            "\n\tFDR: ", fdr, "%",
            "\n\tInflation: ", inflation,
            "\n\tMedian center rows: ", median.center,
            "\n\tVisualize standard outputs from both mcl and cluster", 
                        "commands: ", silent,
            "\n\tMemory used : ", memory.used, "\n\n") 
        }
        ## Normalization (by rank, gaussian, quantiles)
        if(normalization == "rank"){
            if(verbose) cat("Normalizing on ranks... ")
            m <- doRankTransformation(m)
            if(verbose) cat("Done.\n\n")
        }

        else if(normalization == "quantiles"){
            ## normalization by Quantiles method
            if(verbose) cat("Performing quantiles normalization... ")
            m <- normalizeQuantiles(m, ties=TRUE)
            if(verbose) cat("Done.\n\n")
        }

        else if(normalization == "gaussian"){
            ## Gaussian transformation 
            if(verbose) cat("Performing gaussian normalization... ")
            m <- doNormalScore(m, set.seed=set.seed)
            if(verbose) cat("Done.\n\n")
        }

        ## Creating a DBFMCLresult object
        obj <- DBFMCL(data=m, name=name,
                    distance.method=distance.method, 
                    clustering=TRUE, silent=silent, verbose=FALSE, k=k,
                    random=random, memory.used=memory.used, fdr=fdr, 
                    inflation=inflation, set.seed=set.seed, 
                    returnRank=returnRank)

        obj@parameters <- list(
            normalizationMethod=normalization,
            distanceMethod=distance.method,
            k=k, 
            random=random, 
            fdr=fdr,
            set.seed=set.seed, 
            inflation=inflation)

        if(length(obj@size) > 0){
            ## SortBySignatures
            if(verbose) cat("Preparing file for treeview... \n")
            writeDBFMCLresult(obj, paste(name, ".dataMods.txt", sep=""),
                path=".")
            if(verbose) cat("Done file...\n\n")
    
            ## cluster
            if(verbose) cat("clustering columns (samples)... ")
            clusterEisen(paste(name, ".dataMods.txt", sep=""),
                median.center=median.center)
            if(verbose) cat("Clustering done...\n\n")
        }
        ## coming back to the working directory used at start
        setwd(dirAtStart)
    }
    else{
        stop("mcl and cluster programs requires a unix-like OS.")
    }
    return(obj)
}

#################################################################
##    DBF-MCL
#################################################################

DBFMCL <- function(data=NULL, filename=NULL, path=".", name=NULL,
           distance.method=c("pearson", "spearman", "euclidean", "spm", "spgm"), 
           clustering=TRUE, silent=FALSE, verbose=TRUE, k=150,
           random=3, memory.used=1024, fdr=10, inflation=2.0,
           set.seed=123, returnRank=FALSE){

    ## testing the system
    if(.Platform$OS.type != "windows"){

        ## getting parameters 
        dataSource <- getData4DBFMCL(data=data, filename=filename, path=path)
        m <- dataSource$data
        if(is.null(name)) name <- dataSource$name
        if(is.null(name)) name <- "exprs"

        distance.method <- match.arg(distance.method)
        if(clustering){
            if(is.null(inflation)) inflation <- 2.0
            txt <- paste("\n\tInflation: ", inflation, sep="")
        }
        else{
            txt <- ""
        }
        ## writting all parameters 
        if(verbose){
#             cat("\nThe working directory is: ", getwd(), "\n\n")
#             cat("The following parameters will be used:",
#                 "\n\tName: ", name,
#                 "\n\tCorrelation method: ", distance.method,
#                 "\n\tNumber of neighbors: ", k,
#                 "\n\tNumber of randomizations: ", random,
#                 "\n\tFDR: ", fdr, "%",
#                 "\n\tVisualize standard outputs from mcl ",
#                         "and cluster commands: ", silent,
#                 "\n\tMemory used: ", memory.used,
#                 "\n\n")
            
            cat("The following parameters will be used :",
                "\n\tWorking directory: ", getwd(),
                "\n\tName: ", name,
                "\n\tDistance method: ", distance.method,
                "\n\tNumber of neighbors: ", k,
                "\n\tNumber of randomizations: ", random,
                "\n\tFDR: ", fdr, "%",
                "\n\tPerform clustering: ", clustering, txt,
                "\n\tVisualize standard outputs from both mcl and cluster", 
                            "commands: ", silent,
                "\n\tMemory used : ", memory.used, "\n\n") 
        }

        ##DBF algorithm, returns a DBFMCLresult object
        obj <- DBF(m, name,
                    distance.method=distance.method, 
                    silent=silent, 
                    k=k, 
                    random=random, 
                    memory.used=memory.used,
                    fdr=fdr,
                    set.seed=set.seed, returnRank=returnRank)
        if(length(readLines(paste(name,".dbf_out.txt", sep=""))) > 0){
            ## RUN MCL ???
            if(clustering){
                ## Launching mcl
                if(is.null(inflation)) inflation <- 2.0
                MCL(name, inflation=inflation, silent=silent)
                
    
                ## getting mcl results into the DBFMCLresult object
                mclCluster <- readLines(paste(name, ".mcl_out.txt", sep=""))
                sondeList <- NULL
                clusters <- NULL
                size <- NULL
                nb <- 0
    
                for( i in 1:length(mclCluster)){
                    h <- unlist(strsplit(mclCluster[i], "\t"))
                    if(length(h) >= 10){ 
                        nb <- nb + 1
                        sondeList <- c(sondeList, h)
                        clusters <- c(clusters, rep(nb, length(h)))
                        if(is.null(size)){
                            size <- length(h)
                        }
                        else{
                            size <- c(size, length(h))
                        }
                    }
                }
                if(verbose) cat(nb, " signatures containing at least ",
                    "10 probes will be conserved\n\n")
    
                ## build DBFMCLresult object
                if(nb > 0){
		            if(returnRank){
		                m <- round(doRankTransformation(m))
		            }
                    obj@name <- name
                    obj@data <- as.matrix(m[sondeList, ])
                    obj@cluster <- clusters
                    obj@size <- size
    
                    centers <- matrix(nc=ncol(m), nr=nb)
                    ## calcul of the mean profils
                    for(i in 1:nb){
                        centers[i, ] <- apply(obj@data[obj@cluster == i, ], 
                                            2, mean, na.rm=TRUE)
                    }
                    obj@center <- centers
    
                    ## add DBFMCL parameters used to build this object
                    obj@parameters <- list(
                        distanceMethod=distance.method,
                        k=k, 
                        random=random, 
                        fdr=fdr,
                        set.seed=set.seed, 
                        inflation=inflation)
                }
            }
            else{
                ## add only DBF parameters used to build this object
                obj@parameters <- list(
                        distanceMethod=distance.method,
                        k=k, 
                        random=random, 
                        fdr=fdr,
                        set.seed=set.seed)
            }
        }
        else{
            stop("There is no conserved gene.\n\n")
        }
        return(obj)
    }
    else{
        stop("A unix-like OS is required to launch mcl and cluster programs.")
    }
}


#################################################################
##    Getting data
#################################################################

getData4DBFMCL <- function(data=NULL, filename=NULL, path="."){

    ## getting matrix (probesID vs SamplesID) 
    if(!is.null(data)){
        if(inherits(data, "ExpressionSet")){
            data <- as.matrix(exprs(data))
        }
        else if(is.data.frame(data)){
            data <- as.matrix(data)
        }
        if(!is.matrix(data)){
            stop("Please provide an ExpressionSet, a data.frame", 
                " or a matrix.\n")
        }
        name <- NULL
    }
    else{
        if(!is.null(filename)){
            data <- as.matrix(read.table(file.path(path, filename), 
                        sep="\t", head=TRUE, row.names=1, quote=""))
            name <- unlist(strsplit(filename, "\\."))[1]
        }
        else{
            stop("Please provide an ExpressionSet, a data.frame, ", 
                "a matrix or a tabular file\n")
        }
    }
    ## adding dimnames if not provided
    if(is.null(rownames(data)))
       rownames(data) <- paste("gene", 1:nrow(data), sep="")
    if(is.null(colnames(data)))
       colnames(data) <- paste("sample", 1:ncol(data), sep="")

    return(list(data=data, name=name))
}

#################################################################
##    Transformation by rank of the data
#################################################################

doRankTransformation <- function(data=NULL){
    if(is.matrix(data)){
        row <- rownames(data)
        data.rank <- data
        data.rank <- apply(data.rank, 2, rank)
        data <- apply(data.rank, 2, as.double)
        rownames(data) <- row
    }
    else{
        stop("Please provide a matrix.\n")
    }
    return(data)
}


###############################################################
##    COMPUTE DBF algorithm 
###############################################################

DBF <- function (data, name=NULL, 
            distance.method=c("spearman", "pearson", "euclidean","spm", "spgm"), 
            silent=FALSE, 
            k=100, 
            random=3, 
            fdr=10, 
            memory.used=1024,
            set.seed=123, 
            returnRank=FALSE){

    ## testing the system
    if(.Platform$OS.type != "windows"){

        if(!is.null(data)){
            ## getting data and parameters 
            if(is.null(name)) name <- "exprs"
            data <- getData4DBFMCL(data=data)$data
            row <- rownames(data)
            col <- colnames(data)
            distance.method <- match.arg(distance.method)

            ## transforming data into double
            data <- apply(data, 2, as.double)

            if(silent) cat("Computing distances to the kth-nearest neighbors",
                " and associated FDR values... \n")
            outfile <- paste(name, ".dbf_out.txt", sep="")

            

            ## launching DBF
            a <- .C("DBF", 
                    data, 
                    as.integer(nrow(data)),
                    as.integer(ncol(data)),
                    row,
                    col, 
                    distance.method, 
                    as.integer(k),
                    as.integer(random), 
                    as.integer(!silent),
                    as.integer(memory.used), 
                    as.integer(fdr),
                    as.integer(!silent), 
                    m2=vector(length=nrow(data), mode="character"), 
                    outfile,
                    as.integer(set.seed))

            ## creation of the DBFMCLresult object
            obj <- new("DBFMCLresult")

            informative <- a$m2[a$m2 != ""]
            if(length(informative) > 0){

                if(returnRank){
                    data <- round(doRankTransformation(data))
                }
                obj@data <- as.matrix(data[a[[4]]%in%informative, ])
                obj@cluster <- rep(1, nrow(obj@data))
                obj@size <- nrow(obj@data)
                obj@center <- matrix(
                                    apply(obj@data[obj@cluster == 1, ], 
                                        2,
                                        mean, 
                                        na.rm=TRUE),
                                    nr=1)
            }
            return(obj)

        }else{
            stop("Please provide a matrix...\n\n")
        }
    }
    else{
        stop("A unix-like OS is required to launch mcl and cluster programs.")
    }
}



##############################################################
##    MCL
##############################################################

MCL <- function(name, inflation=2.0, silent=FALSE){
    ## testing the system
    if(.Platform$OS.type != "windows"){

        ## Testing mcl installation
        if(system("mcl --version | grep 'Stijn van Dongen'", intern=TRUE) > 0){
            if(!silent){
                cat("Running mcl (graph partitioning)... ")
                verb <- ""
            }
            else{
                verb <- "-V all"
            }
            if(inflation!=2){
                i <-  paste("-I ", as.character(round(inflation, 1)), sep="")
            }else{
                i <- "-I 2.0"
            }
            ## launching mcl program
            system(paste("mcl ", name, ".dbf_out.txt ", i, " --abc -o ",
                    name, ".mcl_out.txt ", verb, sep=""))

            if(!silent) {
                cat("Done.\n\n")
                cat("\t--> creating file : ",
                file.path(getwd(), paste(name, ".mcl_out.txt", sep="")),
                "\n\n")
            }
        }else{
            stop("Please install mcl on your computer...\n",
                "You can download it from : 'http://www.micans.org/mcl/'\n\n")
        }
    }
    else{
        stop("A unix-like OS is required to launch mcl and cluster programs.")
    }
}



################################################################
## This is the normalizeQuantiles function from the marray Bioconductor library.
################################################################
# 
# 
# normalizeQuantiles <- function ( A, ties=TRUE){
# 
#     n <- dim(A)
#     if (is.null(n))
#         return(A)
#     if (n[2] == 1)
#         return(A)
#     O <- S <- array(, n)
#     nobs <- rep(n[1], n[2])
#     i <- (0:(n[1] - 1))/(n[1] - 1)
#     for (j in 1:n[2]) {
#         Si <- sort(A[, j], method="quick", index.return=TRUE)
#         nobsj <- length(Si$x)
#         if (nobsj < n[1]) {
#             nobs[j] <- nobsj
#             isna <- is.na(A[, j])
#             S[, j] <- approx((0:(nobsj - 1))/(nobsj - 1), Si$x,
#                 i, ties="ordered")$y
#             O[!isna, j] <- ((1:n[1])[!isna])[Si$ix]
#         }
#         else {
#             S[, j] <- Si$x
#             O[, j] <- Si$ix
#         }
#     }
#     m <- rowMeans(S)
#     for (j in 1:n[2]) {
#         if (ties)
#             r <- rank(A[, j])
#         if (nobs[j] < n[1]) {
#             isna <- is.na(A[, j])
#             if (ties)
#                 A[!isna, j] <- approx(i, m, (r[!isna] - 1)/(nobs[j] -
#                   1), ties="ordered")$y
#             else A[O[!isna, j], j] <- approx(i, m, (0:(nobs[j] -
#                 1))/(nobs[j] - 1), ties="ordered")$y
#         }
#         else {
#             if (ties)
#                 A[, j] <- approx(i, m, (r - 1)/(n[1] - 1), ties="ordered")$y
#             else A[O[, j], j] <- m
#         }
#     }
#     return(A)
# }

#######################################################################
##    Normal score transformation.
#######################################################################

doNormalScore <- function (sdata, set.seed=123){
    row <- row.names(sdata)
    set.seed(set.seed)
    reference <- sort(rnorm(nrow(sdata)))
    for(i in 1:ncol(sdata)){
        r <- sdata[,i]
        rsort <- sort(r)
        ref <- reference
        if(length(unique(r[!isUnique(r)])) >= 1){
	          seq <- unique(r[!isUnique(r)])
	          for(s in seq){
		            posS <- which(rsort%in%s)
	              ref[posS] <- mean(ref[posS])
            }
        }
        sdata[, i] <- ref[match(r, rsort)]
    }
    sdata <- as.data.frame(sdata)
    row.names(sdata) <-  row
    return(sdata)
}

# #######################################################################
# #     MCL2biolayout
# #######################################################################
# 
# 
# MCL2biolayout <- function(name){
# 
#     ## lecture des deux fichiers a analyser
#     mclCluster <- as.vector(readLines(paste(name,".MCLclusters.txt",sep="")))
#     heartsData <- as.vector(readLines(paste(name,".hearts.txt", sep="")))
#     dataout <- heartsData
#     ## formatage des donnees que l'on souhaite avoir en sortie
#     ## cad mettre a la fin du fichier contenant les cluster de MCL des lignes telles que :
#     ## //NODECLASS    sonde_id    num_cluster
# 
#     for( i in 1:length(mclCluster)){
#         h <- unlist(strsplit(mclCluster[i], "\t"))
#         if(length(h)>=10){
#             for( j in 1: length(h)){
#                 dataout <- c(dataout, paste("//NODECLASS\t",h[j],"\t",i,"\n", sep=""))
#             }
#         }
#     }
# 
#     ## inserer ici les cmd perl pour le traitement des nouveaux resultats
# 
# 
#     ## fichier de sortie des donnees processees
#     writeLines(dataout, paste(name,".heartsClustered.txt",sep=""))
#     cat("\t--> creating file : ",paste(getwd(),"/",name, ".heartsClustered.txt", sep=""),"\n\n")
# 
# }
#######################################################################
##        sortBySignatures
#######################################################################

writeDBFMCLresult <- function(object, filename.out=NULL, path=".",
                verbose=TRUE){

    if(path == ".") path <- getwd()
    if(inherits(object, "DBFMCLresult")){
        if(is.null(filename.out)) filename.out <- "exprs.dataMods.txt"

        data <- object@data
        nb <- 0 
        dataT <- c("clusters", colnames(data))

        ## processing data
        for(i in 1:length(object@size)){
            if(verbose) 
                cat("\n\tCluster ", i, " --> ", object@size[i], " probes")
            if(object@size[i] > 9){
                subData <- data[object@cluster == i, ]
                subData <- cbind(rownames(subData), subData)
                intLine <- matrix(c("NA", rep(NA,ncol(data))), nr=1)
                dataT <- rbind(dataT, subData, intLine)
                nb <- nb + 1 
            } 
        }
        if(verbose) cat("\n\n ", nb, 
            " signatures containing at least 10 probes",
            "will be kept.\n\n")

        ## exporting results
        write.table(dataT,  file.path(path, filename.out),
                    col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
        if(verbose) cat("\t--> creating file : ",
            file.path(path, filename.out), "\n\n")

    }else{
        stop("An object of class DBFMCLresult is required...")
    }
}


#########################################################
##        CLUSTER (en ligne de commande bash)
#########################################################

clusterEisen <- function (filename, median.center=FALSE, silent=FALSE){
   ## testing the system
    if(.Platform$OS.type != "windows"){
        ## Testing cluster 3.0 installation
        if(system("cluster -v |grep 'Michael Eisen'",intern=TRUE) > 0){
            ## two possibilities
            if(median.center == TRUE){
                system(paste("cluster -f ", filename,
                                " -g 0 -e 2 -m a -cg m", sep=""))
            }
            else{
                system(paste("cluster -f ", filename,
                                    " -g 0 -e 2 -m a", sep=""))
            }
            name <- unlist(strsplit(filename, "/"))
            name <- name[length(name)]
            name <- unlist(strsplit(name, ".txt"))[1]
            if(!silent) cat("\n\t--> creating files :\n\t\t",
                    file.path(getwd(), paste(name, ".cdt",sep="")),
                    "\n\t\t",
                    file.path(getwd(), paste(name, ".atr", sep="")),
                    "\n\n")
        }else{
            stop("Please install cluster (command-line version)",
                " on your computer...\nYou can download it from :",
                "'http://rana.lbl.gov/EisenSoftware.htm'\n\n")
        }
    }
    else{
        stop("A unix-like OS is required to launch mcl and cluster programs.")
    }
}


#######################################################################
##        Parser *.cdt files to allow a visualization by the heatmap method
#######################################################################


heatmapFromCDT <- function(cdt.filename, signature=NULL, fac=NULL){
    
    if(!is.null(cdt.filename)){
        ## color palette
        pal <- colorScale(low="green", high="red", mid="black", k=1000)

        ## blank line positions
        bl <- grep("NA\tNA\t1.000000\t", readLines(cdt.filename)[-c(1, 2, 3)])
        nbSignatures <- length(bl)
        cat("This file contains ", nbSignatures, " signatures...\n\n")

        ## getting data
        cdt <- as.matrix(read.table(cdt.filename, quote="", skip=3,
                                    blank.lines.skip=FALSE))
        rownames <- cdt[, 1]
        cdt <- cdt[,-c(1, 2, 3)]
        colnames(cdt) <- unlist(strsplit(
                                readLines(cdt.filename)[1], "\t"))[-c(1, 2, 3)]
        
        cdt <- apply(cdt, 2, as.numeric, re=FALSE) 

        ## add color corresponding to phenotype
        if(!is.null(fac)){
            pheno <- as.factor(fac)
            levels(pheno) <- 1:length(levels(pheno))
            cc <- rainbow(length(levels(pheno)), start=0, end=.3)
            cc <- cc[pheno]
        }else{
            cc <- rep(rainbow(1, start=0, end=.3), length(colnames(cdt)))
        }

        ## creating heatmap according to the requested signatures
        if(!is.null(signature)){
            if(length(signature) == 1){
                cat("Creating a heatmap for the signature : ", signature, 
                    " ...\n\n")
                if(signature == 1){
                    x <- cdt[1:(bl[1] - 1),]
                    r <- rownames[1:(bl[1] - 1)]
                    row.names(x) <- r
                    X11()
                    heatmap(x, Rowv = NA, 
                            main=paste(cdt.filename, "-- Signature 1"), 
                            col=pal, ColSideColors=cc)
                    output <- x
                }
                else{
                    x <- cdt[(bl[signature - 1] + 1):(bl[signature] - 1),]
                    r <- rownames[(bl[signature-1] + 1):(bl[signature] - 1)]
                    row.names(x) <- r
                    X11()
                    heatmap(x, Rowv = NA,
                            main=paste(cdt.filename,"-- Signature", signature),
                            col=pal, ColSideColors=cc)
                    output <- x
                }
            }
            else{
                for(i in signature){
                    cat("Creating a heatmap for the signature : ",
                        signature[i], " ...\n\n")
                    if(signature[i] == "1"){
                        x <- cdt[1:bl[1], ]
                        r <- rownames[1:bl[1]]
                    }
                    else{
                        x <- cdt[(bl[signature[i] - 1] +
1):(bl[signature[i]] - 1), ]
                        r <- rownames[(bl[signature[i] - 1] + 1) :
(bl[signature[i]] - 1)]    
                    }
                    X11()
                    row.names(x) <- r
                    heatmap(x, Rowv = NA,
                            main=paste(cdt.filename, "-- Signature",
                                    signature[i]),
                            col=pal, ColSideColors=cc)
                }
                output <- cdt[-bl, ]
            }
        }
        else{
            cat("Creating heatmaps for all signatures...\n\n")
            X11()
            output <- cdt[-bl, ]
            heatmap(cdt[-bl, ], Rowv = NA,
                    main=paste(cdt.filename, "-- All signatures"), 
                    col=pal, ColSideColors=cc)
        }
    }
    else{
        stop("Please provide a cdt filename.\n\n")
    }
    return(output)
}

###############################################################################
###############################################################################

colorScale <- function (low="white", high=c("green", "red"), mid=NULL, k=50)
{
    low <- col2rgb(low) / 255
    high <- col2rgb(high) / 255
    if (is.null(mid)) {
        r <- seq(low[1], high[1], len=k)
        g <- seq(low[2], high[2], len=k)
        b <- seq(low[3], high[3], len=k)
    }
    if (!is.null(mid)) {
        k2 <- round(k / 2)
        mid <- col2rgb(mid) / 255
        r <- c(seq(low[1], mid[1], len=k2), seq(mid[1], high[1], len=k2))
        g <- c(seq(low[2], mid[2], len=k2), seq(mid[2], high[2], len=k2))
        b <- c(seq(low[3], mid[3], len=k2), seq(mid[3], high[3], len=k2))
    }
    rgb(r, g, b)
}
#########################################################
#########################################################

plotGeneExpProfiles <- function(data=NULL, filename=NULL, path=".", 
                                signatures=NULL, saveHTML=FALSE,
                                filename.out=NULL, X11=TRUE, verbose=FALSE){
    
    ## getting matrix
    if(inherits(data, "DBFMCLresult")){
        nb <- length(data@size)
        m <- data@data
        name <- NULL
    }
    else{ 
        mdata <- getData4DBFMCL(data=data, filename=filename)
        m <- mdata$data
        name <- mdata$name
        nb <- ncol(m)
    }
    if(is.null(signatures)){
        signatures <- 1:nb
    }    

    ## median-centering of row
    med <- apply(m, 1, median)
    mc <- m - med
    ##  hclust
    dis <- dist(t(mc))
    h <- hclust(dis, meth="av")
    m <- mc[, h$order]
    myOut <- ""

    ##Exporting to HTML format
    if(saveHTML){
        if(is.null(filename.out)){
            if(is.null(name)){
                filename.out <- "profile_of_each_signature.html"
            }
            else{
                filename.out <- paste(name, 
                                      "profile_of_each_signature.html", sep="_")
            }
        }
        if(is.null(name)){
                imgDir <- "Profiles_images" 
        }
        else{
           imgDir <- paste("Profiles_images", name, sep="")
        }
        setwd(path)
        DirAtStart <- getwd()
        dir.create(imgDir)
        setwd(imgDir)
        myOut <- "<html><table>\n\t"
    }


    ## For each signature (obtained from a file or from a DBFMCLresult object).

    if(inherits(data, "DBFMCLresult")){
        for(i in 1:nb){
            if(length(which(i == signatures))>0){
                myPlot <- matplotProfiles(m[data@cluster == i, ], 
                        imgName=paste("profile_signature_", i, ".png", sep=""),
                        saveHTML=saveHTML, verbose=verbose, X11=X11)
                myOut <- paste(myOut, myPlot, sep="") 
            }
        }
    }
    else{ 
        ind <- 1:nrow(m)
        ind <- ind[m[,1]%in%NA == TRUE]
        if(length(ind) > 0){
            
            for(i in 1:(length(ind))){
                if(length(which(i == signatures))>0){
                    if(i == 1){
                        mSub <- m[1:(ind[i] - 1), ]
                    }
                    else{
                        mSub <- m[(ind[i-1] + 1):(ind[i] - 1), ]
                    }
                    myPlot <- matplotProfiles(mSub,
                            imgName=paste("profile_signature_", i, ".png",
sep=""),
                            saveHTML=saveHTML, verbose=verbose, X11=X11) 
                    myOut <- paste(myOut, myPlot, sep="") 
                }
            }
            ##last profile
            myPlot <- matplotProfiles( m[(ind[length(ind)] + 1):nrow(m), ],    
                   imgName=paste("profile_signature_", length(ind) + 1, ".png",
sep=""),
                    saveHTML=saveHTML, verbose=verbose, X11=X11)
            myOut <- paste(myOut, myPlot, sep="")
        }
        else{
    
            myPlot <- matplotProfiles(m,
                    imgName="profile_signature.png",
                    saveHTML=saveHTML, verbose=verbose, X11=X11)
            myOut <- paste(myOut, myPlot, sep="")
        }
    }

    
    if(saveHTML){
        myOut <- paste(myOut, "</table></html>", sep="")
        fn <- file.path(path, filename.out)
        writeLines(myOut, fn)
        browseURL(fn)
        setwd(DirAtStart)

    }

}


matplotProfiles <- function(data, imgName, saveHTML=FALSE,
        X11=TRUE, verbose=FALSE){
    ##generation of a new graphic window
    if(saveHTML == FALSE){
        if(X11){
            X11()
        }
        myOut <- 0
    }
    else{
        fnImg <- file.path(getwd(), imgName)
        png(file=fnImg, width=650, height=650)
        myOut <- "<tr align=center>\n"    
    }
    ##creating plot
    matplot(t(data), type="l", xlab="Samples", ylab="Intensity", pch=".")
    mean <- apply(data, 2, mean, na.rm=TRUE)
    lines(mean, lwd=3, col="green")
    title(imgName)    

    if(verbose){
        cat("Mean expression profile : ", imgName, "\n")
        cat("Num. \t SampleID \t Mean intensity\n")
        for(i in 1:ncol(data)){
            cat("[", i, "]\t", colnames(data)[i], "\t", mean[i], "\n")
        }
        cat("Genes contained in this signature :\n")
        cat(rownames(data), "\n", fill=TRUE, sep="\t")
    }
    ## Add this plot to the html file and save it as png
    if(saveHTML){
        dev.off()
        myOut <- paste(myOut, "\t<img src=", imgName, 
                    ">\n</tr>\t\n", sep="")
    }
    return(myOut)
}

#########################################################
##      getSignatures
#########################################################
getSignatures <- function(
            field=c("gene", "probe", "platform", "experiment", "annotation"),
            value=NULL, qValue=NULL, nbMin=NULL, verbose=TRUE,
            save=FALSE){
    ## Web service address
    url <- "http://tagc.univ-mrs.fr/tB/TBWS/servlets/TBWS?"
    if(!is.null(value)){
        ## data processing before request
        tab <- NULL
        field <- match.arg(field)
        if(field == "platform") field <- "GPL"
        if(field == "experiment") field <- "GSE"
        request <- value
        if((field == "gene")||(field == "probe")||(field == "GPL")||
            (field == "GSE")||(field == "annotation")){
            if((field == "gene")||(field == "probe")){
                if(length(unlist(strsplit(value, " "))) > 1){
                    str <- unlist(strsplit(value, " "))
                    value <- str[1]
                    for(i in 2: length(str)){
                        if(length(unlist(strsplit(value, "[&|!]"))) > 1){
                                    value <- paste(value, str[i], sep="")
                        }else{
                            value <- paste(value, str[i], sep="%20")
                        }
                    } 
                }
                if(length(unlist(strsplit(value, "&"))) > 1){
                    str <- unlist(strsplit(value, "&"))
                    value <- str[1]
                    for(i in 2: length(str)){
                        value <- paste(value, "*", str[i], sep="")
                    } 
                }
            }
            if(field == "annotation"){
                while(length(unlist(strsplit(value, "&"))) > 1){
                    value <- sub("&", "^*^", value)
                }
                value <- paste("^", value, "^", sep="")
                value <- unlist(strsplit(value, "|", fixed=TRUE))
                text <- value[1]
                if(length(value) > 1){
                    for(i in 2:length(value)) 
                        text <- paste(text, value[i], sep="^|^")
                }
                value <- text
                while(length(unlist(strsplit(value, " "))) > 1){
                    value <- sub(" ", "%20", value)
                }
                if(is.null(qValue)) qValue <- 0
                filename <- paste("TBRequest__", field, "_", request, sep="")
                url <- paste(url, "type=", field, "&request=", value, 
                "&qValue=", qValue, sep="")
            }
            else{
                if(!is.null(nbMin) && field == "gene"){
                    filename <- paste("TBRequest__", field, "_", value,
                                    "nbMin", nbMin, sep="")
                    url <- paste(url, "type=", field, "&request=", value,
                                "&nbMin=", nbMin, sep="")
                }
                else{
                    filename <- paste("TBRequest__", field, "_", value, sep="")
                    url <- paste(url, "type=", field, "&request=", value,sep="")
                }
            }
            if(save){
                download.file(url, destfile=filename, cacheOK=FALSE,
                            quiet=!verbose)
            }
            if(length(readLines(url)) > 0){
                ## perform request
                if(is.null(nbMin))tab <- readLines(url)
                else {
                    tab <- read.table(url, sep="\t", quote="")
                    colnames(tab) <- c("Signature", "nb.Genes")
                }
                if(verbose){
                    if(is.null(nbMin)) cat(length(tab))
                    else cat(nrow(tab))
                    cat(" signatures were found for the request:\n",
                            field, " = ", request)
                    if(!is.null(nbMin)) cat(" and nbMin = ", nbMin)
                    cat("\n\n")
                }
            }
            else{
                if(verbose){
                    cat("No result were found for :",
                    field, " = ", request, "\n\n")
                }
            }
        }
        else{
            stop("WARNING...Field should be one of : 'gene', 'probe', 'data', 
                'platform', 'experiment' and 'annotation' \n\n")
        }
    }
    else{
        stop("WARNING... ",
            "Please provide a query to the TBrowser database...\n\n")
    }
    return(tab)
}

#########################################################
##      getData
#########################################################

getExpressionMatrix <- function(signatureID=NULL, verbose=TRUE, save=FALSE){
    ## Web service address
    url <- "http://tagc.univ-mrs.fr/tB/TBWS/servlets/TBWS?"
    if(!is.null(signatureID)){
        ## data processing before request
        tab <- NULL
        filename <- paste("Expression_Matrix_signatureID_", signatureID, sep="")
        url <- paste(url, "type=data&request=", signatureID, sep="")
        if(save){
            download.file(url, destfile=filename, cacheOK=FALSE, quiet=!verbose)
        }
        if(readLines(url)[1]!= "ERROR : Invalid request"){
            ## request
            PG <- readLines(url)[1:2]
            probes <- unlist(strsplit(PG[1], "\t"))
            genes <- unlist(strsplit(PG[2], "\t"))
            tab <- cbind(probes, genes, 
                        read.table(url, quote="", sep="\t",
                                    head=TRUE, skip=2))
            if(verbose){
                cat("Downloading expression matrix for",
                "transcriptional signature: ", signatureID, 
                "(",
                length(colnames(tab))-2, "samples x",
                nrow(tab), "probes)\n")
            }
        }
        else{ 
            if(verbose){
                cat("The transcriptional signature ",signatureID, 
                "doesn\'t exist\n")
            }
        }
    }
    else{
        stop("WARNING... ",
            "Please provide a transcriptional signature ID",
            " from the TBrowser database...\n")
    }
    return(tab)
}

#########################################################
##      getTBInfo
#########################################################
getTBInfo <- function(field=c("platform", "experiment", "signature", "samples"),
    value=NULL, verbose=TRUE, save=FALSE){
    ## Web service address
    url <- "http://tagc.univ-mrs.fr/tB/TBWS/servlets/TBWS?"
    if(!is.null(value)){
        ## Processing data before request
        tab <- NULL
        field <- match.arg(field)
        filename <- paste("Info__", field, "_", value, sep="")
        url <- paste(url, "type=", field, "&request=", value,sep="")
        if(save){
            download.file(url, destfile=filename, cacheOK=FALSE, quiet=!verbose)
        }
        if(readLines(url)[1]!= "ERROR : Invalid request"){
            ## Request
            if(field == "samples"){
                tab <- read.table(url, quote="", sep="\t", head=FALSE)
                tab <- as.matrix(tab)
                colnames(tab) <- c("sampleID", "Title", "Description", "Source")
            }
            else{
                tab <- read.table(url, quote="", sep="\t", head=FALSE)
                view <- paste(tab[, 1], " = ", tab[, 2])
                if(verbose){
                    cat("A result was found for :",
                            field, " = ", value, "\n")
                    writeLines(strwrap(view, width=80, exdent=5))
                }
                tab <- as.matrix(tab)
                row.names(tab) <- tab[, 1]
                tab <- as.matrix(tab[, 2])
            }
        }
        else{
            if(verbose){
                cat("This", field, "is not in the TBrowser database\n\n")
            }
        }
    }
    else{
        stop("WARNING... ",
            "Please provide a query to the TBrowser database...\n\n")
    }
    return(tab)
}


#########################################################
##                    createGraph4BioC
#########################################################

createGraph4BioC <- function(request=NULL, prop=50){

    if(!is.null(request)){
        
        ## Getting signature list corresponding to the request
        reqSign <- getSignatures(field="gene", value=request)
        nbSign <- length(reqSign)
        listsGenes <- new.env()
        totaleList <- NULL
        ## Composition of these signatures
        for(i in 1:nbSign){
            ## Getting genes falling into a signature
            reqGenes <- getExpressionMatrix(signatureID = reqSign[i], 
                            verbose = FALSE, save = FALSE)
            listsGenes[[as.character(reqSign[i])]] <- unique(reqGenes[, 2])
            totaleList <- unique(c(totaleList, reqGenes))
        }
        listsGenes <- as.list(listsGenes)
        cat(length(totaleList), 
            " differents genes.\n")

        ## list of the most frequent genes
        ## to build the adjacent matrix
        genes <- unlist(listsGenes)
        tabCount <- sort(
                        table(genes)[table(genes) >= round(nbSign*prop/100)])
        cat(length(tabCount), 
            "genes fall in more than ", 
            round(nbSign*prop / 100), 
            " signatures...\n")
        conservedGenes <- rownames(tabCount)
        nbGenes <- length(conservedGenes)

        ## build the boolean matrix : 
        ## conserved genes in the signatures? (1==TRUE et 0==FALSE)
        corrSign <- matrix(0, nc=nbSign, nr=nbGenes)
        rownames(corrSign) <- conservedGenes
        colnames(corrSign) <- reqSign
        for(i in 1:nbSign){
            corrSign[conservedGenes%in%listsGenes[[i]], i] <- 1
        }

        ##build adjacence matrix
        adjMat <- matrix(NA, nc=nbGenes, nr=nbGenes)
        rownames(adjMat) <- conservedGenes
        colnames(adjMat) <- conservedGenes
        for(i in 1:nbGenes){
            for(j in 1:nbGenes){
                if(i == j){
                    adjMat[i, j] <- 0 
                }else{
                    tabC <- table(corrSign[i, ] == corrSign[j, ])
                    if(length(grep("TRUE", unlist(dimnames(tabC))))>0)
                        adjMat[i, j] <-
                            tabC[[grep("TRUE", unlist(dimnames(tabC)))]]
                    else adjMat[i, j] <- 0
                }
            }
        }
        ## delete gene which correspond to "NA"
        adjMat <- adjMat[-grep("N/A", rownames(adjMat)),
                        -grep("N/A", colnames(adjMat))]
        return(adjMat)
    }else{
        stop("Please provide a query to the TBrowser database...\n")
    }
}


#########################################################
##      END PACKAGE RTools4TB
#########################################################
