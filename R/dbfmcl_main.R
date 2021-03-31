################################################################
##        Main script of the R PACKAGE : DBFMCL
##
## Authors : BERGON A, J. BAVAIS
##  with the collaboration of LOPEZ F., TEXTORIS J. and PUTHIER D.
##
##
## R CMD SHLIB dbf.c -o dbf
##
#################################################################
##    DEFINITION OF A SPECIFIC CLASS OBJECT : DBFMCLresult
#################################################################
library(Biobase)
library(ggplot2)
library(reshape2)

setClass("DBFMCLresult",

   representation=list(
       name="character",
       data="matrix",
       cluster="vector",
       size="vector",
       center="matrix",
       parameters="list"
   ),
   prototype=list(
       name=character(),
       data=matrix(nr=0, nc=0),
       cluster=numeric(),
       size=numeric(),
       center=matrix(nc=0, nr=0),
       parameters=list()
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
##    Define the plot_profile for class DBFMCLresult
#################################################################

setGeneric("plot_profile", function(object, type=c("line", "tile")) {
  standardGeneric("plot_profile")
})

setMethod("plot_profile", signature(object = "DBFMCLresult"), function(object, type=c("line", "tile")) {

        ## getting matrix
        nb <- length(object@size)
        m <- object@data

        ## median-centering of row
        med <- apply(m, 1, median)
        mc <- m - med

        ##  hclust on samples
        dis <- dist(t(mc))
        h <- hclust(dis, meth="av")
        m <- mc[, h$order]
        myOut <- ""


        ## melting

        m_melt <- as.data.frame(m)
        m_melt$cluster <- object@cluster
        m_melt$gene <- row.names(object@data)

        m_melt <- melt(m_melt,
                       id.vars=c("cluster", "gene"),
                       variable.name="samples")


        ## ploting

        if(type == "line"){
        p <- ggplot(data=m_melt, aes(x=samples,
                                     y=value))
        p <- p + geom_line(color="gray", aes(group=gene))
        p <- p + facet_grid(cluster~.)
        p <- p +  geom_line(data=m_melt %>%
                             group_by(cluster, samples)  %>%
                             summarise(cluster_mean=mean(value)),
                        aes(x = samples, y = cluster_mean, group=cluster),
                        color="black")
        }else if(type == "tile"){
            p <- ggplot(data=m_melt, aes(x=samples,
                             y=gene, fill=value))
            p <- p + geom_tile()
        }

        p <- p + facet_grid(cluster~., scales = "free_y")
        p <- p + theme_bw()
        p <- p + theme(strip.text.y = element_text(angle=0),
                       axis.text.x = element_text(angle=90),
                       axis.text.y=element_blank())

        return(p)
})



#################################################################
##    Define the write function for class DBFMCLresult
#################################################################

setGeneric("write_clusters", function(object, filename.out=NULL, path=".",
                verbose=TRUE) {
  standardGeneric("write_clusters")
})

setMethod("write_clusters", signature(object = "DBFMCLresult"), function(object, filename.out=NULL, path=".",
                verbose=TRUE) {
    if(path == ".") path <- getwd()


    if(is.null(filename.out))
      filename.out <- "exprs.dataMods.txt"

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


})



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



#################################################################
##    Getting data
#################################################################

getData4DBFMCL <- function(data=NULL, filename=NULL, path="."){

    ## getting matrix (probesID vs SamplesID)
    if(!is.null(data)){
        if(inherits(data, "Seurat")){
            data <- as.matrix(data@assays$RNA@data)
        }
        else if(is.data.frame(data)){
            data <- as.matrix(data)
        }
        if(!is.matrix(data)){
            stop("Please provide a Seurat Object, a data.frame",
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
#########################################################
##      END PACKAGE DBFMCL
#########################################################
