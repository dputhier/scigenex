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
library(dplyr)

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

setGeneric("plot_dbf", function(object,
                                    type=c("line", "tile"),
                                    to_log2=TRUE,
                                    color_palette=NULL,
                                    standardizing=FALSE,
                                    centering=TRUE) {
  standardGeneric("plot_dbf")
})

setMethod("plot_dbf", signature(object = "DBFMCLresult"), function(object,
                                                                       type=c("line", "tile"),
                                                                       to_log2=TRUE,
                                                                       color_palette="#0000BF,#0000FF,#0080FF,#00FFFF,#40FFBF,#80FF80,#BFFF40,#FFFF00,#FF8000,#FF0000,#BF0000",

                                                                       standardizing=FALSE,
                                                                       centering=TRUE) {

        # The type of diagram
        type = type[1]


        # color_palette_list = color_palette.split(",")
        # if len(color_palette_list) < 2:
        #    message("Need more than 2 colors for heatmap color palette.",
        #            type="ERROR")

        ## getting matrix
        nb <- length(object@size)
        m <- object@data

        if(to_log2){
            m <- log2(m + 1)
        }

        ## median-centering of row
        if (centering){
            mean_row <- apply(m, 1, mean)
            m <- sweep(m, MARGIN = 1, STATS=mean_row, FUN="-")
        }

            ## median-centering of row
        if (standardizing){
            sd_row <- apply(m, 1, sd)
            m <- sweep(m, MARGIN = 1, STATS=sd_row, FUN="/")
        }

        ##  hclust on samples
        dis <- dist(t(m))
        h <- hclust(dis, method="av")
        m <- m[, h$order]
        myOut <- ""


        ## melting

        m_melt <- as.data.frame(m)
        m_melt$cluster <- object@cluster
        m_melt$gene <- row.names(object@data)

        m_melt <- melt(m_melt,
                       id.vars=c("cluster", "gene"),
                       variable.name="samples")


        ## ploting

        col <- unlist(strsplit('#67001f,#b2182b,#d6604d,#f4a582,#fddbc7,#f7f7f7,#d1e5f0,#92c5de,#4393c3,#2166ac,#053061', ","))
        color.ramp <- colorRampPalette(col)(10)
        # Note that samples, value, gene, cluster
        # may appear as undefined variable to R check.
        # A workaround is to define them as NULL first...
        samples <- value <- gene <- cluster <- cluster_mean <- NULL
        if(type == "line"){
        p <- ggplot(data=m_melt, aes(x=samples,
                                     y=value))
        p <- p + geom_line(color="gray", aes(group=gene))
        p <- p + theme_bw()
        p <- p + facet_grid(cluster~.)
        p <- p +  geom_line(data=m_melt %>%
                                    group_by(cluster, samples)  %>%
                                        summarise(cluster_mean=mean(value)),
                            aes(x = samples,
                                y = cluster_mean,
                                group=cluster),
                            color="black")

        p <- p + theme(strip.text.y = element_text(angle=0),
                       axis.text.x = element_blank())

        }else if(type == "tile"){
            p <- ggplot(data=m_melt,
                        aes(x=samples,
                            y=gene,
                            fill=value))

            p <- p + geom_tile()
            p <- p + theme_bw()
            p <- p + scale_fill_gradientn(colours = color.ramp,
                              name="Signal")
            p <- p + theme(strip.text.y = element_text(angle=0),
                           axis.text.x  = element_blank(),
                           axis.text.y  = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.ticks.x = element_blank())
        }

        p <- p + facet_grid(cluster~., scales = "free_y")

        return(p)
})

#################################################################
##    Define the write function for class DBFMCLresult
#################################################################

setGeneric("write_dbf", function(object, filename.out=NULL, path=".",
                verbose=TRUE) {
  standardGeneric("write_dbf")
})

setMethod("write_dbf", signature(object = "DBFMCLresult"), function(object, filename.out=NULL, path=".",
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
            intLine <- matrix(c("NA", rep(NA,ncol(data))), nrow=1)
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



#' The "Density Based Filtering and Markov CLustering" algorithm (DBF-MCL).
#' 
#' DBF-MCL is a tree-steps adaptative algorithm that \emph{(i)}find elements
#' located in dense areas (DBF), \emph{(ii)}uses selected items to construct a
#' graph, \emph{(iii)}performs graph partitioning using the Markov CLustering
#' Algorithm (MCL).
#' 
#' This function requires installation of the mcl program
#' (\url{http://www.micans.org/mcl}). See "Warnings" section for more
#' informations.
#' 
#' When analyzing a noisy dataset, one is interested in isolating dense regions
#' as they are populated with genes/elements that display weak distances to
#' their nearest neighbors (i.e. strong profile similarities). To isolate these
#' regions DBF-MCL computes, for each gene/element, the distance with its kth
#' nearest neighbor (DKNN).In order to define a critical DKNN value that will
#' depend on the dataset and below which a gene/element will be considered as
#' falling in a dense area, DBF-MCL computes simulated DKNN values by using an
#' empirical randomization procedure. Given a dataset containing n genes and p
#' samples, a simulated DKNN value is obtained by sampling n distance values
#' from the gene-gene distance matrix D and by extracting the kth-smallest
#' value. This procedure is repeated n times to obtain a set of simulated DKNN
#' values S. Computed distributions of simulated DKNN are used to compute a FDR
#' value for each observed DKNN value. The critical value of DKNN is the one
#' for which a user-defined FDR value (typically 10\%) is observed. Genes with
#' DKNN value below this threshold are selected and used to construct a graph.
#' In this graph, edges are constructed between two genes (nodes) if one of
#' them belongs to the k-nearest neighbors of the other. Edges are weighted
#' based on the respective coefficient of correlation (\emph{i.e.}, similarity)
#' and the graph obtained is partitioned using the Markov CLustering Algorithm
#' (MCL).
#' 
#' @param data a \code{matrix}, \code{data.frame} or \code{Seurat} object.
#' @param filename a character string representing the file name.
#' @param name a prefix for the names of the intermediary files created by DBF
#' and MCL.
#' @param path a character string representing the data directory where
#' intermediary files are to be stored. Default to current working directory.
#' @param distance.method a method to compute the distance to the k-th nearest
#' neighbor. One of "pearson" (Pearson's correlation coefficient-based
#' distance), "spearman" (Spearman's rho-based distance), "euclidean". The
#' "spm" distance corresponds to the arithmetic mean :("pearson"+"spearman")/2
#' whereas "spgm" is the geometric mean : sqrt("pearson"*"spearman).
#' @param clustering indicates whether partitioning step (MCL) should be
#' applied to the data. If \code{clustering = FALSE}, the function returns a
#' \code{DBFMCLresult} object that contains informative elements (as detected
#' by the DBF step) coerced into a single cluster.
#' @param silent if set to TRUE, the progression of distance matrix calculation
#' is not displayed.
#' @param verbose if set to TRUE the function runs verbosely.
#' @param k the neighborhood size.
#' @param random the number of simulated distributions S to compute. By default
#' \code{random = 3}.
#' @param memory.used size of the memory used to store part of the distance
#' matrix. The subsequent sub-matrix is used to computed simulated distances to
#' the k-th nearest neighbor (see detail section).
#' @param fdr an integer value corresponding to the false discovery rate
#' (range: 0 to 100).
#' @param inflation the main control of MCL. Inflation affects cluster
#' granularity. It is usually chosen somewhere in the range \code{[1.2-5.0]}.
#' \code{inflation = 5.0} will tend to result in fine-grained clusterings, and
#' whereas \code{inflation = 1.2} will tend to result in very coarse grained
#' clusterings. By default, \code{inflation = 2.0}. Default setting gives very
#' good results for microarray data when k is set between 70 and 250.
#' @param set.seed specify seeds for random number generator.
#' @return a DBFMCLresults class object.
#' @section Warnings: With the current implementation, this function only works
#' only on UNIX-like plateforms.
#' 
#' MCL should be installed. One can used the following command lines in a
#' terminal:
#' 
#' \code{# Download the latest version of mcl (the script has been tested
#' successfully with the 06-058 version).}
#' \code{wget http://micans.org/mcl/src/mcl-latest.tar.gz}
#' \code{# Uncompress and install mcl}
#' \code{tar xvfz mcl-latest.tar.gz}
#' \code{cd mcl-xx-xxx}
#' \code{./configure}
#' \code{make}
#' \code{sudo make install}
#' \code{# You should get mcl in your path}
#' \code{mcl -h}
#' @author Bergon A., Bavais J., Textoris J., Granjeaud S., Lopez F and Puthier
#' D.
#' @seealso \code{\link{createSignatures4TB}}
#' @references Van Dongen S. (2000) A cluster algorithm for graphs. National
#' Research Institute for Mathematics and Computer Science in the 1386-3681.
#' @keywords manip
#' @examples
#' 
#' \dontrun{
#' ## with an artificial dataset
#' 
#' m <- matrix(rnorm(80000), nc=20)
#' m[1:100,1:10] <- m[1:100,1:10] + 4
#' m[101:200,11:20] <- m[101:200,11:20] + 3
#' m[201:300,5:15] <- m[201:300,5:15] + -2
#' plot_profile(res)
#' write_clusters(res, filename.out="ALL.sign.txt")
#' 
#' }
#' 
#' @export DBFMCL
DBFMCL <- function(data=NULL, filename=NULL, path=".", name=NULL,
           distance.method=c("pearson", "spearman", "euclidean", "spm", "spgm"),
           clustering=TRUE, silent=FALSE, verbose=TRUE, k=150,
           random=3, memory.used=1024, fdr=10, inflation=2.0,
           set.seed=123){

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
                    set.seed=set.seed)
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

                    obj@name <- name
                    obj@data <- as.matrix(m[sondeList, ])
                    obj@cluster <- clusters
                    obj@size <- size

                    centers <- matrix(ncol=ncol(m), nrow=nb)
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


#' Density-Based Filtering.
#' 
#' This function is an internal function used by \code{\link{DBFMCL}} to detect
#' informative elements (\emph{i.e.}, those that belong to dense regions). User
#' should not use this function. Instead they can use the \code{\link{DBFMCL}}
#' function with \code{clustering} argument set to \code{FALSE}.
#' 
#' See \code{\link{DBFMCL}}
#' 
#' @param data a matrix or data.frame
#' @param name a prefix for the file name
#' @param distance.method a method to compute the distance to the k-th nearest
#' neighbor. One of "pearson" (Pearson's correlation coefficient-based
#' distance), "spearman" (Spearman's rho-based distance) or "euclidean".
#' @param silent if set to TRUE (default), the progression of distance matrix
#' calculation is not displayed.
#' @param k the neighborhood size.
#' @param random the number of simulated distributions S to compute. By default
#' \code{random = FALSE}.
#' @param fdr a value for the false discovery rate.
#' @param memory.used size of the memory used to store part of the distance
#' matrix. The subsequent sub-matrix is used to computed simulated distances to
#' the k-th nearest neighbor (see detail section).
#' @param set.seed specify seeds for random number generator.
#' @section Warnings: Works only on UNIX-alikes platforms.
#' @author Bergon A., Bavais J., Textoris J., Granjeaud S., Lopez F and Puthier
#' D.
#' @seealso \code{\link{DBFMCL}}, \code{\link{createSignatures4TB}}
#' @references Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud
#' S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful
#' and flexible toolbox to explore productively the transcriptional landscape
#' of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.
#' @keywords manip
#' @export DBF
DBF <- function (data, name=NULL,
            distance.method=c("spearman", "pearson", "euclidean","spm", "spgm"),
            silent=FALSE,
            k=100,
            random=3,
            fdr=10,
            memory.used=1024,
            set.seed=123){

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


                obj@data <- as.matrix(data[a[[4]]%in%informative, ])
                obj@cluster <- rep(1, nrow(obj@data))
                obj@size <- nrow(obj@data)
                obj@center <- matrix(
                                    apply(obj@data[obj@cluster == 1, ],
                                        2,
                                        mean,
                                        na.rm=TRUE),
                                    nrow=1)
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

#' Invokes the Markov CLustering algorithm (MCL).
#' 
#' This function invokes the mcl system command. MCL is a clustering algorithm
#' for graphs that was developped by Stijn van Dongen (see references for
#' further informations).
#' 
#' 
#' @param name a character string corresponding to the file name.
#' @param inflation the main control of MCL. Inflation affects cluster
#' granularity. It is usually chosen somewhere in the range \code{[1.2-5.0]}.
#' \code{inflation = 5.0} will tend to result in fine-grained clusterings, and
#' whereas \code{inflation = 1.2} will tend to result in very coarse grained
#' clusterings. By default, \code{inflation = 2.0}. Default setting gives very
#' good results for microarray data when k is set around 100.
#' @param silent if set to TRUE, the progression of the MCL partitionning is
#' not displayed.
#' @return Returns a file with the ".mcl\_out.txt" extension.
#' @section warning: Works only on UNIX-like plateforms. MCL should be
#' installed. The following command lines can be used for installation.
#' 
#' \code{# Download the latest version of mcl (RTools4TB has been tested
#' successfully with the 06-058 version).}
#' \code{wget http://micans.org/mcl/src/mcl-latest.tar.gz}
#' \code{# Uncompress and install mcl}
#' \code{tar xvfz mcl-latest.tar.gz}
#' \code{cd mcl-xx-xxx}
#' \code{./configure}
#' \code{make}
#' \code{sudo make install}
#' \code{# You should get mcl in your path}
#' \code{mcl -h}
#' @author Bergon A., Lopez F., Textoris J., Granjeaud S. and Puthier D.
#' @references Stijn van Dongen. A cluster algorithm for graphs.  Technical
#' Report INS-R0010, National Research Institute for Mathematics and Computer
#' Science in the Netherlands, Amsterdam, May 2000.
#' \url{http://www.cwi.nl/ftp/CWIreports/INS/INS-R0010.ps.Z}
#' @keywords manip
#' @export MCL
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
                        sep="\t", header=TRUE, row.names=1, quote=""))
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
