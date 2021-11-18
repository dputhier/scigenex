#################################################################
##    Define the write function for class ClusterSet
#################################################################


#' @title
#' write_clust
#' @description
#' Write a ClusterSet into a flat file.
#' @param object ClusterSet. 
#' @param filename_out The outfile name.
#' @param out_path The path to the file.
#' @param nb_na_row Number of separating rows (containing NAs).
#' @return Write a file.
#' @export
#'
#' @examples
#' # see online help.
setGeneric("write_clust",
           
           function(object, 
                    filename_out = NULL,
                    out_path = ".",
                    nb_na_row=3) {
             standardGeneric("write_clust")
           })


#' @rdname write_clust
setMethod(
  "write_clust",
  signature(object = "ClusterSet"),
  function(object,
           filename_out = NULL,
           out_path = ".",
           nb_na_row=5) {
    
    if (out_path == ".") out_path <- getwd()
    
    if (is.null(filename_out)) {
      filename_out <- "exprs.dataMods.txt"
    }
    
    data <- object@data
    nb <- 0
    dataT <- c("clusters", colnames(data))
    
    ## processing data
    for (i in 1:length(object@size)) {
      
      print_msg(paste("Cluster ", i, " --> ", object@size[i], " probes"), msg_type="INFO")
      
      subData <- data[object@cluster == i, ]
      subData <- cbind(rownames(subData), subData)
      if (nb_na_row > 0){
        intLine <- matrix(rep(NA, (ncol(data) + 1)*nb_na_row), nrow = nb_na_row)
        dataT <- rbind(dataT, subData, intLine)
      }
      nb <- nb + 1
      
    }
    
    
    ## exporting results
    print_msg("Exporting results", msg_type="DEBUG")
    write.table(dataT, file.path(out_path, filename_out),
                col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE
    )
    
    print_msg(paste("--> Creating file : ",
                    file.path(out_path, filename_out)),
              msg_type="DEBUG")
  }
)