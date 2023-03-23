#' @title out_image
#' 
#' @description Export plot using either "export" or "ggsave".
#' 
#' @usage out_image(x=NULL,
#' filename = filename,
#' width = width,
#' height = height,
#' method = "export")
#' 
#' @param x The plot. Default is `NULL` (`last_plot()`).
#' @param filename The name of the exporting figure.
#' @param width The width of figure.
#' @param height The height of figure.
#' @param method The method used to export the plot. Choose "export"(default) or "ggsave".
#' 
#' @return `NULL`
#' 
#' @export
#' 
#' @author GM.W
#' 
out_image <- function(x=NULL,filename,width,height,method="export"){
  if(method=="export"){
    require(export)
    graph2tif(x=x,file=filename,width=width,height=height)
    graph2pdf(x=x,file=filename,width=width,height=height)
  }

  if(method=="ggsave"){
    require(ggplot2)
    if(is.null(x)){x=last_plot()}
    ggsave(plot=x,filename=paste0(filename,".tiff"),width = width,height = height)
    ggsave(plot=x,filename=paste0(filename,".pdf"),width = width,height = height,device = "pdf")
  }
}

#' @title out_plots
#' 
#' @description Export plot using either "export" or "ggsave" to the "plots" dir.
#' 
#' @usage out_image(x=NULL,
#' filename = filename,
#' width = width,
#' height = height,
#' method = "export")
#' 
#' @param x The plot. Default is `NULL` (`last_plot()`).
#' @param filename The name of the exporting figure. Then `filename = paste0("plots/", filename)`.
#' @param width The width of figure.
#' @param height The height of figure.
#' @param method The method used to export the plot. Choose "export"(default) or "ggsave".
#' 
#' @return `NULL`
#' 
#' @export
#' 
#' @author GM.W
#' 
out_plots <- function(x=NULL, filename, width, height, method = "export"){
  if(!dir.exists("/plots")){dir.create("/plots")}
  out_image(x=x, filename = paste0("plots/", filename),
            width = width, height = height, method = method)
}