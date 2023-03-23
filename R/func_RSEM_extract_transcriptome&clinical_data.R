#' @title RESM_preview
#' 
#' @description Cast a glance at the type of disease-samples.
#' 
#' @usage RESM_preview()
#' 
#' @param prefix_read The path of Xena data "TcgaTargetGTEX_phenotype.txt.gz".
#' 
#' @return A data.frame containing the type and number of diseases.
#' 
#' @export
#' 
#' @example RESM_preview()
#' 
#' @author GM.W
#' 
RESM_preview <- function(prefix_read = "D:/BaiduSyncdisk/Public_Database/TCGA_info/raw_xena/")
{
  # A function to preview the data category
  path_ph = paste0(prefix_read,"TcgaTargetGTEX_phenotype.txt.gz")
  ph = data.table::fread(path_ph)
  ph = ph[ph$`_study` %in% c("GTEX","TCGA"),]
  cate = as.data.frame(table(ph$`_primary_site`))
  return(cate)
}


#' @title RSEM_extract_whole
#' 
#' @description Extract and save the data of selected disease type from Xena data.
#' 
#' @usage RSEM_extract_whole(project = "STAD",
#' keep_cate = "stomach")
#' 
#' @param project The name of your project, which will be the prefix of the saved data.
#' @param keep_cate The type of disease. Use `RESM_preview()` to Cast a glance.
#' @param prefix_read The path your "TcgaTargetGTEX_phenotype.txt.gz" file stores.
#' @param prefix_out The path to save the extract data.
#' 
#' @return NULL
#' 
#' @export
#' 
#' @example RSEM_extract_whole(project = "STAD",
#' keep_cate = "stomach")
#' 
#' @author GM.W
#' 
RSEM_extract_whole <- function(project, keep_cate,
                               prefix_read="D:/BaiduSyncdisk/Public_Database/TCGA_info/raw_xena/",
                               prefix_out="D:/BaiduSyncdisk/Public_Database/TCGA_info/processed/")
{
  # A main function to extract certain transcriptome&clinical data from the RSEM dataset.
  out_file = paste0(prefix_out,project,"_transcriptome_whole_TCGA&GTEX",".Rdata")
  extract_data <- function(project=parent.frame()$project, keep_cate=parent.frame()$keep_cate,
                           prefix_read=parent.frame()$prefix_read, out_file=parent.frame()$out_file)
  {
    # A function to extract data
    cat(crayon::blue(">>> Processing the data:","\n"))
    count = data.table::fread(paste0(prefix_read,"TcgaTargetGtex_gene_expected_count.gz"),data.table = F)
    rownames(count) = count$sample
    count = count[,-1]

    fpkm = data.table::fread(paste0(prefix_read,"TcgaTargetGtex_rsem_gene_fpkm.gz"),data.table = F)
    rownames(fpkm) = fpkm$sample
    fpkm = fpkm[,-1]

    ph = data.table::fread(paste0(prefix_read,"TcgaTargetGTEX_phenotype.txt.gz"),data.table = F)

    survival_data = data.table::fread(paste0(prefix_read,"TCGA_survival_data"),data.table = F)

    probe_map = data.table::fread(paste0(prefix_read,"probeMap_gencode.v23.annotation.gene.probemap"),data.table = F)

    ph = ph[ph$`_study` %in% c("GTEX","TCGA"),]
    keep = ph[ph$`_primary_site` %in% keep_cate,1]
    # 2^exprset - b
    count = count[,intersect(keep,colnames(count))]
    count = round(2^count - 1)
    fpkm = fpkm[,intersect(keep,colnames(fpkm))]
    fpkm = 2^fpkm - 0.001
    fpkm[fpkm < 0] = 0
    ph = ph[ph$sample %in% keep,]
    survival_data = survival_data[survival_data$sample %in% keep,]
    count = IOBR::anno_eset(eset = count,annotation = probe_map,symbol = "gene",
                            probe = "id",method = "mean")
    fpkm = IOBR::anno_eset(eset = fpkm,annotation = probe_map,symbol = "gene",
                           probe = "id",method = "mean")
    save(count,fpkm,ph,survival_data,file = out_file)
    cat(crayon::blue(paste0(">>> The data of project: ",project,"have been successfully extracted!"), "\n"))
  }

  if(file.exists(out_file))
  {
    cat(crayon::blue(">>> The processed file is already exist!", "\n"))
    if(file.size(out_file) < 1e6)
    {
      cat(crayon::blue("But the file size is too small, re-extract now!", "\n"))
      extract_data()
    }else{
      cat(crayon::blue(">>> The processed file seems to have a normal size, terminate this procedure automatically."))
    }
  }else
  {
    extract_data()
  }

}

