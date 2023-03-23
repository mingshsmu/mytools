#' set Class
setClass(Class = "enrich",
         slots = list(DEG="data.frame", dys_gene="data.frame",enrich="list",
                      GSEA="list",para="list"))
#' @title enrich_analysis
#' 
#' @description A function to conduct enrichment analysis including Hypergeometric distribution and GSEA.
#' 
#' @usage enrich_analysis(DEG_res, 
#' n_top = 200,
#' method = "DESeq2",
#' signatures = c("GO","KEGG","HALLMARK","REACTOME","TF","miRNA"),
#' padj_threshold = 0.05)
#' 
#' @param DEG_res The result of DEG_calculate
#' @param n_top The top n genes were selected for Hypergeometric distribution (GO & KEGG).
#' @param method The method of `DEG_calculate`. It should be one of c("DESeq2","edgeR","limma_voom","limma","wilcox").
#' @param signatures The signatures for GSEA. Default is all of c("GO","KEGG","HALLMARK","REACTOME","TF","miRNA").
#' @param padj_threshold The p-adj threshold of DEGs. 
#' 
#' @return A object of "enrich" containing the results of Hypergeometric distribution and GSEA, as well as parameters used.
#' 
#' @export
#' 
#' @example enrich_analysis(DEG_res = DEG, n_top = 200, method = "DESeq2",
#' signatures=c("HALLMARK"),
#' padj_threshold=0.05)
#' 
#' @author GM.W
#' 
enrich_analysis <- function(DEG_res, n_top=200,method=c("DESeq2","edgeR","limma_voom","limma","wilcox"),
                            signatures=c("GO","KEGG","HALLMARK","REACTOME","TF","miRNA"),
                            padj_threshold=0.05){
  # 1. check parameter
  if(!all(signatures %in% c("GO","KEGG","HALLMARK","REACTOME","TF","miRNA"))){
    stop("Error: parameter error! The signature names you provide are not right!")
  }
  if(!is.integer(n_top)){tn_top <- as.integer(n_top)}
  method <- method[1]
  if(!method %in% c("DESeq2","edgeR","limma_voom","limma","wilcox")){
    stop("Error: parameter error! The method you provide is not allowed!")
  }
  
  # 2. preprocessing
  log2fc <- switch (method,
                    "DESeq2" = "log2FoldChange",
                    "edgeR" = "logFC",
                    "limma_voom" = "logFC",
                    "limma" = "logFC",
                    "wilcox" = "logFC"
  )
  
  padj <- switch (method,
                  "DESeq2" = "padj",
                  "edgeR" = "FDR",
                  "limma_voom" = "adj.P.Val",
                  "limma" = "adj.P.Val",
                  "wilcox" = "P.Value"
  )
  cat(crayon::green(">>> Performing enrichment analysis...","\n"))
  # 3. enrichment analysis
  require(clusterProfiler)
  require(org.Hs.eg.db)
  require(dplyr)
  dys_gene <- DEG_res[which(DEG_res[,padj] < padj_threshold),]
  dys_gene$entrez <- mapIds(x=org.Hs.eg.db, keys = rownames(dys_gene),
                            column = "ENTREZID", keytype = "SYMBOL")
  dys_gene <- na.omit(dys_gene)
  dys_gene <- dys_gene[order(dys_gene[,log2fc],decreasing = T),]
  gene_up <- head(dys_gene, n_top)$entrez
  gene_down <- tail(dys_gene, n_top)$entrez
  gene_dys <- unique(c(gene_up,gene_down))
  # GO
  cat(crayon::green(">>> GO enrichment analysis...","\n"))
  GO.up <- enrichGO(gene_up,OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID", ont = "ALL")
  GO.up <- setReadable(GO.up, org.Hs.eg.db)
  GO.up <- GO.up@result
  GO.down <- enrichGO(gene_down,OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID", ont = "ALL")
  GO.down <- setReadable(GO.down, org.Hs.eg.db)
  GO.down <- GO.down@result
  GO.all <- enrichGO(gene_dys,OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID", ont = "ALL")
  GO.all <- setReadable(GO.all, org.Hs.eg.db)
  GO.all <- GO.all@result
  GO.list <- list(GO.up, GO.down, GO.all)
  names(GO.list) <- c("GO.up","GO.down","GO.all")
  # KEGG
  cat(crayon::green(">>> KEGG enrichment analysis...","\n"))
  R.utils::setOption("clusterProfiler.download.method",'auto')
  KEGG.up <- enrichKEGG(gene_up,organism = "hsa",keyType = "kegg")
  KEGG.up <- setReadable(KEGG.up,OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID")
  KEGG.up <- KEGG.up@result
  KEGG.down <- enrichKEGG(gene_down,organism = "hsa",keyType = "kegg")
  KEGG.down <- setReadable(KEGG.down,OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID")
  KEGG.down <- KEGG.down@result
  KEGG.all <- enrichKEGG(gene_dys,organism = "hsa",keyType = "kegg")
  KEGG.all <- setReadable(KEGG.all,OrgDb = org.Hs.eg.db,
                           keyType = "ENTREZID")
  KEGG.all <- KEGG.all@result
  KEGG.list <- list(KEGG.up, KEGG.down, KEGG.all)
  names(KEGG.list) <- c("KEGG.up","KEGG.down","KEGG.all")
  # 4. GSEA
  gene_list <- DEG_res[,log2fc]
  names(gene_list) <- rownames(DEG_res)
  gene_list <- sort(gene_list,decreasing = T)
  gene_list <- na.omit(gene_list)
  GSEA_res <- lapply(signatures, function(x){
    cat(crayon::green(">>> GSEA",x,"...","\n"))
    path_dir <- "D:/BaiduSyncdisk/Public_Database/GSEA_signature/"
    path_sig <- switch(x,
                       "GO" = "c5.go.v2022.1.Hs.symbols_GO.gmt",
                       "KEGG" = "c2.cp.kegg.v2022.1.Hs.symbols_KEGG.gmt",
                       "HALLMARK" = "h.all.v2022.1.Hs.symbols_hallmark.gmt",
                       "REACTOME" = "c2.cp.reactome.v2022.1.Hs.symbols.gmt",
                       "TF" = "c3.tft.v2022.1.Hs.symbols_TF.gmt",
                       "miRNA" = "c3.mir.v2022.1.Hs.symbols_miRNA.gmt")
    term2gene <- read.gmt(paste0(path_dir, path_sig))
    gsea.res <- GSEA(geneList = gene_list, TERM2GENE = term2gene)
    return(gsea.res)
  })
  names(GSEA_res) <- signatures
  
  # 5. object
  object <- new("enrich",DEG = DEG_res,dys_gene = dys_gene,
                enrich = list(GO = GO.list,
                              KEGG = KEGG.list),
                GSEA = GSEA_res,
                para = list(n_top = n_top, method = method,
                            signatures = signatures, padj_threshold = padj_threshold))
  cat(crayon::green(">>> Done!","\n"))
  return(object)
}


#' @title plot_enrich
#' 
#' @description A function to visualize the results of GO and KEGG enrichment analysis.
#' 
#' @usage plot_enrich(enrich_obj, 
#' type=c("dot","bar","bubble"),
#' enrich=c("GO","KEGG"),
#' facet=FALSE,
#' regulate=c("up","down","dys","wrap"),
#' max_term=18,
#' title=NULL,
#' color=NULL)
#' 
#' @param enrich_obj The result of `enrich_analysis`.
#' @param type The type of plot.
#' @param enrich Either "GO" or "KEGG".
#' @param facet Logical.
#' @param regulation The regulation of genes used in enrichment analysis.
#' @param max_term The maximum number of terms to draw.
#' @param title The title of plot.
#' @param color The color used in plot. Default is "Set1".
#' 
#' @return A ggplot object.
#' 
#' @export
#' 
#' @examples plot_enrich(enrich_obj = enrich_obj,type = "dot",
#' enrich="GO",regulate = "up",title = "GO.up")
#' 
#' @author GM.W
#' 
plot_enrich <- function(enrich_obj, type=c("dot","bar","bubble"),enrich=c("GO","KEGG"),facet=FALSE,
                        regulate=c("up","down","dys","wrap"),max_term=18,title=NULL,color=NULL){
  # 1.Initialization
  type <- type[1]
  enrich <- enrich[1]
  regulate <- regulate[1]
  require(ggplot2)
  if(is.null(color)){
    color <- RColorBrewer::brewer.pal(9,"Set1")
  }
  
  #' @title prepare_enrich_df
  #' 
  prepare_enrich_df <- function(i,enrich_obj=parent.frame()$enrich_obj,enrich=parent.frame()$enrich,max_term=parent.frame()$max_term){
    enrich.df <- enrich_obj@enrich[[enrich]][[i]]
    enrich.df <- enrich.df[order(enrich.df$p.adjust),]
    if(enrich=="KEGG"){
      kegg.class <- read.csv("D:/BaiduSyncdisk/Public_Database/GSEA_signature/KEGG_pathway_classifications_ok.csv",row.names = 2)
      kegg.class <- kegg.class[enrich.df$ID,]
      enrich.df$ONTOLOGY <- kegg.class$Pathway.Class.1
    }
    enrich.df <- na.omit(enrich.df)
    n <- dim(enrich.df)[1]
    n <- ifelse(n > max_term,max_term,n)
    enrich.df <- enrich.df[1:n,]
    regulation <- switch (i,
      "1" = "Activated Genes",
      "2" = "Suppressed Genes",
      "3" = "Dysregulated Genes")
    enrich.df$regulation <- factor(regulation,levels = c("Activated Genes","Suppressed Genes","Dysregulated Genes"))
    # Ratio
    tmp <- stringr::str_split(enrich.df$GeneRatio,"/",simplify = T)
    enrich.df$GeneRatio <- round(as.numeric(tmp[,1])/as.numeric(tmp[,2]),2)
    return(enrich.df)
  }
  
  enrich.up <- prepare_enrich_df(1)
  enrich.down <- prepare_enrich_df(2)
  enrich.dys <- prepare_enrich_df(3)
  
  half_max_term <- round(max_term/2)
  third_max_term <- round(max_term/3)
  enrich.wrap <- rbind(enrich.up[1:half_max_term,],enrich.down[1:half_max_term,])
  enrich.wrap <- rbind(enrich.wrap,enrich.dys[1:half_max_term,])
  
  # factor
  #' @title factor_df
  #' 
  factor_df <- function(enrich.df){
    enrich.df$label <- NA
    enrich.df[1:third_max_term,"label"] <- enrich.df[1:third_max_term,"ID"]
    enrich.df <- dplyr::arrange(enrich.df,ONTOLOGY,desc(GeneRatio))
    enrich.df$Description <- factor(enrich.df$Description,levels = rev(enrich.df$Description))
    return(enrich.df)
  }
  enrich.up <- factor_df(enrich.up)
  enrich.down <- factor_df(enrich.down)
  enrich.dys <- factor_df(enrich.dys)
  
  plot_data <- switch (regulate,
                       "up" = enrich.up,
                       "down" = enrich.down,
                       "dys" = enrich.dys,
                       "wrap" = enrich.wrap
  )
  
  # 2. plot type
  if(regulate != "wrap"){
    p <- ggplot(data = plot_data)+
      labs(y="",title=title)+
      scale_color_manual(values = color)+
      scale_fill_manual(values = color)+
      theme_bw()+
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5))
      
    if(type == "dot"){ 
      p <- p + 
          geom_point(aes(x=GeneRatio,y=Description,color=ONTOLOGY,fill=ONTOLOGY,size=-log10(p.adjust)),shape=21,alpha=0.8)+
          scale_size(range = c(3,9))
    }else if(type=="bar"){
      p <- p +
          geom_bar(aes(x=GeneRatio,y=Description,color=ONTOLOGY,fill=ONTOLOGY),stat = "identity")
    }else if(type=="bubble"){
      p <- p+
          geom_jitter(aes(x=GeneRatio,y=-log10(p.adjust),color=ONTOLOGY,fill=ONTOLOGY,size=-log10(p.adjust)),shape=21,alpha=0.8)+
          scale_size(range = c(3,9))
    }
    if(facet==TRUE){
      p <- p+
          facet_wrap(.~ONTOLOGY)
    }
    if(enrich == "KEGG"){
      p <- p+
          guides(color=guide_legend(title = "Pathway Class"),
                 fill=guide_legend(title = "Pathway Class"))
    }
    
  }else{
    # regulation == "wrap"
    p <- ggplot(plot_data,aes(x=GeneRatio,y=Description,fill=regulation,color=regulation))+
      geom_point(aes(size=-log10(p.adjust)),shape=21,alpha=0.8)+
      scale_size(range = c(3,9))+
      labs(y="",title=title)+
      facet_wrap(.~regulation)+
      theme_bw()+
      theme(panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5))+
      scale_color_manual(values = color)+
      scale_fill_manual(values = color)
  }
 return(p)
}


#' @title plot_GSEA
#'
#' @description A function to visualize the results of GO and KEGG enrichment analysis.
#' 
#' @usage plot_GSEA(enrich_obj, 
#' type=c("GSEAcurve","bar","dot"),
#' geneSetID=NULL,
#' regulate=NULL, 
#' max_signatures=6, 
#' addGene=NULL, 
#' sig=c("GO","KEGG","HALLMARK","REACTOME","TF","miRNA"),
#' color=NULL, 
#' subplot=3, 
#' legend.position=NULL, 
#' title=NULL)
#' 
#' @param enrich_obj The result of `enrich_analysis`.
#' @param type The type of plot.
#' @param geneSetID The id of geneset.
#' @param regulate The regulation of signatures.
#' @param max_signatures The maximum number of signatures.
#' @param addGene The gene of interest to add on graph.
#' @param sig The signature choiced to visualize
#' @param color The color used for curve.
#' @param subplot The number of plots to draw.
#' @param legend.position The position of legend.
#' @param title The title of plot.
#' 
#' @return A plot.
#' 
#' @export
#' 
#' @example plot_GSEA(enrich_obj = enrich_obj,type = "GSEAcurve",
#' regulate = "pos",sig = "HALLMARK",subplot = 2,
#' max_signatures = 10)
#' 
#' @author GM.W
#' 
plot_GSEA <- function(enrich_obj, type=c("GSEAcurve","bar","dot"), geneSetID=NULL,
                      regulate=NULL, max_signatures=6, addGene=NULL, 
                      sig=c("GO","KEGG","HALLMARK","REACTOME","TF","miRNA"),
                      color=NULL, subplot=3, legend.position=NULL, title=NULL){
  # 1. Parameters
  type <- type[1]
  regulate <- regulate[1]
  sig <- sig[1]
  if(is.null(color)){
    color <- c("#D20A13","#088247","#FFD121","#223D6C","#91612D","#58CDD9",
               "#11AA4D","#223D6C","#E0367A","#7A142C","#5D90BA","#6E568C")
    }
  
  if(!is.null(geneSetID) & !is.null(regulate)){
    cat(crayon::green("Both params geneSetID and regulate are detected. Plot will be drawn accroding to geneSetID"))
    }
  
  ### GSEAcurve
  if(type == "GSEAcurve"){
    if(length(geneSetID) > 1 & !is.null(addGene)){
      stop("The param addGene should be passed only when the geneSetID was given 1.")
      }
    if(is.null(legend.position)){
      legend.position <- c(0.8, 0.8)
      }
      
  }
  
  ### bar and dot plot
  if(type %in% c("bar","dot")){
    if(!is.null(addGene)){
      cat(crayon::green("The params addGene and subplot will be not used unless the type is GSEAcurve."))
      }
    if(is.null(legend.position)){
      legend.position <- "bottom"
      }
  }
  
  # 2. Type
  obj <- enrich_obj@GSEA[[sig]]
  obj.df <- obj@result
  if(type=="GSEAcurve"){
    if(!is.null(regulate)){
      up.index <- which(obj.df$NES > 0)
      n <- length(up.index)
      n <- ifelse(n > max_signatures, max_signatures, n)
      up.index <- up.index[1:n]
      
      down.index <- which(obj.df$NES < 0)
      n <- length(down.index)
      n <- ifelse(n > max_signatures, max_signatures, n)
      down.index <- down.index[1:n]
      
      both.index <- c(up.index,down.index)
      if(regulate=="pos"){
        if(length(up.index)==1){
          p <- .gseaNb(object = obj,subPlot = subplot,geneSetID = up.index,
                               legend.position = legend.position)
        }else{
          curveCol <-  color[1:length(up.index)]
          p <- .gseaNb(object = obj,subPlot = subplot,geneSetID = up.index,
                               curveCol = curveCol,legend.position = legend.position)
        }
        
      }else if(regulate=="neg"){
        if(length(down.index)==1){
          p <- .gseaNb(object = obj,subPlot = subplot,geneSetID = down.index,
                               legend.position = legend.position)
        }else{
          curveCol <- color[(length(color)-length(down.index)+1):length(color)]
          p <- .gseaNb(object = obj,subPlot = subplot,geneSetID = down.index,
                               curveCol = curveCol,legend.position = legend.position)
        }
        
      }else if(regulate=="both"){
        curveCol <- c(color[1:length(up.index)],color[(length(color)-length(down.index)+1):length(color)])
        p <- .gseaNb(object = obj,subPlot = subplot,geneSetID = both.index,
                             curveCol = curveCol,legend.position = legend.position)
        
      }
      
      
    }
    
    if(!is.null(geneSetID)){
      if(length(geneSetID)==1){
        p <- .gseaNb(object = obj,geneSetID = geneSetID,subPlot = subplot,addGene = addGene)
      }else{
        p <- .gseaNb(object = obj,geneSetID = geneSetID,subPlot = subplot,
                        curveCol = color[1:length(geneSetID)],legend.position = legend.position)}
    }
  }
  
  if(type %in% c("bar","dot")){
    id_split <- stringr::str_split(obj.df$ID,"_",n=2,simplify = T)
    obj.df$term <- gsub("_"," ",id_split[,2])
    if(sig %in% c("TF","miRNA")){
      obj.df$term <- obj.df$ID
    }
    legend.title <- sig
    obj.df$regulate <- ifelse(obj.df$NES>0,"Positive","Negative")
    obj.df$regulate <- factor(obj.df$regulate, levels = c("Positive","Negative"))
    obj.df <- obj.df[order(obj.df$NES),]
    plot_data <- rbind(head(obj.df,max_signatures),tail(obj.df,max_signatures))
    plot_data <- plot_data[unique(plot_data$ID),]
    plot_data$term <- factor(plot_data$term,levels = unique(plot_data$term))
    plot_data_pos <- plot_data[plot_data$regulate=="Positive",]
    plot_data_neg <- plot_data[plot_data$regulate=="Negative",]
    
    require(ggplot2)
    if(type == "bar"){
      p <- ggplot(plot_data,aes(x=NES,y=term,fill=regulate),)+
        geom_bar(stat = "identity")+
        geom_text(data=plot_data_pos,aes(x=0,y=term,label=term),hjust=1,size=3)+
        geom_text(data=plot_data_neg,aes(x=0,y=term,label=term),hjust=0,size=3)+
        xlim((min(plot_data$NES)-2),(max(plot_data$NES)+2))+
        labs(y="",title = title)+
        theme_bw()+
        theme(legend.position = legend.position,
              panel.grid = element_blank(),
              plot.title = element_text(hjust = 0.5),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())+
        scale_fill_manual(values = color)+
        guides(fill=guide_legend(title = legend.title))
    }else if(type=="dot"){
      p <- ggplot(plot_data,aes(x=NES,y=term,fill=NES,color=NES),)+
        geom_point(aes(size=-log10(p.adjust)),shape=21)+
        labs(y="",title = title,subtitle = legend.title)+
        theme_bw()+
        theme(legend.position = legend.position,
              panel.grid = element_blank(),
              plot.title = element_text(hjust = 0.5))+
        scale_color_gradient2(high = color[1],low = color[2],mid = "grey")+
        scale_fill_gradient2(high = color[1],low = color[2],mid = "grey")+
        scale_size(range = c(4,10))
    }
    
  }
    
  # 3. return
  return(p)
}

#' @title .gseaNb
#' @description Plot GSEAcurves
#' 
.gseaNb <- function (object = NULL, subPlot = 3, lineSize = 0.8, geneSetID = NULL, 
                     rmSegment = FALSE, termWidth = 40, segCol = "red", addGene = NULL, 
                     geneCol = NULL, arrowAngle = 20, arrowLength = 0.2, arrowEnd = "first", 
                     arrowType = "closed", curveCol = c("#76BA99", "#EB4747", 
                                                        "#996699"), htCol = c("#08519C", "#A50F15"), rankCol = c("#08519C", 
                                                                                                                 "white", "#A50F15"), rankSeq = 5000, htHeight = 0.3, 
                     force = 20, max.overlaps = 50, geneSize = 4, newGsea = FALSE, 
                     addPoint = TRUE, newCurveCol = c("#336699", "white", "#993399"), 
                     newHtCol = c("#336699", "white", "#993399"), rmHt = FALSE, 
                     addPval = FALSE, pvalX = 0.9, pvalY = 0.9, pvalSize = 4, 
                     pCol = "grey30", pHjust = 1, rmPrefix = TRUE, nesDigit = 2, 
                     pDigit = 2, markTopgene = FALSE, topGeneN = 5, kegg = FALSE, 
                     legend.position = "right", add.geneExpHt = FALSE, exp = NULL, 
                     scale.exp = TRUE, sample.order = NULL, exp.col = c("blue", 
                                                                        "white", "red"), ht.legend = TRUE, ght.relHight = 0.4, 
                     ght.geneText.size = 6, ght.facet = FALSE, ght.facet.scale = "free", 
                     termID.order = NULL, rank.gene = NULL, rank.gene.nudgey = 2) 
{
  require(dplyr)
  gsdata <- purrr::map_df(geneSetID, function(setid) {
    GseaVis::gsInfo(object, geneSetID = setid) %>% dplyr::mutate(id = setid)
  })
  if (kegg == FALSE) {
    gsdata1 <- purrr::map_df(unique(gsdata$Description), 
                             function(setid) {
                               tmp <- gsdata %>% dplyr::filter(Description == 
                                                                 setid) %>% dplyr::mutate(gene_name = names(object@geneList)) %>% 
                                 dplyr::filter(position == 1)
                             })
  }
  else {
    gene2Symbol <- object@gene2Symbol %>% data.frame()
    gsdata1 <- purrr::map_df(unique(gsdata$Description), 
                             function(setid) {
                               tmp <- gsdata %>% dplyr::filter(Description == 
                                                                 setid) %>% dplyr::mutate(gene_name = gene2Symbol$.) %>% 
                                 dplyr::filter(position == 1)
                             })
  }
  data_ga <- data.frame(object) %>% dplyr::filter(ID %in% 
                                                    geneSetID)
  data_ga <- data_ga[unique(gsdata$id), ]
  niceTit <- purrr::map_chr(unique(gsdata$Description), function(x) {
    tit <- unlist(strsplit(x, split = "_"))
    if (length(tit) == 1) {
      niceTit <- paste(stringr::str_to_title(tit[1:length(tit)]), 
                       collapse = " ") %>% stringr::str_wrap(., width = termWidth)
    }
    else {
      if (rmPrefix == TRUE) {
        niceTit <- paste(tit[2:length(tit)], collapse = " ") %>% stringr::str_wrap(., width = termWidth)
      }
      else {
        niceTit <- paste(tit[1:length(tit)], collapse = " ") %>% stringr::str_wrap(., width = termWidth)
      }
    }
  })
  if (length(geneSetID) != 1) {
    ledend.t <- niceTit
    niceTit <- ""
  }
  if (length(geneSetID) == 1) {
    line <- ggplot2::geom_line(ggplot2::aes_(color = ~runningScore), 
                               size = lineSize)
    line.col <- ggplot2::scale_color_gradient(low = curveCol[1], 
                                              high = curveCol[2])
    legend.position = "none"
  }
  else {
    mulcol <- curveCol
    names(mulcol) <- unique(gsdata$Description)
    line <- ggplot2::geom_line(ggplot2::aes_(color = ~Description), 
                               size = lineSize)
    line.col <- ggplot2::scale_color_manual(values = mulcol, 
                                            labels = ledend.t, name = parent.frame()$sig)
    legend.position = legend.position
  }
  pcurve <- ggplot2::ggplot(gsdata, ggplot2::aes_(x = ~x, 
                                                  y = ~runningScore)) + line + line.col + ggplot2::geom_hline(yintercept = 0, 
                                                                                                              size = lineSize, color = "black", lty = "dashed") + 
    ggplot2::theme_bw(base_size = 14) + 
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme(legend.position = legend.position, 
                                                                                                    legend.box.background = ggplot2::element_blank(), plot.title = ggplot2::element_text(hjust = 0.5), 
                                                                                                    panel.grid = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), 
                                                                                                    axis.text.x = ggplot2::element_blank(), axis.line.x = ggplot2::element_blank(), 
                                                                                                    axis.title.x = ggplot2::element_blank(), legend.background = ggplot2::element_rect(fill = "transparent"), 
                                                                                                    plot.margin = ggplot2::margin(t = 0.2, r = 0.2, b = 0, 
                                                                                                                                  l = 0.2, unit = "cm")) + ggplot2::ylab("Running Enrichment Score") + 
    ggplot2::ggtitle(niceTit)
  midpoint <- sum(range(gsdata$runningScore))/2
  pnew <- ggplot2::ggplot(gsdata, ggplot2::aes_(x = ~x, y = ~runningScore, 
                                                color = ~runningScore)) + ggplot2::geom_hline(yintercept = 0, 
                                                                                              size = lineSize, color = "black", lty = "dashed") + 
    ggplot2::geom_line(size = lineSize) + ggplot2::geom_segment(data = gsdata1, 
                                                                ggplot2::aes_(xend = ~x, yend = 0)) + ggplot2::theme_bw(base_size = 14) + 
    ggplot2::scale_color_gradient2(low = newCurveCol[1], 
                                   mid = newCurveCol[2], high = newCurveCol[3], midpoint = midpoint) + 
    ggplot2::scale_x_continuous(expand = c(0, 0)) + ggplot2::theme(legend.position = "none", 
                                                                   plot.title = ggplot2::element_text(hjust = 0.5), axis.ticks.x = ggplot2::element_blank(), 
                                                                   axis.text.x = ggplot2::element_blank(), axis.line.x = ggplot2::element_blank(), 
                                                                   axis.title.x = ggplot2::element_blank(), legend.background = ggplot2::element_rect(fill = "transparent"), 
                                                                   plot.margin = ggplot2::margin(t = 0.2, r = 0.2, b = 0, 
                                                                                                 l = 0.2, unit = "cm")) + ggplot2::ylab("Running Enrichment Score") + 
    ggplot2::ggtitle(niceTit)
  if (addPoint == TRUE) {
    panother <- pnew + ggplot2::geom_point()
  }
  else {
    panother <- pnew
  }
  if (newGsea == FALSE) {
    pcurveRes <- pcurve
  }
  else {
    pcurveRes <- panother
  }
  if (is.null(addGene)) {
    plabel <- pcurveRes
  }
  else {
    if (markTopgene == TRUE) {
      geneLabel <- gsdata1 %>% dplyr::arrange(x) %>% dplyr::slice_head(n = topGeneN)
    }
    else {
      geneLabel <- gsdata1 %>% dplyr::filter(gene_name %in% 
                                               addGene)
    }
    if (nrow(geneLabel) == 0) {
      print("Your gene is not in this pathway! Please choose again!")
    }
    else {
      if (rmSegment == TRUE) {
        if (is.null(geneCol)) {
          plabel <- pcurveRes + ggrepel::geom_text_repel(data = geneLabel, 
                                                         ggplot2::aes_(label = ~gene_name), force = force, 
                                                         max.overlaps = max.overlaps, size = geneSize, 
                                                         fontface = "italic", arrow = ggplot2::arrow(angle = arrowAngle, 
                                                                                                     length = ggplot2::unit(arrowLength, "cm"), 
                                                                                                     ends = arrowEnd, type = arrowType))
        }
        else {
          plabel <- pcurveRes + ggrepel::geom_text_repel(data = geneLabel, 
                                                         ggplot2::aes_(label = ~gene_name), force = force, 
                                                         max.overlaps = max.overlaps, size = geneSize, 
                                                         fontface = "italic", color = geneCol, arrow = ggplot2::arrow(angle = arrowAngle, 
                                                                                                                      length = ggplot2::unit(arrowLength, "cm"), 
                                                                                                                      ends = arrowEnd, type = arrowType))
        }
      }
      else {
        if (is.null(geneCol)) {
          plabel <- pcurveRes + ggplot2::geom_segment(data = geneLabel, 
                                                      ggplot2::aes_(xend = ~x, yend = 0), color = segCol) + 
            ggrepel::geom_text_repel(data = geneLabel, 
                                     ggplot2::aes_(label = ~gene_name), force = force, 
                                     max.overlaps = max.overlaps, size = geneSize, 
                                     fontface = "italic", arrow = ggplot2::arrow(angle = arrowAngle, 
                                                                                 length = ggplot2::unit(arrowLength, 
                                                                                                        "cm"), ends = arrowEnd, type = arrowType))
        }
        else {
          plabel <- pcurveRes + ggplot2::geom_segment(data = geneLabel, 
                                                      ggplot2::aes_(xend = ~x, yend = 0), color = segCol) + 
            ggrepel::geom_text_repel(data = geneLabel, 
                                     ggplot2::aes_(label = ~gene_name), force = force, 
                                     max.overlaps = max.overlaps, size = geneSize, 
                                     fontface = "italic", color = geneCol, 
                                     arrow = ggplot2::arrow(angle = arrowAngle, 
                                                            length = ggplot2::unit(arrowLength, 
                                                                                   "cm"), ends = arrowEnd, type = arrowType))
        }
      }
    }
  }
  if (addPval == TRUE) {
    pLabel <- paste0("NES: ", round(data_ga$NES, digits = nesDigit), 
                     "\n", "Pvalue: ", ifelse(data_ga$pvalue < 0.001, 
                                              "< 0.001", round(data_ga$pvalue, digits = pDigit)), 
                     "\n", "Ajusted Pvalue: ", ifelse(data_ga$p.adjust < 
                                                        0.001, "< 0.001", round(data_ga$p.adjust, digits = pDigit)), 
                     "\n", sep = " ")
    px <- pvalX * nrow(gsdata[which(gsdata$id == geneSetID[1]), 
    ])
    py <- pvalY * sum(abs(range(gsdata$runningScore))) + 
      min(gsdata$runningScore)
    if (length(geneSetID) == 1) {
      pLabelOut <- plabel + ggplot2::annotate(geom = "text", 
                                              x = px, y = py, label = pLabel, size = pvalSize, 
                                              color = pCol, fontface = "italic", hjust = pHjust)
    }
    else {
      mytable <- tibble::tibble(x = px, y = py, table = list(tibble::tibble(NES = round(data_ga$NES, 
                                                                                        digits = nesDigit), Pvalue = ifelse(data_ga$pvalue < 
                                                                                                                              0.001, "< 0.001", round(data_ga$pvalue, digits = pDigit)), 
                                                                            `Ajusted Pvalue` = ifelse(data_ga$p.adjust < 
                                                                                                        0.001, "< 0.001", round(data_ga$p.adjust, 
                                                                                                                                digits = pDigit)))))
      pLabelOut <- plabel + ggpp::geom_table(data = mytable, 
                                             ggplot2::aes(px, py, label = table))
    }
  }
  else {
    pLabelOut <- plabel
  }
  if (length(geneSetID) == 1) {
    line.col <- ggplot2::scale_color_manual(values = "black")
  }
  if (add.geneExpHt == TRUE) {
    pseg.b = 0
  }
  else {
    pseg.b = 0.2
  }
  pseg <- ggplot2::ggplot(gsdata, ggplot2::aes_(x = ~x, y = ~runningScore, 
                                                color = ~Description)) + ggplot2::geom_segment(data = gsdata1, 
                                                                                               ggplot2::aes_(x = ~x, xend = ~x, y = 0, yend = 1), show.legend = F) + 
    line.col + ggplot2::scale_x_continuous(expand = c(0, 
                                                      0)) + ggplot2::scale_y_continuous(expand = c(0, 0)) + 
    ggplot2::theme_bw(base_size = 14) + ggplot2::theme(axis.ticks = ggplot2::element_blank(), 
                                                       axis.text = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(), 
                                                       panel.grid = ggplot2::element_blank(), axis.line.x = ggplot2::element_blank(), 
                                                       strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_blank(), 
                                                       panel.spacing = ggplot2::unit(0.1, "cm"), plot.margin = ggplot2::margin(t = 0, 
                                                                                                                               r = 0.2, b = pseg.b, l = 0.2, unit = "cm")) + ggplot2::xlab("Rank in Ordered Dataset") + 
    ggplot2::facet_wrap(~Description, ncol = 1)
  if (subPlot > 2) {
    pseg <- pseg + ggplot2::theme(axis.title.x = ggplot2::element_blank())
  }
  else {
    pseg <- pseg
  }
  d <- purrr::map_df(unique(gsdata$Description), function(setid) {
    tmp <- gsdata %>% dplyr::filter(Description == setid)
    v <- seq(1, sum(tmp$position), length.out = 9)
    inv <- findInterval(rev(cumsum(tmp$position)), v)
    if (min(inv) == 0) {
      inv <- inv + 1
    }
    color <- (grDevices::colorRampPalette(c(htCol[1], "white", 
                                            htCol[2])))(10)
    ymin <- 0
    yy <- htHeight
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, 
                    xmax = xmax, col = color[unique(inv)], Description = setid)
  })
  pseg_ht <- pseg + ggplot2::geom_rect(ggplot2::aes_(xmin = ~xmin, 
                                                     xmax = ~xmax, ymin = ~ymin, ymax = ~ymax, fill = ~I(col)), 
                                       data = d, alpha = 0.8, inherit.aes = FALSE)
  pseg_ht1 <- pseg_ht + ggplot2::xlab("") + ggplot2::theme(axis.title.x = ggplot2::element_blank(), 
                                                           plot.margin = ggplot2::margin(t = -0.1, r = 0.2, b = 0, 
                                                                                         l = 0.2, unit = "cm"))
  if (add.geneExpHt == TRUE) {
    prank.b = 0
  }
  else {
    prank.b = 0.2
  }
  prank <- ggplot2::ggplot(gsdata[which(gsdata$Description == 
                                          unique(gsdata$Description)[1]), ], ggplot2::aes_(x = ~x, 
                                                                                           y = ~geneList)) + ggplot2::geom_col(ggplot2::aes_(fill = ~geneList), 
                                                                                                                               width = 1, color = NA, show.legend = F) + ggplot2::scale_fill_gradient2(low = rankCol[1], 
                                                                                                                                                                                                       mid = rankCol[2], high = rankCol[3], midpoint = 0) + 
    ggplot2::geom_hline(yintercept = 0, size = 0.8, color = "black", 
                        lty = "dashed") + ggplot2::scale_x_continuous(breaks = seq(0, 
                                                                                   nrow(gsdata), rankSeq)) + ggplot2::theme_bw(base_size = 14) + 
    ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                   plot.margin = ggplot2::margin(t = -0.1, r = 0.2, 
                                                 b = prank.b, l = 0.2, unit = "cm")) + ggplot2::coord_cartesian(expand = 0) + 
    ggplot2::ylab("Ranked List") + ggplot2::xlab("Rank in Ordered Dataset")
  if (add.geneExpHt == TRUE) {
    prank <- prank + ggplot2::theme(axis.title.x = ggplot2::element_blank())
  }
  else {
    prank <- prank
  }
  if (kegg == FALSE) {
    rank.g <- data.frame(logfc = object@geneList, gene_name = names(object@geneList)) %>% 
      dplyr::mutate(x = 1:length(object@geneList))
  }
  else {
    rank.g <- data.frame(logfc = object@geneList, gene_name = object@gene2Symbol) %>% 
      dplyr::mutate(x = 1:length(object@geneList))
  }
  if (!is.null(rank.gene)) {
    target.rank.g <- rank.g %>% dplyr::filter(gene_name %in% 
                                                rank.gene) %>% dplyr::mutate(vjust = ifelse(logfc > 
                                                                                              0, "bottom", "top"), nudge_y = ifelse(logfc > 0, 
                                                                                                                                    -rank.gene.nudgey, rank.gene.nudgey))
    prank <- prank + ggrepel::geom_text_repel(data = target.rank.g, 
                                              ggplot2::aes(x = as.numeric(x), y = 0, label = gene_name, 
                                                           vjust = vjust, nudge_y = nudge_y), max.overlaps = 200, 
                                              direction = "x", angle = 90, fontface = "italic", 
                                              size = geneSize)
  }
  d <- purrr::map_df(unique(d$Description), function(x) {
    tmp <- d %>% dplyr::filter(Description == x)
    htcolor <- (grDevices::colorRampPalette(newHtCol))(nrow(tmp))
    tmp <- tmp %>% dplyr::mutate(htcol = htcolor)
  })
  ht <- ggplot2::ggplot(gsdata, ggplot2::aes_(x = ~x, y = ~runningScore)) + 
    ggplot2::geom_rect(ggplot2::aes_(xmin = ~xmin, xmax = ~xmax, 
                                     ymin = ~ymin, ymax = ~ymax, fill = ~I(htcol)), data = d, 
                       alpha = 0.8, inherit.aes = FALSE) + ggplot2::scale_x_continuous(expand = c(0, 
                                                                                                  0)) + ggplot2::scale_y_continuous(expand = c(0, 0)) + 
    ggplot2::theme_bw(base_size = 14) + ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                                                       axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), 
                                                       axis.title = ggplot2::element_blank(), plot.margin = ggplot2::margin(t = 0, 
                                                                                                                            r = 0.2, b = 0.2, l = 0.2, unit = "cm"))
  if (add.geneExpHt == TRUE) {
    target.g <- purrr::map_df(data_ga$ID, function(x) {
      tmp <- data_ga %>% dplyr::filter(ID == x)
      coregene <- unique(unlist(strsplit(tmp$core_enrichment, 
                                         split = "\\/")))
      output <- data.frame(gene_name = coregene, ID = x, 
                           Description = tmp$Description) %>% dplyr::distinct(., 
                                                                              gene_name, .keep_all = TRUE)
    })
    gpos <- if (kegg == TRUE) {
      match(target.g$gene_name, object@gene2Symbol)
    }
    else {
      match(target.g$gene_name, names(object@geneList))
    }
    ginfo <- target.g %>% dplyr::mutate(gpos = gpos) %>% 
      dplyr::arrange(gpos)
    if (scale.exp == TRUE) {
      gexp <- t(scale(t(exp[, 2:ncol(exp)]), scale = TRUE, 
                      center = TRUE)) %>% data.frame()
      gexp$gene_name <- exp[, 1]
    }
    else {
      gexp <- exp
      colnames(gexp)[1] <- "gene_name"
    }
    exp.long <- gexp %>% dplyr::filter(gene_name %in% unique(ginfo$gene_name)) %>% 
      dplyr::left_join(., ginfo[, 1:2], by = "gene_name") %>% 
      reshape2::melt(., id.vars = c("gene_name", "ID"))
    exp.long$gene_name <- factor(exp.long$gene_name, levels = unique(ginfo$gene_name))
    if (!is.null(sample.order)) {
      exp.long$variable <- factor(exp.long$variable, levels = sample.order)
    }
    if (!is.null(termID.order)) {
      exp.long$ID <- factor(exp.long$ID, levels = termID.order)
    }
    ght <- ggplot2::ggplot(exp.long) + ggplot2::geom_tile(ggplot2::aes(x = gene_name, 
                                                                       y = variable, fill = value), color = NA, show.legend = ht.legend) + 
      ggplot2::theme_bw(base_size = 14) + ggplot2::coord_cartesian(expand = 0) + 
      ggplot2::scale_fill_gradient2(low = exp.col[1], 
                                    mid = exp.col[2], high = exp.col[3], midpoint = 0, 
                                    name = "Z-Score") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                                                                           vjust = 0.5, hjust = 1, size = ght.geneText.size), 
                                                                       axis.text = ggplot2::element_text(color = "black"), 
                                                                       axis.ticks.x = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), 
                                                                       plot.margin = ggplot2::margin(t = -0.1, r = 0.2, 
                                                                                                     b = 0.2, l = 0.2, unit = "cm")) + ggplot2::scale_y_discrete(position = "right") + 
      ggplot2::xlab("") + ggplot2::ylab("")
    if (ght.facet == TRUE) {
      fght <- ght + ggplot2::facet_wrap(~ID, ncol = 1, 
                                        scales = ght.facet.scale, strip.position = "left") + 
        ggplot2::theme(strip.background = ggplot2::element_rect(color = NA, 
                                                                fill = "grey90"), strip.placement = "outside")
    }
    else {
      fght <- ght
    }
  }
  if (newGsea == FALSE) {
    if (subPlot == 1) {
      pres <- pLabelOut
    }
    else if (subPlot == 2) {
      if (rmHt == FALSE) {
        pres <- aplot::plot_list(gglist = list(pLabelOut, 
                                               pseg_ht), ncol = 1, heights = c(0.8, 0.2))
      }
      else {
        pres <- aplot::plot_list(gglist = list(pLabelOut, 
                                               pseg), ncol = 1, heights = c(0.8, 0.2))
      }
    }
    else if (subPlot == 3) {
      if (rmHt == FALSE) {
        pres <- aplot::plot_list(gglist = list(pLabelOut, 
                                               pseg_ht1, prank), ncol = 1, heights = c(0.5, 
                                                                                       0.2, 0.3))
      }
      else {
        pres <- aplot::plot_list(gglist = list(pLabelOut, 
                                               pseg, prank), ncol = 1, heights = c(0.5, 0.2, 
                                                                                   0.3))
      }
    }
    else {
      print("Please give 1/2/3 parameters!")
    }
  }
  else {
    if (rmHt == FALSE) {
      pres <- aplot::plot_list(gglist = list(pLabelOut, 
                                             ht), ncol = 1, heights = c(0.9, 0.1))
    }
    else {
      pres <- pLabelOut
    }
  }
  if (add.geneExpHt == TRUE) {
    pfinal <- aplot::plot_list(gglist = list(pres + ggplot2::xlab(""), 
                                             fght), ncol = 1, heights = c(1 - ght.relHight, ght.relHight))
  }
  else {
    pfinal <- pres
  }
  return(pfinal)
}
