#' @title DEG_calculate
#' 
#' @description Calculate Differentially expressed genes between 
#'     high- and low- gene expressed samples or two specific groups.
#'
#' @usage DEG_calculate(exprset = count,
#'  dtype = "count",
#'  method = "DESeq2",
#'  gene = "NEIL3")
#'
#' @param exprset The exprset to perform DEG analysis
#' @param dtype The type of exprset, either "count", "fpkm", "tpm", or "microarray"
#' @param method The method used to perform DEG analysis, either "DESeq2", "edgeR", or "limma"
#' @param gene The gene of interest. Param `gene` and param `samples_*` should not given at the same time.
#' @param percent Float, (0,1]. The percent of samples used in this analysis. e.g. percent=0.8 equals to top 40% and tail 40% samples. Default is 1.
#' @param samples_con The vector of samples treated as control.
#' @param samples_case The vector of samples treated as case.
#' 
#' @return The data.frame of DEGs
#' 
#' @export
#' 
#' @example DEG_calculate_gene(exprset=exprset, dtype="count",method="DESeq2",
#'    gene="PTK7", group_con = NULL, group_case = NULL)
#'    
#' @author GM.W
#'
DEG_calculate <- function(exprset, dtype=c("count","fpkm","tpm","microarray"),
                          method=c("DESeq2","edgeR","limma_voom","limma","wilcox"),
                          gene=NULL, percent=1,samples_con=NULL, samples_case=NULL)
{
  # 1. 判断是否传递分组参数
  if(is.null(samples_con) & is.null(samples_case) & !is.null(gene)){
    # 传入了某个基因
    rt_tmp <- data.frame(samples = colnames(exprset),t(exprset[gene,]))
    rt_tmp <- rt_tmp[order(rt_tmp[,2]),,drop=F]
    if(percent <=0 | percent > 1){
      stop("Error! The percent should be in (0,1]!")
    }
    num_samples <- round(percent * nrow(rt_tmp) / 2)
    rt_tmp <- rbind(head(rt_tmp, num_samples), tail(rt_tmp, num_samples))
    rt_tmp <- rt_tmp[!duplicated(rt_tmp$samples),]
    exprset <- exprset[,rt_tmp$samples]
    group_list <- ifelse(rt_tmp[,2] > median(rt_tmp[,2]),"case","con")
    
  }else if(is.null(gene) & !is.null(samples_con) & !is.null(samples_case)){
    # 传入了con & case 的样本
    exprset <- exprset[,c(samples_con,samples_case)]
    group_list <- c(rep("con",length(samples_con)), rep("case",length(samples_case)))
  }else{
    # 传入参数错误
      stop("Error: params error. Param gene and params samples_* should not given at the same time.")
  }
  group_list <- factor(group_list,levels = c("con","case"))
  
  # 2. 预处理
  exprset <- exprset[rowSums(exprset) > 0,]
  # 3. 判断method
  dtype <- dtype[1]
  method <- method[1]
  if(dtype == "count" & method %in% c("DESeq2","edgeR","limma_voom")){
    # count
    # anti-log, if necessary
    if(max(exprset) < 100){
      exprset <- 2^exprset - 1
      exprset <- round(exprset)
    }
    res <- switch (method,
      "DESeq2" = .DEG_DESeq2(exprset,group_list),
      "edgeR" = .DEG_edgeR(exprset,group_list),
      "limma_voom" = .DEG_limma_voom(exprset,group_list),
      stop("Unmatched method: ", method)
    )
    
  }else if(dtype %in% c("fpkm","tpm") & method %in% c("limma","wilcox")){
    # fpkm
    # log transfrom, if necessary
    if(max(exprset) > 100) exprset <- log2(exprset+1)
    res <- switch (method,
      "limma" = .DEG_limma(exprset,group_list),
      "wilcox" = .DEG_wilcox(exprset,group_list),
      stop("Unmatched method: ", method)
    )
    
  }else if(dtype == "microarray" & method =="limma"){
    # microarray
    res <- .DEG_limma(exprset,group_list)
    
  }else{
    stop("Error: params error. The dtype and method you select don't match.")
  }
  
  return(res)
}

#' @title DEG_volcano_plot
#' 
#' @description Plot a volcano graph based on DEG result using tools in c("DESeq2","edgeR","limma_voom","limma","wilcox")
#' 
#' @usage DEG_volcano_plot(DEG_res,
#'  method="DESeq2",
#'  logFC_threshold=1, 
#'  pvalue_threshold=0.05, 
#'  label_num_per=10,
#'  title=NULL)
#' 
#' @param DEG_res The result of `DEG_calculate`
#' @param method The method used to perform DEG analysis
#' @param logFC_threshold The threshold of log2FoldChange. Default is 1.
#' @param pvalue_threshold The threshold of pvalue. Default is 0.05.
#' @param label_num_per The number of gene to be shown in each side.
#' @param title The title of the plot.
#' 
#' @return A ggplot object.
#' 
#' @export
#' 
#' @examples DEG_volcano_plot(DEG_res,method=c"DESeq2")
#' 
#' @author GM.W
#' 
DEG_volcano_plot <- function(DEG_res, method=c("DESeq2","edgeR","limma_voom","limma","wilcox"),
                             logFC_threshold=1, pvalue_threshold=0.05, label_num_per=10,title=NULL){
  # 4. 分类 Up & Down
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
  
  DEG_res$change <- ifelse(DEG_res[,padj] < pvalue_threshold & abs(DEG_res[,log2fc]) > logFC_threshold,
                       ifelse(DEG_res[,log2fc] > logFC_threshold,"Up","Down"),"NoSig")
  DEG_res$change <- factor(DEG_res$change, levels = c("Up","NoSig","Down"))
  anno <- DEG_res[!is.na(DEG_res[,padj]),]
  # anno <- na.omit(anno)
  anno$weight <- anno[,log2fc] * (-log10(anno[,padj]))
  anno <- anno[order(anno$weight,decreasing = T),]
  if(!is.null(label_num_per)){
    label_gene <- rownames(anno)[c(1:label_num_per,(nrow(anno)-label_num_per+1):nrow(anno))]
    anno$label <- NA
    anno[label_gene,"label"] <- label_gene
  }
  n_regulate <- table(anno$change)
  anno_label <- paste0("Up-regulated: ",n_regulate[1],"\nNot significant: ",n_regulate[2],
                       "\nDown-regulated: ",n_regulate[3])
  anno[,"x"] <- anno[,log2fc]
  anno[,"y"] <- -log10(anno[,padj])
  
  # 作图
  require(ggplot2)
  require(ggrepel)
  p <- ggplot(anno,aes(x=x, y=y,size= y,color=change,fill=change))+
    geom_point(shape=21,alpha=0.7)+
    geom_hline(yintercept = -log10(pvalue_threshold),lty=3,color="#999999")+
    geom_vline(xintercept = c(-logFC_threshold,logFC_threshold),lty=3,color="#999999")+
    labs(x="log2FoldChange",y="-log10(FDR)")+
    annotate("text",x=0.5*range(anno[,"x"])[2],y=0,
             label=anno_label,fontface = 'italic',hjust=0)+
    scale_color_manual(values = c("#ed0000","#bebebe","#00468b"))+
    scale_fill_manual(values = c("#ed0000","#bebebe","#00468b"))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")
  if(!is.null(title)){  p <- p+ggtitle(title)}
  if(!is.null(label_num_per)){ p <- p+geom_text_repel(aes(label=label))}
  print(p)
  return(p)
}

# .DEG_DESeq2
.DEG_DESeq2 <- function(exprset, group_list){
  require(DESeq2)
  colData <- data.frame(row.names = colnames(exprset),condition=group_list)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = exprset,colData = colData,
                                design = ~condition)
  dds <- DESeq2::DESeq(dds)
  res <- as.data.frame(results(dds))
  return(res)
}

# .DEG_edgeR
.DEG_edgeR <- function(exprset, group_list){
  require(edgeR)
  dge <- DGEList(counts = exprset,group = group_list)
  dge$samples$lib.size <- colSums(dge$counts)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~0+group_list)
  rownames(design) <- colnames(dge)
  colnames(design) <- levels(group_list)
  dge <- estimateGLMCommonDisp(dge,design)
  dge <- estimateGLMTrendedDisp(dge,design)
  dge <- estimateGLMTagwiseDisp(dge,design)
  fit <- glmFit(dge,design)
  fit2 <- glmLRT(fit,contrast = c(-1,1))
  res <- as.data.frame(topTags(fit2,n=Inf))
  return(res)
}

# .DEG_limma_voom
.DEG_limma_voom <- function(exprset, group_list){
  require(limma)
  design <- model.matrix(~0+group_list)
  colnames(design) <- levels(group_list)
  rownames(design) <- colnames(exprset)
  dge <- DGEList(counts = exprset)
  dge <- calcNormFactors(dge)
  v <- voom(dge,design,normalize.method = "quantile")
  fit <- lmFit(v,design)
  contrasts <- paste(rev(levels(group_list)),collapse = "-")
  cont.matrix <- makeContrasts(contrasts = contrasts,levels = design)
  fit2 <- contrasts.fit(fit,cont.matrix)
  fit2 <- eBayes(fit2)
  res <- topTable(fit2,n=Inf)
  return(res)
}

# .DEG_limma
.DEG_limma <- function(exprset, group_list){
  require(limma)
  design <- model.matrix(~0+group_list)
  colnames(design) <- levels(group_list)
  rownames(design) <- colnames(exprset)
  fit <- lmFit(exprset, design)  
  contrasts <- paste(rev(levels(group_list)),collapse = "-")
  cont.matrix <- makeContrasts(contrasts=contrasts, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  res <- topTable(fit2, n=Inf)
  return(res)
}

# .DEG_wilcox
.DEG_wilcox <- function(exprset, group_list){
  n <- table(group_list)
  res <- data.frame()
  for (i in rownames(exprset)) {
    rt=data.frame(expression=t(exprset[i,]),group=group_list)
    colnames(rt)[1] <- "expression"
    wilcoxTest <- wilcox.test(expression~group, data = rt)
    baseMeans <- mean(rt[, 1])
    conGeneMeans <- mean(rt[1:n[1], 1])
    caseGeneMeans <- mean(rt[(n[1]+1):nrow(rt), 1])
    logFC <- log2(caseGeneMeans+0.001) - log2(conGeneMeans+0.001)
    pvalue <- wilcoxTest$p.value
    res_tmp <- cbind(gene=i,logFC=logFC,Meanbase=baseMeans,MeanCon=conGeneMeans,MeanCase=caseGeneMeans,P.Value=pvalue)
    res <- rbind(res,res_tmp)
  }
  rownames(res) <- res$gene
  res <- res[,-1]
  for (i in colnames(res)) {
    res[,i] <- as.numeric(res[,i])
    res[,i] <- round(res[,i], 6)
  }
  res[which(res$P.Value == "NaN"),"P.Value"] <- NA
  return(res)
}
