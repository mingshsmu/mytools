# rm(list = ls())


# lrpairs <- SpaTalk::lrpairs
# pathways <- SpaTalk::pathways

#' @title filter_lr_paths
#' @description find out the LR pairs with the receptors have
#' at least 1 downstream TF and pathways.
#'
#' @param st_dec_data Deconvoluted Spatial matrix.
#' @param species Character. "Human" or "Mouse". Default is "Human".
#' @param lrpairs Ligand-Receptor interaction pairs. Default is `SpaTalk::lrpairs`.
#' @param pathways Pathways with source and destination. Default is `SpaTalk::pathways`.
#' @param max_hop Numeric. The maximum hop to calculate neighboring cell-cell interaction. Default is 4.
#' @param if_doParallel Logical. Whether to use multi-cores to accelerate the process. Default is True.
#' @param use_n_cores Numeric. How many cores to use for this process. Default is (all cores - 2).
#'
#' @return A list consisting of (lrpairs, pathways, max_hop)
#' @export
#'
#' @examples lr_paths_filtered <- filter_lr_paths(st_dec_data)
#' lrpairs_filtered <- lr_paths_filtered$lrpairs
#' pathways_filtered <- lr_paths_filtered$pathways
#' max_hop_used <- lr_paths_filtered$max_hop
#'
#' @author GM.W
#'
filter_lr_paths <- function(st_dec_data, species = "Human",
                          lrpairs = NULL, pathways = NULL, max_hop = NULL,
                          if_doParallel = T,use_n_cores = NULL){
  # Parameters ----------------------------------
  if(!species %in% c("Human","Mouse")){
    stop("Please provide a correct species name! 'Human' and 'Mouse' are avaliable.")
  }
  if(is.null(lrpairs)){
    lrpairs <- SpaTalk::lrpairs
  }
  if(is.null(pathways)){
    pathways <- SpaTalk::pathways
  }
  if (!all(c("ligand", "receptor", "species") %in% names(lrpairs))) {
    stop("Please provide a correct lrpairs data.frame! See demo_lrpairs()!")
  }
  if (!all(c("src", "dest", "pathway", "source", "type", "src_tf",
             "dest_tf", "species") %in% names(pathways))) {
    stop("Please provide a correct pathways data.frame! See demo_pathways()!")
  }
  if(is.null(max_hop)){
    max_hop <- 4
  }
  if (is.null(use_n_cores)) {
    n_cores <- parallel::detectCores()
    n_cores <- n_cores - 2
  }else {
    n_cores <- use_n_cores
  }
  if (if_doParallel) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
  }

  # Data Organization ---------------------------
  # Filter the non-zero expression genes
  gene_expressed_ratio <- rowSums(st_dec_data)
  st_dec_data <- st_dec_data[which(gene_expressed_ratio > 0),]
  if (nrow(st_dec_data) == 0) {
    stop("No expressed genes in st_dec_data!")
  }
  lrpair <- lrpairs[lrpairs$species == species, ]
  lrpair <- lrpair[lrpair$ligand %in% rownames(st_dec_data) &
                     lrpair$receptor %in% rownames(st_dec_data), ]
  if (nrow(lrpair) == 0) {
    stop("No ligand-recepotor pairs found in st_dec_data!")
  }
  pathways <- pathways[pathways$species == species, ]
  pathways <- pathways[pathways$src %in% rownames(st_dec_data) &
                         pathways$dest %in% rownames(st_dec_data), ]
  ggi_tf <- pathways[, c("src", "dest", "src_tf", "dest_tf")]
  ggi_tf <- unique(ggi_tf)
  lrpair <- lrpair[lrpair$receptor %in% ggi_tf$src, ]
  if (nrow(lrpair) == 0) {
    stop("No downstream target genes found for receptors!")
  }

  # Processing ----------------------------------
  # Find the downstream activated TF and targets of receptors
  cat(crayon::blue(">>>Filter lrpairs and pathways",
                   "\n"))
  res_ggi <- NULL
  receptor_name <- unique(lrpair$receptor)
  if (if_doParallel) {
    res_ggi <- foreach::foreach(i = 1:length(receptor_name),
                                .combine = "c", .packages = "Matrix") %dopar% {
                                  ggi_res <- NULL
                                  lr_receptor <- receptor_name[i]
                                  ggi_tf1 <- ggi_tf[ggi_tf$src == lr_receptor, ]
                                  ggi_tf1 <- unique(ggi_tf1[ggi_tf1$dest %in% rownames(st_dec_data),])
                                  if (nrow(ggi_tf1) > 0) {
                                    n <- 0
                                    ggi_tf1$hop <- n + 1
                                    while (n <= max_hop) {
                                      ggi_res <- rbind(ggi_res, ggi_tf1)
                                      ggi_tf1 <- ggi_tf[ggi_tf$src %in% ggi_tf1$dest,
                                      ]
                                      ggi_tf1 <- unique(ggi_tf1[ggi_tf1$dest %in% rownames(st_dec_data), ])
                                      if (nrow(ggi_tf1) == 0) {
                                        break
                                      }
                                      ggi_tf1$hop <- n + 2
                                      n <- n + 1
                                    }
                                    ggi_res <- unique(ggi_res)
                                    ggi_res_yes <- ggi_res[ggi_res$src_tf == "YES" |
                                                             ggi_res$dest_tf == "YES", ]
                                    if (nrow(ggi_res_yes) > 0) {
                                      lr_receptor
                                    }
                                  }
                                }
  }else {
    for (i in 1:length(receptor_name)) {
      ggi_res <- NULL
      lr_receptor <- receptor_name[i]
      ggi_tf1 <- ggi_tf[ggi_tf$src == lr_receptor, ]
      ggi_tf1 <- unique(ggi_tf1[ggi_tf1$dest %in% rownames(st_dec_data),
      ])
      if (nrow(ggi_tf1) > 0) {
        n <- 0
        ggi_tf1$hop <- n + 1
        while (n <= max_hop) {
          ggi_res <- rbind(ggi_res, ggi_tf1)
          ggi_tf1 <- ggi_tf[ggi_tf$src %in% ggi_tf1$dest, ]
          ggi_tf1 <- unique(ggi_tf1[ggi_tf1$dest %in%
                                      rownames(st_dec_data), ])
          if (nrow(ggi_tf1) == 0) {break}
          ggi_tf1$hop <- n + 2
          n <- n + 1
        }
        ggi_res <- unique(ggi_res)
        ggi_res_yes <- ggi_res[ggi_res$src_tf == "YES" |
                                 ggi_res$dest_tf == "YES", ]
        if (nrow(ggi_res_yes) > 0) {
          res_ggi <- c(res_ggi, lr_receptor)
        }
      }
    }
  }
  if (length(res_ggi) == 0) {
    stop("No downstream transcriptional factors found for receptors!")
  }
  if (if_doParallel) {
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
  }
  lrpair <- lrpair[lrpair$receptor %in% res_ggi, ]
  if (nrow(lrpair) == 0) {
    stop("No ligand-recepotor pairs found!")
  }
  cat(crayon::green("***Done***", "\n"))

  # Return --------------------------------------
  return(list(lrpairs = lrpair, pathways = pathways, max_hop = max_hop))
}


#' @title find_CCI
#'
#' @description Find out the
#' @usage find_CCI(st_dec_data, st_dec_meta,
#'   lr_paths_filtered,
#'   celltype_sender, celltype_receiver,
#'   n_neighbor = 10, min_pairs = 5, min_pairs_ratio = 0,
#'   per_num = 1000, pvalue = 0.05, co_exp_ratio = 0.1,
#'   if_doParallel = TRUE, use_n_cores = NULL)
#'
#' @param st_dec_data Deconvoluted spatial matrix.
#' @param st_dec_meta The metadata whose row names are cells and contains cell, x, y, and celltype.
#' @param lr_paths_filtered The result of `filter_lr_paths`, which contains filtered lrpairs, pathways, and max_hop.
#' @param celltype_sender Character. The sender cell type.
#' @param celltype_receiver Character. The receiver cell type.
#' @param n_neighbor Numeric. Number of neighbor cells to select as the proximal cell-cell pair. Default is 10.
#' @param min_pairs Numeric. Min proximal cell-cell pairs between for sending and receiving cell types. Default is 5.
#' @param min_pairs_ratio Min proximal cell-cell pairs ratio between for sending and receiving cell types. Default is 0.
#' @param per_num Number of repeat times for permutation test. Default is 1000.
#' @param pvalue Include the significantly proximal LR pairs with this cutoff of p value from permutation test. Default is 0.05.
#' @param co_exp_ratio Min cell ratio in receiving cells with co-expressed source and target genes for predicting the downstream pathway activity.
#' @param if_doParallel Use doParallel. Default is TRUE.
#' @param use_n_cores Number of CPU cores to use. Default is all cores - 2.
#'
#'
#' @return A list consisting of (lrpair, tf, cell_pair, sender_receiver, params)
#'
#' @export
#'
#' @examples CCI_res <- find_CCI(st_dec_data = SpaTalk_dec_sc,st_dec_meta = SpaTalk_cellmeta,
#'   lr_paths_filtered = lr_paths_filtered,
#'   celltype_sender = "Fibroblasts", celltype_receiver = "Cancer_clone_A",
#'   n_neighbor = 10)
#'   lrpair <- CCI_res$lrpair
#'   TF <- CCI_res$tf
#'   cell_pair <- CCI_res$cell_pair
#'
#'   @author GM.W
#'
find_CCI <- function(st_dec_data, st_dec_meta, lr_paths_filtered,
                     celltype_sender, celltype_receiver,
                     n_neighbor = 10, min_pairs = 5, min_pairs_ratio = 0,
                     per_num = 1000, pvalue = 0.05, co_exp_ratio = 0.1,
                     if_doParallel = TRUE, use_n_cores = NULL){
  # Functions -----------------------------------
  ###############################################
  .co_exp <- function(x) {
    x_1 <- x[1:(length(x)/2)]
    x_2 <- x[(length(x)/2 + 1):length(x)]
    x_12 <- x_1 * x_2
    x_12_ratio <- length(x_12[x_12 > 0])/length(x_12)
    return(x_12_ratio)
  }

  .co_exp_p <- function(x) {
    x_real <- x[length(x)]
    x_per <- x[-length(x)]
    x_p <- length(x_per[x_per >= x_real])/length(x_per)
    return(x_p)
  }

  .co_exp_batch <- function(st_dec_data, ggi_res, cell_pair) {
    ggi_res_temp <- unique(ggi_res[, c("src", "dest")])
    cell_receiver <- unique(cell_pair$cell_receiver)
    m <- floor(nrow(ggi_res_temp)/5000)
    i <- 1
    res <- NULL
    while (i <= (m + 1)) {
      m_int <- 5000 * i
      if (m_int < nrow(ggi_res_temp)) {
        ggi_res_temp1 <- ggi_res_temp[((i - 1) * 5000 +
                                         1):(5000 * i), ]
      }
      else {
        if (m_int == nrow(ggi_res_temp)) {
          ggi_res_temp1 <- ggi_res_temp[((i - 1) * 5000 +
                                           1):(5000 * i), ]
          i <- i + 1
        }
        else {
          ggi_res_temp1 <- ggi_res_temp[((i - 1) * 5000 +
                                           1):nrow(ggi_res_temp), ]
        }
      }
      ndata_src <- st_dec_data[ggi_res_temp1$src, cell_receiver]
      ndata_dest <- st_dec_data[ggi_res_temp1$dest, cell_receiver]
      ndata_gg <- cbind(ndata_src, ndata_dest)
      ggi_res_temp1$co_ratio <- apply(ndata_gg, 1, .co_exp)
      res <- rbind(res, ggi_res_temp1)
      i <- i + 1
    }
    res$merge_key <- paste0(res$src, res$dest)
    ggi_res$merge_key <- paste0(ggi_res$src, ggi_res$dest)
    ggi_res <- merge(ggi_res, res, by = "merge_key", all.x = T,
                     sort = F)
    ggi_res <- ggi_res[, c(2:6, 9)]
    colnames(ggi_res) <- c("src", "dest", "src_tf", "dest_tf",
                           "hop", "co_ratio")
    return(ggi_res)
  }

  #' @title .generate_ggi_res
  #' @description 找到receptor中共表达的TF和target
  .generate_ggi_res <- function(ggi_tf, cell_pair, receptor_name,
                                st_dec_data, max_hop, co_exp_ratio) {
    ggi_res <- NULL
    ggi_tf1 <- ggi_tf[ggi_tf$src == receptor_name, ]
    ggi_tf1 <- unique(ggi_tf1[ggi_tf1$dest %in% rownames(st_dec_data),])
    n <- 0
    ggi_tf1$hop <- n + 1
    while (n <= max_hop) {
      ggi_res <- rbind(ggi_res, ggi_tf1)
      ggi_tf1 <- ggi_tf[ggi_tf$src %in% ggi_tf1$dest,
      ]
      ggi_tf1 <- unique(ggi_tf1[ggi_tf1$dest %in% rownames(st_dec_data),
      ])
      if (nrow(ggi_tf1) == 0) {
        break
      }
      ggi_tf1$hop <- n + 2
      n <- n + 1
    }
    ggi_res <- unique(ggi_res)
    ggi_res_temp <- unique(ggi_res[, c("src", "dest")])
    if (nrow(ggi_res_temp) >= 5000) {
      ggi_res <- .co_exp_batch(st_dec_data, ggi_res, cell_pair)
    }
    else {
      ndata_src <- st_dec_data[ggi_res$src, cell_pair$cell_receiver]
      ndata_dest <- st_dec_data[ggi_res$dest, cell_pair$cell_receiver]
      ndata_gg <- cbind(ndata_src, ndata_dest)
      ggi_res$co_ratio <- apply(ndata_gg, 1, .co_exp)
    }
    ggi_res <- ggi_res[ggi_res$co_ratio > co_exp_ratio,
    ]
    return(ggi_res)
  }

  #' @title .generate_tf_gene_all
  #' @description 找到所有的TF及其所在的hop层数
  .generate_tf_gene_all <- function(ggi_res, max_hop) {
    tf_gene_all <- NULL
    ggi_hop <- ggi_res[ggi_res$hop == 1, ]
    for (k in 1:max_hop) {
      ggi_hop_yes <- ggi_hop[ggi_hop$dest_tf == "YES",]
      if (nrow(ggi_hop_yes) > 0) {
        ggi_hop_tf <- ggi_res[ggi_res$hop == k + 1,]
        if (nrow(ggi_hop_tf) > 0) {
          ggi_hop_yes <- ggi_hop_yes[ggi_hop_yes$dest %in% ggi_hop_tf$src, ]
          if (nrow(ggi_hop_yes) > 0) {
            tf_gene <- ggi_hop_yes$hop
            names(tf_gene) <- ggi_hop_yes$dest
            tf_gene_all <- c(tf_gene_all, tf_gene)
          }
        }
      }
      ggi_hop_no <- ggi_hop[ggi_hop$dest_tf == "NO", ]
      ggi_hop <- ggi_res[ggi_res$hop == k + 1, ]
      ggi_hop <- ggi_hop[ggi_hop$src %in% ggi_hop_no$dest, ]
    }
    return(tf_gene_all)
  }
  #' @title .get_mediator
  .get_mediator <- function (ggi_res, tf_gene, tf_hop, receptor)
  {
    tf_path <- NULL
    ggi_res1 <- ggi_res[ggi_res$dest == tf_gene & ggi_res$hop ==
                          tf_hop, ]
    if (tf_hop > 1) {
      tf_hop_new <- tf_hop - 1
      for (i in tf_hop_new:1) {
        ggi_res2 <- ggi_res[ggi_res$dest %in% ggi_res1$src &
                              ggi_res$hop == i, ]
        ggi_res1 <- ggi_res1[ggi_res1$src %in% ggi_res2$dest,
        ]
        ggi_res2 <- ggi_res2[ggi_res2$dest %in% ggi_res1$src,
        ]
        if (i == tf_hop_new) {
          tf_path <- rbind(tf_path, ggi_res1, ggi_res2)
        }
        else {
          tf_path <- rbind(tf_path, ggi_res2)
        }
        ggi_res1 <- ggi_res2
      }
    }
    else {
      tf_path <- ggi_res1
    }
    tf_path_new <- NULL
    ggi_res1 <- tf_path[tf_path$src == receptor & tf_path$hop ==
                          1, ]
    if (tf_hop > 1) {
      for (i in 2:tf_hop) {
        ggi_res2 <- tf_path[tf_path$src %in% ggi_res1$dest &
                              tf_path$hop == i, ]
        # ggi_res2 <- ggi_res2[ggi_res2$src %in% ggi_res1$dest, ]
        if (i == 2) {
          tf_path_new <- rbind(tf_path_new, ggi_res1,
                               ggi_res2)
        }
        else {
          tf_path_new <- rbind(tf_path_new, ggi_res2)
        }
        ggi_res1 <- ggi_res2
      }
    }
    else {
      tf_path_new <- ggi_res1
    }

    mediators <- paste0(tf_path_new$src,"->",tf_path_new$dest,"__",tf_path_new$hop)
    mediators <- paste0(mediators, collapse = ",")
    return(mediators)
  }

  #' @title .generate_tf_res
  #' @description 找出每个TF下游激活的target数目
  .generate_tf_res <- function(tf_gene_all, celltype_sender,
                               celltype_receiver, receptor_name, ggi_res) {
    receptor_tf_temp <- data.frame(celltype_sender = celltype_sender,
                                   celltype_receiver = celltype_receiver, receptor = receptor_name,
                                   tf = names(tf_gene_all), mediator = "", n_hop = as.numeric(tf_gene_all),
                                   n_target = 0, target="", stringsAsFactors = F)
    tf_names <- names(tf_gene_all)
    tf_n_hop <- as.numeric(tf_gene_all)
    for (i in 1:length(tf_names)) {
      # target
      ggi_res_tf <- ggi_res[ggi_res$src == tf_names[i] &
                              ggi_res$hop == tf_n_hop[i] + 1, ]
      receptor_tf_temp$n_target[i] <- length(unique(ggi_res_tf$dest))
      receptor_tf_temp$target[i] <- paste0(unique(ggi_res_tf$dest),collapse = ",")
      # mediator
      receptor_tf_temp$mediator[i] <- .get_mediator(ggi_res, tf_names[i], tf_n_hop[i], receptor_name)
    }
    return(receptor_tf_temp)
  }

  #' @title .random_walk
  #' @description random walk检验，计算TF score.
  .random_walk <- function(receptor_tf, ggi_res) {
    receptor_name <- unique(receptor_tf$receptor)
    tf_score <- rep(0, nrow(receptor_tf))
    names(tf_score) <- receptor_tf$tf
    max_n_hop <- 10
    for (i in 1:10000) {
      ggi_res1 <- ggi_res[ggi_res$src == receptor_name, ]
      n_hop <- 1
      while (nrow(ggi_res1) > 0 & n_hop <= max_n_hop) {
        target_name <- sample(x = 1:nrow(ggi_res1),
                              size = 1)
        ggi_res1 <- ggi_res1[target_name, ]
        if (ggi_res1$dest %in% names(tf_score)) {
          tf_score[ggi_res1$dest] <- tf_score[ggi_res1$dest] +
            1
        }
        ggi_res1 <- ggi_res[ggi_res$src == ggi_res1$dest,
        ]
        n_hop <- n_hop + 1
      }
    }
    tf_score <- as.numeric(tf_score/10000)
    return(tf_score)
  }
  ###############################################

  #' @title st_dist
  #' @description Calculate the distance between cells
  st_dist <- function (st_dec_meta)
  {
    st_dist <- as.matrix(stats::dist(x = cbind(st_dec_meta$x, st_dec_meta$y)))
    rownames(st_dist) <- st_dec_meta[, "cell"]
    colnames(st_dist) <- st_dec_meta[, "cell"]
    return(st_dist)
  }

  #' @title get_cellpair
  #' @description Find cell-cell pairs in neighborhood.
  get_cellpair <- function (celltype_dist, st_dec_meta, celltype_sender, celltype_receiver,
                            n_neighbor)
  {
    cell_sender <- st_dec_meta[st_dec_meta$celltype == celltype_sender,
    ]
    cell_receiver <- st_dec_meta[st_dec_meta$celltype == celltype_receiver,
    ]
    cell_pair <- list()
    for (j in 1:nrow(cell_sender)) {
      celltype_dist1 <- celltype_dist[, cell_sender$cell[j]]
      celltype_dist1 <- celltype_dist1[celltype_dist1 > 0]
      celltype_dist1 <- celltype_dist1[order(celltype_dist1)]
      cell_pair[[j]] <- names(celltype_dist1[1:n_neighbor])
    }
    names(cell_pair) <- cell_sender$cell
    cell_pair <- as.data.frame(cell_pair, stringsAsFactors = F)
    cell_pair <- reshape2::melt(cell_pair, measure.vars = colnames(cell_pair),
                                variable.name = "cell_sender", value.name = "cell_receiver")
    cell_pair$cell_sender <- as.character(cell_pair$cell_sender)
    cell_pair <- cell_pair[cell_pair$cell_receiver %in% cell_receiver$cell,]
    return(cell_pair)
  }

  #' @title lr_distance
  #' @description Find coexpressed Ligand-Receptor pairs in sender-receiver cells and calculate P-value.
  lr_distance <- function (st_dec_data, cell_pair, lrdb, celltype_sender, celltype_receiver,
                           per_num, pvalue)
  {
    lrdb$celltype_sender <- celltype_sender
    lrdb$celltype_receiver <- celltype_receiver
    lrdb$lr_co_exp_num <- 0
    lrdb$lr_co_ratio <- 0
    lrdb$lr_co_ratio_pvalue <- 1
    rownames(lrdb) <- 1:nrow(lrdb)
    ndata_ligand <- st_dec_data[lrdb$ligand, cell_pair$cell_sender]
    ndata_receptor <- st_dec_data[lrdb$receptor, cell_pair$cell_receiver]
    if (nrow(lrdb) == 1) {
      ndata_lr <- ndata_ligand * ndata_receptor
      lrdb$lr_co_exp_num <- length(ndata_lr[ndata_lr > 0])
      lrdb$lr_co_ratio <- length(ndata_lr[ndata_lr > 0])/length(ndata_lr)
    }
    else {
      ndata_lr <- cbind(ndata_ligand, ndata_receptor)
      lrdb$lr_co_ratio <- apply(ndata_lr, 1, .co_exp)
      lrdb$lr_co_exp_num <- apply(ndata_lr, 1, .co_exp) *
        nrow(cell_pair)
    }
    res_per <- list()
    for (j in 1:per_num) {
      set.seed(j)
      cell_id <- sample(x = 1:ncol(st_dec_data), size = nrow(cell_pair) *
                          2, replace = T)
      ndata_ligand <- st_dec_data[lrdb$ligand, cell_id[1:(length(cell_id)/2)]]
      ndata_receptor <- st_dec_data[lrdb$receptor, cell_id[(length(cell_id)/2 +
                                                          1):length(cell_id)]]
      if (nrow(lrdb) == 1) {
        ndata_lr <- ndata_ligand * ndata_receptor
        res_per[[j]] <- length(ndata_lr[ndata_lr > 0])/length(ndata_lr)
      }
      else {
        ndata_lr <- cbind(ndata_ligand, ndata_receptor)
        res_per[[j]] <- apply(ndata_lr, 1, .co_exp)
      }
    }
    names(res_per) <- paste0("P", 1:length(res_per))
    res_per <- as.data.frame(res_per)
    res_per$real <- lrdb$lr_co_ratio
    lrdb$lr_co_ratio_pvalue <- apply(res_per, 1, .co_exp_p)
    lrdb <- lrdb[lrdb$lr_co_ratio_pvalue < pvalue, ]
    return(lrdb)
  }

  #' @title get_tf_res
  #' @description Find downstream activated TF and targets of receptors.
  get_tf_res <- function (celltype_sender, celltype_receiver, lrdb, ggi_tf,
                          cell_pair, st_dec_data, max_hop, co_exp_ratio)
  {
    receptor_tf <- NULL
    receptor_name <- unique(lrdb$receptor)
    for (j in 1:length(receptor_name)) {
      ggi_res <- .generate_ggi_res(ggi_tf, cell_pair, receptor_name[j],
                                   st_dec_data, max_hop, co_exp_ratio)
      if (nrow(ggi_res) > 0) {
        tf_gene_all <- .generate_tf_gene_all(ggi_res, max_hop)
        tf_gene_all <- data.frame(gene = names(tf_gene_all),
                                  hop = tf_gene_all, stringsAsFactors = F)
        tf_gene_all_new <- unique(tf_gene_all)
        tf_gene_all <- tf_gene_all_new$hop
        names(tf_gene_all) <- tf_gene_all_new$gene  # @@@ unique一下要废这么大劲？
        ggi_res$dest_tf_enrich <- "NO"
        if (!is.null(tf_gene_all)) {
          ggi_res[ggi_res$dest %in% names(tf_gene_all),
          ]$dest_tf_enrich <- "YES"
          receptor_tf_temp <- .generate_tf_res(tf_gene_all,
                                               celltype_sender, celltype_receiver, receptor_name[j],
                                               ggi_res)
          receptor_tf_temp$score <- .random_walk(receptor_tf_temp,
                                                 ggi_res)
          receptor_tf <- rbind(receptor_tf, receptor_tf_temp)
        }
      }
    }
    return(receptor_tf)
  }

  #' @title get_score
  #' @description Calculate the communication score between sender and receiver.
  get_score <- function (lrdb, receptor_tf)
  {
    lrdb$score <- 0
    for (j in 1:nrow(lrdb)) {
      receptor_name <- lrdb$receptor[j]
      score_lr <- 1 - lrdb$lr_co_ratio_pvalue[j]

      if (receptor_name %in% receptor_tf$receptor) {
        receptor_tf_temp <- receptor_tf[receptor_tf$receptor ==
                                          receptor_name, ]
        receptor_tf_temp$score_rt <- receptor_tf_temp$n_target *
          receptor_tf_temp$score/receptor_tf_temp$n_hop
        score_rt <- sum(receptor_tf_temp$score_rt) * (-1)
        score_rt <- 1/(1 + exp(score_rt))
        lrdb$score[j] <- sqrt(score_lr * score_rt)
        lrdb$score_lr[j] <- score_lr
        lrdb$score_rt[j] <- score_rt
      }
    }
    lrdb <- lrdb[lrdb$score > 0, ]
    if (nrow(lrdb) == 0) {
      return("NA")
    }
    else {
      return(lrdb)
    }
  }

  #' @title lr_distance_doParallel
  #'
  lr_distance_doParallel <- function (st_dec_data, cell_pair, lrdb, celltype_sender, celltype_receiver,
                                      per_num, pvalue)
  {
    require(doParallel)
    require(foreach)
    lrdb$celltype_sender <- celltype_sender
    lrdb$celltype_receiver <- celltype_receiver
    lrdb$lr_co_exp_num <- 0
    lrdb$lr_co_ratio <- 0
    lrdb$lr_co_ratio_pvalue <- 1
    rownames(lrdb) <- 1:nrow(lrdb)
    ndata_ligand <- st_dec_data[lrdb$ligand, cell_pair$cell_sender]
    ndata_receptor <- st_dec_data[lrdb$receptor, cell_pair$cell_receiver]
    if (nrow(lrdb) == 1) {
      ndata_lr <- ndata_ligand * ndata_receptor
      lrdb$lr_co_exp_num <- length(ndata_lr[ndata_lr > 0])
      lrdb$lr_co_ratio <- length(ndata_lr[ndata_lr > 0])/length(ndata_lr)
    }
    else {
      ndata_lr <- cbind(ndata_ligand, ndata_receptor)
      lrdb$lr_co_ratio <- apply(ndata_lr, 1, .co_exp)
      lrdb$lr_co_exp_num <- apply(ndata_lr, 1, .co_exp) *
        nrow(cell_pair)
    }
    res_per <- foreach::foreach(j = 1:per_num, .packages = "Matrix",
                                .export = ".co_exp") %dopar% {
                                  set.seed(j)
                                  cell_id <- sample(x = 1:ncol(st_dec_data), size = nrow(cell_pair) *
                                                      2, replace = T)
                                  ndata_ligand <- st_dec_data[lrdb$ligand, cell_id[1:(length(cell_id)/2)]]
                                  ndata_receptor <- st_dec_data[lrdb$receptor, cell_id[(length(cell_id)/2 +
                                                                                          1):length(cell_id)]]
                                  if (nrow(lrdb) == 1) {
                                    ndata_lr <- ndata_ligand * ndata_receptor
                                    length(ndata_lr[ndata_lr > 0])/length(ndata_lr)
                                  }
                                  else {
                                    ndata_lr <- cbind(ndata_ligand, ndata_receptor)
                                    as.numeric(apply(ndata_lr, 1, .co_exp))
                                  }
                                }
    names(res_per) <- paste0("P", 1:length(res_per))
    res_per <- as.data.frame(res_per)
    res_per$real <- lrdb$lr_co_ratio
    lrdb$lr_co_ratio_pvalue <- apply(res_per, 1, .co_exp_p)
    lrdb <- lrdb[lrdb$lr_co_ratio_pvalue < pvalue, ]
    return(lrdb)
  }

  #' @title get_tf_res_doParallel
  #' @description
  get_tf_res_doParallel <- function (celltype_sender, celltype_receiver, lrdb, ggi_tf,
            cell_pair, st_dec_data, max_hop, co_exp_ratio)
  {
    receptor_name <- unique(lrdb$receptor)
    receptor_tf <- NULL
    receptor_tf <- foreach::foreach(j = 1:length(receptor_name),
                                    .packages = "Matrix", .combine = "rbind",
                                    .export = c(".generate_ggi_res", ".generate_tf_gene_all",
                                                ".generate_tf_res", ".random_walk")) %dopar%
      {
        ggi_res <- .generate_ggi_res(ggi_tf, cell_pair,
                                     receptor_name[j], st_dec_data, max_hop, co_exp_ratio)
        if (nrow(ggi_res) > 0) {
          tf_gene_all <- .generate_tf_gene_all(ggi_res,
                                               max_hop)
          tf_gene_all <- data.frame(gene = names(tf_gene_all),
                                    hop = tf_gene_all, stringsAsFactors = F)
          tf_gene_all_new <- unique(tf_gene_all)
          tf_gene_all <- tf_gene_all_new$hop
          names(tf_gene_all) <- tf_gene_all_new$gene
          ggi_res$dest_tf_enrich <- "NO"
          if (!is.null(tf_gene_all)) {
            ggi_res[ggi_res$dest %in% names(tf_gene_all),
            ]$dest_tf_enrich <- "YES"
            receptor_tf_temp <- .generate_tf_res(tf_gene_all,
                                                 celltype_sender, celltype_receiver, receptor_name[j],
                                                 ggi_res)
            receptor_tf_temp$score <- .random_walk(receptor_tf_temp,
                                                   ggi_res)
            receptor_tf_temp
          }
        }
      }
    return(receptor_tf)
  }


  # Parameters ----------------------------------
  # st_dec_data
  st_dec_data <- st_dec_data[,rownames(st_dec_meta)]
  gene_expressed_ratio <- rowSums(st_dec_data)
  st_dec_data <- st_dec_data[which(gene_expressed_ratio > 0),]
  if (nrow(st_dec_data) == 0) {
    stop("No expressed genes in st_dec_data!")
  }

  # parallel
  if (is.null(use_n_cores)) {
    n_cores <- parallel::detectCores()
    n_cores <- n_cores - 2
  }else {
    n_cores <- use_n_cores
  }
  if (if_doParallel) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
  }

  # sender & receiver
  if (!celltype_sender %in% st_dec_meta$celltype) {
    stop("Please provide a correct celltype_sender")
  }
  if (!celltype_receiver %in% st_dec_meta$celltype) {
    stop("Please provide a correct celltype_receiver")
  }

  # Data organization ---------------------------
  celltype_dist <- st_dist(st_dec_meta)
  lrdb <- lr_paths_filtered$lrpairs
  pathways <- lr_paths_filtered$pathways
  max_hop <- lr_paths_filtered$max_hop

  # Processing ----------------------------------
  cell_pair <- get_cellpair(celltype_dist, st_dec_meta, celltype_sender,
                            celltype_receiver, n_neighbor)
  cell_sender <- st_dec_meta[st_dec_meta$celltype == celltype_sender,]
  cell_receiver <- st_dec_meta[st_dec_meta$celltype == celltype_receiver,]
  cell_pair_all <- nrow(cell_sender) * nrow(cell_receiver)/2

  if (nrow(cell_pair) <= min_pairs) {
    stop(paste0("Cell pairs found between ", celltype_sender,
                " and ", celltype_receiver, " less than min_pairs!"))
  }
  if (nrow(cell_pair) <= cell_pair_all * min_pairs_ratio) {
    stop(paste0("Cell pairs found between ", celltype_sender,
                " and ", celltype_receiver, " less than min_pairs_ratio!"))
  }

  ggi_tf <- unique(pathways[, c("src", "dest", "src_tf", "dest_tf")])
  cat(crayon::cyan(">>>Begin to find LR pairs", "\n"))
  if (if_doParallel) {
    lrdb <- lr_distance_doParallel(st_dec_data, cell_pair,
                                    lrdb, celltype_sender, celltype_receiver, per_num,
                                    pvalue)
  }else {
    lrdb <- lr_distance(st_dec_data, cell_pair, lrdb, celltype_sender,
                        celltype_receiver, per_num, pvalue)
  }

  # The downstream of receptors
  cat(crayon::cyan(">>>Begin to score downstream TFs and targets of receptors", "\n"))
  receptor_tf <- NULL
  if (nrow(lrdb) > 0) {
    if (if_doParallel) {
      receptor_tf <- get_tf_res_doParallel(celltype_sender,
                                            celltype_receiver, lrdb, ggi_tf, cell_pair,
                                            st_dec_data, max_hop, co_exp_ratio)
    }else {
      receptor_tf <- get_tf_res(celltype_sender, celltype_receiver,
                                lrdb, ggi_tf, cell_pair, st_dec_data, max_hop, co_exp_ratio)
    }
    if (is.null(receptor_tf)) {
      stop(paste0("No LR pairs found between ", celltype_sender,
                  " and ", celltype_receiver))
    }

    cat(crayon::cyan(">>>Begin to score inter- and intra-cellular interactions", "\n"))
    lrdb <- get_score(lrdb, receptor_tf)
  }else {
    stop(paste0("No LR pairs found between ", celltype_sender,
                " and ", celltype_receiver))
  }

  # stop parallel
  if (if_doParallel) {
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
  }

  cat(crayon::green("***Done***", "\n"))
  # Return --------------------------------------
  # 这里没有返回st_dec_data和st_dec_meta，以后会整合一下
  cellpairs <- list(cell_pair)
  names(cellpairs) <- paste0(celltype_sender, " -> ", celltype_receiver)
  res <- list(lrpair = lrdb,
              tf = receptor_tf,
              cell_pair = cellpairs,
              params = list(min_pairs = min_pairs, min_pairs_ratio = min_pairs_ratio,
                            per_num = per_num, pvalue = pvalue, co_exp_ratio = co_exp_ratio)
  )
  return(res)
}


# plot ------------------------------------------

#' @title plot_CCI_ccdist
#' @description Plot sender and receiver cells in a scatter plot using ggplot2.
#'
#' @usage plot_CCI_ccdist(st_dec_meta,
#' CCI_res,
#' celltype_sender,
#' celltype_receiver,
#' title = NULL,
#' color = NULL,
#' size = 1,
#' if_plot_others = TRUE,
#' if_plot_density = TRUE,
#' if_plot_edge = TRUE,
#' if_show_arrow = TRUE,
#' if_order_cell = TRUE,
#' legend.position = "bottom",
#' arrow_width = 0.02,
#' arrow_col = "#000000",
#' plot_cells = NULL)
#'
#' @param st_dec_meta The metadata whose row names are cells and contains cell, x, y, and celltype.
#' @param CCI_res The result of `find_CCI`, which contains data.frame including lrpair, tf, and cell_pair.
#' @param celltype_sender Character. The type of sender cells.
#' @param celltype_receiver Character. The type of receiver cells. Default is `NULL`.
#' @param title Character. The title of the plot. Default is `NULL`.
#' @param color Character. The color to be used. Default is `NULL`.
#' @param size Numeric. The size of the points. Default is 1.
#' @param if_plot_others Logical. Whether to plot other cell types? Default is `TRUE`.
#' @param if_plot_density Logical. Whether to plot the density of sender and receiver cells? Default is `TRUE`.
#' @param if_plot_edge Logical. Whether to plot segment between cell pairs? Default is `TRUE`.
#' @param if_show_arrow Logical. Whether to plot arrow from sender to receiver cells? Default is `TRUE`.
#' @param if_order_cell Logical. Whether to order cell before plotting? Default is `TRUE`.
#' @param legend.position Character or two-dimension vector. The position of legend. Default is "bottom".
#' @param arrow_width Numeric. The width of arrow. Default is 0.02.
#' @param arrow_col Character. The color of arrow. Default is "#000000" (black).
#' @param plot_cells Character. The specified cells to plot. Default is `NULL`.
#'
#' @return A ggplot object (when if_plot_density = FALSE) or ggExtraPlot object (when if_plot_density = TRUE).
#'
#' @export
#'
#' @example plot_CCI_ccdist(st_dec_meta = SpaTalk_cellmeta, CCI_res = CCI_res,
#' celltype_sender = "Fibroblasts", celltype_receiver = "Cancer_clone_A",
#' title = "SpaTalk CCI",color = c("#192a56","#c23616","grey80"))
#'
#' @author GM.W
#'
plot_CCI_ccdist <- function(st_dec_meta, CCI_res, celltype_sender, celltype_receiver, title = NULL,
                            color = NULL, size = 1, if_plot_others = TRUE, if_plot_density = TRUE,
                            if_plot_edge = TRUE, if_show_arrow = TRUE, if_order_cell = TRUE,
                            legend.position = "bottom", arrow_width = 0.02, arrow_col = "#000000",
                            plot_cells = NULL){
  # Check Input ---------------------------------
  if (!celltype_sender %in% st_dec_meta$celltype) {
    stop("Please provide the correct name of celltype_sender!")
  }
  if (!celltype_receiver %in% st_dec_meta$celltype) {
    stop("Please provide the correct name of celltype_receiver!")
  }
  if (!is.null(plot_cells)) {
    if (all(plot_cells %in% st_dec_meta$cell)) {
      st_dec_meta <- st_dec_meta[st_dec_meta$cell %in% plot_cells,
      ]
    }
    else {
      warning("Exclude the cells that is not in st_dec_meta!")
      plot_cells <- plot_cells[plot_cells %in% st_dec_meta$cell]
      if (length(plot_cells) < 10) {
        stop("Number of cells is less than 10!")
      }
      st_dec_meta <- st_dec_meta[st_dec_meta$cell %in% plot_cells,
      ]
    }
  }

  # Data Organization ---------------------------
  cell_pair <- CCI_res$cell_pair[[paste0(celltype_sender, " -> ", celltype_receiver)]]
  if (is.null(cell_pair)) {
    stop("No LR pairs found from the celltype_sender to celltype_receiver!")
  }
  # The coords of arrow
  cell_pair <- cell_pair[cell_pair$cell_sender %in% st_dec_meta$cell &
                           cell_pair$cell_receiver %in% st_dec_meta$cell, ]
  cell_pair$x1 <- 0
  cell_pair$y1 <- 0
  cell_pair$x2 <- 0
  cell_pair$y2 <- 0
  for (i in 1:nrow(cell_pair)) {
    d1 <- cell_pair$cell_sender[i]
    d2 <- st_dec_meta[st_dec_meta$cell == d1, ]
    cell_pair$x1[i] <- d2$x
    cell_pair$y1[i] <- d2$y
    d1 <- cell_pair$cell_receiver[i]
    d2 <- st_dec_meta[st_dec_meta$cell == d1, ]
    cell_pair$x2[i] <- d2$x
    cell_pair$y2[i] <- d2$y
  }
  # Rename & Classification
  st_dec_meta[!st_dec_meta$celltype %in% c(celltype_sender, celltype_receiver), ]$celltype <- "Others"
  st_dec_meta[st_dec_meta$celltype == celltype_sender, ]$celltype <- paste0("Sender: ", celltype_sender)
  st_dec_meta[st_dec_meta$celltype == celltype_receiver, ]$celltype <- paste0("Receiver: ", celltype_receiver)
  # factor celltype
  if (if_plot_others) {
    st_dec_meta$celltype <- factor(st_dec_meta$celltype,
                                   levels = c(paste0("Sender: ", celltype_sender),
                                              paste0("Receiver: ", celltype_receiver),
                                              "Others"))
  }
  else {
    st_dec_meta <- st_dec_meta[st_dec_meta$celltype != "Others", ]
    st_dec_meta$celltype <- factor(st_dec_meta$celltype,
                                   levels = c(paste0("Sender: ", celltype_sender),
                                              paste0("Receiver: ", celltype_receiver)))
  }

  # color
  if (is.null(color)) {
    cellname <- unique(as.character(st_dec_meta$celltype))
    cellname <- cellname[order(cellname)]
    if ("unsure" %in% cellname) {
      cellname <- cellname[-which(cellname == "unsure")]
    }
    col_manual <- ggpubr::get_palette(palette = "Set1",
                                      k = length(cellname)-1)
    if (if_plot_others) {
      color <- c(col_manual, "grey80")
    }

  }

  # order cells before plot
  if(if_order_cell){
    st_dec_meta <- st_dec_meta[order(st_dec_meta$celltype,decreasing = T),]
  }

  # ggplot --------------------------------------
  require(ggplot2)
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = st_dec_meta,
                        ggplot2::aes(x, y, color = celltype), size = size) +
    ggplot2::scale_color_manual(values = color) +
    ggplot2::scale_fill_manual(values = color) +
    ggplot2::theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                   panel.background = element_rect(colour = "white", fill="white"),
                   plot.background = element_rect(colour = "white", fill="white"),
                   panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
                   axis.text =element_blank(),
                   axis.ticks =element_blank(),
                   axis.title =element_blank(),
                   legend.position = legend.position,
                   legend.title=element_text(size = 10,face="bold"),
                   legend.text=element_text(size = 8),
                   legend.key = element_rect(colour = "transparent", fill = "white"),
                   legend.key.size = unit(0.45, 'cm'),
                   strip.text = element_text(size = 10,face="bold"),
                   plot.title = element_text(hjust = 0.5,face = "bold"))+
    ggplot2::guides(color=guide_legend(title="Cell Type"),
                    fill=guide_legend(title="Cell Type"))+
    ggplot2::ggtitle(title)
  if (if_plot_edge) {
    if (if_show_arrow) {
      p <- p +
        ggplot2::geom_segment(data = cell_pair,
                            aes(x = x1, y = y1, xend = x2, yend = y2),
                            arrow = arrow(length = unit(arrow_width,"inches"), type = "closed"),
                            color = arrow_col)
    }
    else {
      p <- p +
        ggplot2::geom_segment(data = cell_pair,
                            aes(x = x1, y = y1, xend = x2, yend = y2),
                            color = arrow_col)
    }
  }
  if (if_plot_density) {
    return(ggExtra::ggMarginal(p = p, type = "density", groupFill = T))
  }
  else {
    return(p)
  }
}


#' @title plot_CCI_lrpair
#'
#' @description Plot ligand of sender cells and receptor of receiver cells in a scatter plot using ggplot2.
#'
#' @usage plot_CCI_lrpair(st_dec_data, st_dec_meta, CCI_res,
#' celltype_sender, celltype_receiver, ligand, receptor,
#' title, color = NULL, size = 1, if_plot_others = TRUE,
#' if_plot_density = TRUE, if_plot_edge = TRUE,
#' if_show_arrow = TRUE, if_order_cell = TRUE, legend.position = "bottom",
#' arrow_width = 0.02, arrow_col = "#000000", plot_cells = NULL)
#'
#' @param st_dec_data Deconvoluted Spatial matrix.
#' @param st_dec_meta The metadata whose row names are cells and contains cell, x, y, and celltype.
#' @param CCI_res The result of `find_CCI`, which contains data.frame including lrpair, tf, and cell_pair.
#' @param celltype_sender Character. The type of sender cells.
#' @param ligand Character. The expressed ligand of celltype_sender cells.
#' @param celltype_receiver Character. The type of receiver cells.
#' @param receptor Character. The expressed receptor of celltype_receiver cells.
#' @param title Character. The title of the plot. Default is `NULL`.
#' @param color Character. The color to be used. Default is `NULL`.
#' @param size Numeric. The size of the points. Default is 1.
#' @param if_plot_others Logical. Whether to plot other cell types? Default is `TRUE`.
#' @param if_plot_density Logical. Whether to plot the density of sender and receiver cells? Default is `TRUE`.
#' @param if_plot_edge Logical. Whether to plot segment between cell pairs? Default is `TRUE`.
#' @param if_show_arrow Logical. Whether to plot arrow from sender to receiver cells? Default is `TRUE`.
#' @param if_order_cell Logical. Whether to order cell before plotting? Default is `TRUE`.
#' @param legend.position Character or two-dimension vector. The position of legend. Default is "bottom".
#' @param arrow_width Numeric. The width of arrow. Default is 0.02.
#' @param arrow_col Character. The color of arrow. Default is "#000000" (black).
#' @param plot_cells Character. The specified cells to plot. Default is `NULL`.
#'
#' @return A ggplot object (when if_plot_density = FALSE) or ggExtraPlot object (when if_plot_density = TRUE).
#'
#' @export
#'
#' @example plot_CCI_lrpair(st_dec_data = SpaTalk_dec_sc,st_dec_meta = SpaTalk_cellmeta,CCI_res = CCI_res,
#' celltype_sender = "Fibroblasts", celltype_receiver = "Cancer_clone_A",
#' ligand = "DCN",receptor = "EGFR",
#' title = "SpaTalk CCI",color = c("#192a56","#c23616","grey80"))
#'
#' @author GM.W
#'
plot_CCI_lrpair <- function(st_dec_data, st_dec_meta, CCI_res,
                            celltype_sender, celltype_receiver, ligand, receptor,
                            title, color = NULL, size = 1, if_plot_others = TRUE,
                            if_plot_density = TRUE, if_plot_edge = TRUE,
                            if_show_arrow = TRUE, if_order_cell = TRUE, legend.position = "bottom",
                            arrow_width = 0.02, arrow_col = "#000000", plot_cells = NULL){
  # Check Input ---------------------------------
  if (!celltype_sender %in% st_dec_meta$celltype) {
    stop("Please provide the correct name of celltype_sender!")
  }
  if (!celltype_receiver %in% st_dec_meta$celltype) {
    stop("Please provide the correct name of celltype_receiver!")
  }
  if (!ligand %in% rownames(st_dec_data)) {
    stop("Please provide the correct name of ligand!")
  }
  if (!receptor %in% rownames(st_dec_data)) {
    stop("Please provide the correct name of receptor!")
  }
  if (!is.null(plot_cells)) {
    if (all(plot_cells %in% st_dec_meta$cell)) {
      st_dec_meta <- st_dec_meta[st_dec_meta$cell %in% plot_cells,
      ]
    }
    else {
      warning("Exclude the cells that is not in st_dec_meta!")
      plot_cells <- plot_cells[plot_cells %in% st_dec_meta$cell]
      if (length(plot_cells) < 10) {
        stop("Number of cells is less than 10!")
      }
      st_dec_meta <- st_dec_meta[st_dec_meta$cell %in% plot_cells,
      ]
    }
  }

  # Data Organization ---------------------------
  cell_pair <- CCI_res$cell_pair[[paste0(celltype_sender, " -> ", celltype_receiver)]]
  if (is.null(cell_pair)) {
    stop("No LR pairs found from the celltype_sender to celltype_receiver!")
  }
  # The coords of arrow
  cell_pair <- cell_pair[cell_pair$cell_sender %in% st_dec_meta$cell &
                           cell_pair$cell_receiver %in% st_dec_meta$cell, ]
  cell_pair$x1 <- 0
  cell_pair$y1 <- 0
  cell_pair$x2 <- 0
  cell_pair$y2 <- 0
  for (i in 1:nrow(cell_pair)) {
    d1 <- cell_pair$cell_sender[i]
    d2 <- st_dec_meta[st_dec_meta$cell == d1, ]
    cell_pair$x1[i] <- d2$x
    cell_pair$y1[i] <- d2$y
    d1 <- cell_pair$cell_receiver[i]
    d2 <- st_dec_meta[st_dec_meta$cell == d1, ]
    cell_pair$x2[i] <- d2$x
    cell_pair$y2[i] <- d2$y
  }

  # Rename & Classification
  st_dec_data <- st_dec_data[, st_dec_meta$cell]
  st_dec_meta$ligand <- as.numeric(st_dec_data[ligand, ])
  st_dec_meta$receptor <- as.numeric(st_dec_data[receptor, ])
  st_dec_meta$Expressed_genes <- "Others"
  st_dec_meta_ligand <- st_dec_meta[st_dec_meta$celltype == celltype_sender &
                              st_dec_meta$ligand > 0, ]
  st_dec_meta_receptor <- st_dec_meta[st_dec_meta$celltype == celltype_receiver &
                                st_dec_meta$receptor > 0, ]
  celltype_sender <- paste0(celltype_sender, ": ", ligand)
  celltype_receiver <- paste0(celltype_receiver, ": ", receptor)
  st_dec_meta[st_dec_meta$cell %in% st_dec_meta_ligand$cell, ]$Expressed_genes <- celltype_sender
  st_dec_meta[st_dec_meta$cell %in% st_dec_meta_receptor$cell, ]$Expressed_genes <- celltype_receiver
  cell_pair <- cell_pair[cell_pair$cell_sender %in% st_dec_meta_ligand$cell, ]
  cell_pair <- cell_pair[cell_pair$cell_receiver %in% st_dec_meta_receptor$cell, ]
  st_dec_meta$Expressed_genes <- factor(st_dec_meta$Expressed_genes,
                                    levels = c(celltype_sender, celltype_receiver, "Others"))

  # color
  if (is.null(color)) {
    cellname <- unique(st_dec_meta$Expressed_genes)
    cellname <- cellname[order(cellname)]
    col_manual <- ggpubr::get_palette(palette = "Set1",
                                      k = length(cellname)-1)
    if (if_plot_others) {
      color <- c(col_manual, "grey80")
    }
  }

  # factor celltype
  if (!if_plot_others) {
    st_dec_meta <- st_dec_meta[st_dec_meta$Expressed_genes %in% c(celltype_sender, celltype_receiver),]
  }

  # order cells before plot
  if(if_order_cell){
    st_dec_meta <- st_dec_meta[order(st_dec_meta$Expressed_genes,decreasing = T),]
  }

  # ggplot --------------------------------------
  require(ggplot2)
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = st_dec_meta,
                        ggplot2::aes(x=x, y=y, color = Expressed_genes), size = size) +
    ggplot2::scale_color_manual(values = color) +
    ggplot2::scale_fill_manual(values = color) +
    ggplot2::theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                   panel.background = element_rect(colour = "white", fill="white"),
                   plot.background = element_rect(colour = "white", fill="white"),
                   panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
                   axis.text =element_blank(),
                   axis.ticks =element_blank(),
                   axis.title =element_blank(),
                   legend.position = legend.position,
                   legend.title=element_text(size = 10,face="bold"),
                   legend.text=element_text(size = 8),
                   legend.key = element_rect(colour = "transparent", fill = "white"),
                   legend.key.size = unit(0.45, 'cm'),
                   strip.text = element_text(size = 10,face="bold"),
                   plot.title = element_text(hjust = 0.5,face = "bold"))+
    ggplot2::guides(color=guide_legend(title="Cell Type"),
                    fill=guide_legend(title="Cell Type"))+
    ggplot2::ggtitle(title)
  if (if_plot_edge) {
    if (if_show_arrow) {
      p <- p +
        ggplot2::geom_segment(data = cell_pair,
                              aes(x = x1, y = y1, xend = x2, yend = y2),
                              arrow = arrow(length = unit(arrow_width,"inches"), type = "closed"),
                              color = arrow_col)
    }
    else {
      p <- p +
        ggplot2::geom_segment(data = cell_pair,
                              aes(x = x1, y = y1, xend = x2, yend = y2),
                              color = arrow_col)
    }
  }
  if (if_plot_density) {
    return(ggExtra::ggMarginal(p = p, type = "density", groupFill = T))
  }
  else {
    return(p)
  }

}


#' @title plot_CCI_lrpath
#'
#' @description Plot a network of Ligand, Receptor, Mediators, TFs, and Targets.
#'
#' @param CCI_res The result of `find_CCI`, which contains data.frames including lrpair, tf, and cell_pair.
#' @param lr_paths_filtered The result of `filter_lr_paths`, which contains data.frames including lrpairs, and pathways.
#' @param celltype_sender Character. The type of sender cells.
#' @param ligand Character. The expressed ligand of celltype_sender cells.
#' @param celltype_receiver Character. The type of receiver cells.
#' @param receptor Character. The expressed receptor of celltype_receiver cells.
#' @param TFs Character. The expressed TF genes of celltype_receiver cells. Default is `NULL`.
#' @param title Character. The title of the plot. Default is `NULL`.
#' @param color_celltype Character. The color used for celltype (sender and receiver).
#' @param color_function Character. The color used for gene function (Ligand, Receptor, Mediator, TF, Target).
#' @param size Numeric. The size of point. Default is 5.
#' @param edge_width Numeric. The width of arrow. Default is 0.1.
#' @param legend.position Character or two-demension vector. The position of legend. Default is "right".
#' @param max_hop Numeric. The maximum hop to plot. Default is `NULL`.
#' @param flip Logical. Whether to flip the plot. Default is `FALSE`.
#' @param ReturnData Logical. Whether to return the lr_path data (for further visualization in cytoscape). Default is FALSE.
#'
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @example plot_CCI_lrpath(lr_paths_filtered = lr_paths_filtered, CCI_res = CCI_res,
#'                celltype_sender = "Fibroblasts", celltype_receiver = "Cancer_clone_A",
#'                ligand = "DCN", receptor = "EGFR",title = "DCN-EGFR",max_hop = 2)
#'
#' @author GM.W
#'
plot_CCI_lrpath <- function(lr_paths_filtered,CCI_res, celltype_sender, celltype_receiver,
                            ligand, receptor, TFs = NULL, title = NULL, color_celltype = NULL,
                            color_function = NULL, size=5, edge_width = 0.1,
                            legend.position = "right", max_hop=NULL, flip=FALSE,
                            ReturnData = FALSE){
  # Check Input ---------------------------
  lrpair <- CCI_res$lrpair
  TF <- CCI_res$tf
  pathways <- lr_paths_filtered$pathways
  max_hop_used <- lr_paths_filtered$max_hop
  if (!celltype_sender %in% lrpair$celltype_sender) {
    stop("Please provide the correct name of celltype_sender!")
  }
  if (!celltype_receiver %in% lrpair$celltype_receiver) {
    stop("Please provide the correct name of celltype_receiver!")
  }
  if (!ligand %in% lrpair$ligand) {
    stop("Please provide the correct name of ligand!")
  }
  if (!receptor %in% lrpair$receptor) {
    stop("Please provide the correct name of receptor!")
  }
  check <- lrpair[lrpair$ligand == ligand & lrpair$receptor == receptor,]
  if(nrow(check) == 0){
    stop("Please provide the correct ligand-receptor pair!")
  }
  if(!is.null(TFs)){
    if(!all(TFs %in% TF$tf)){
      stop("Please provide the correct TFs!")
    }
  }
  if(!is.null(max_hop)){
    if(max_hop > max_hop_used){
      stop("Please provide the max_hop smaller than the max_hop used in `filter_lr_paths()`!")
    }
  }
  # Data Organization ---------------------------
  ### Ligand-Receptor-TF-Target path
  TF <- TF[which(TF$receptor == receptor),]
  if(!is.null(max_hop)){
    TF <- TF[which(TF$n_hop <= max_hop),]
  }
  if(!is.null(TFs)){
    TF <- TF[TF$tf %in% TFs,]
  }
  # L_R
  L_R <- data.frame(src = ligand, dest = receptor, hop = 0)
  # R_TF
  R_TF_tmp <- TF$mediator
  R_TF_tmp <- as.data.frame(stringr::str_split(R_TF_tmp, ",", simplify = T))
  R_TF <- NULL
  for (i in 1:nrow(R_TF_tmp)) {
    R_TF_1 <- as.character(R_TF_tmp[i,])
    R_TF_1 <- R_TF_1[which(R_TF_1 != "")]
    R_TF <- c(R_TF, R_TF_1)
  }
  R_TF <- as.data.frame(stringr::str_split(R_TF, "->|__",simplify = T))
  colnames(R_TF) <- c("src","dest","hop")
  R_TF <- unique(R_TF)
  # TF_Target
  TF_Target_tmp <- TF[, c("tf","n_hop","target")]
  TF_Target_genes <- as.data.frame(stringr::str_split(TF_Target_tmp$target, ",", simplify = T))
  TF_Target <- data.frame()
  for (i in 1:nrow(TF_Target_tmp)) {
    gene1 <- as.character(TF_Target_genes[i,])
    gene1 <- gene1[which(gene1 != "")]
    TF_Target_1 <- data.frame(src = TF_Target_tmp[i,"tf"], dest = gene1, hop = TF_Target_tmp[i,"n_hop"]+1)
    TF_Target <- rbind(TF_Target, TF_Target_1)
  }
  TF_Target <- unique(TF_Target)

  lr_path <- rbind(L_R, R_TF, TF_Target)
  lr_path$hop <- as.numeric(lr_path$hop)
  lr_path$TF <- "No"
  TF_gene <- unique(pathways[which(pathways$dest_tf == "YES"),"dest"])
  lr_path[lr_path$dest %in% TF_gene, "TF"] <- "Yes"
  lr_path$target <- "No"
  lr_path[lr_path$dest %in% TF_Target$dest, "target"] <- "Yes"
  # plot_node
  lr_path$hop <- lr_path$hop + 2
  plot_node <- data.frame(src = c("",lr_path$src),
                          dest = c(ligand, lr_path$dest),
                          x = c(1,lr_path$hop),
                          TF = c("No", lr_path$TF),
                          target = c("No", lr_path$target),
                          stringsAsFactors = F)
  plot_node <- plot_node[!duplicated(paste0(plot_node$dest, "_",plot_node$x)),]
  freq_y <- as.data.frame(table(plot_node$x), stringAsFactors = F)
  freq_y_max <- max(freq_y$Freq)
  plot_node$y <- 0
  # plot_node_new
  plot_node_new <- NULL
  for (i in 1:max(plot_node$x)) {
    plot_node_tmp <- plot_node[plot_node$x == i, ]
    if(nrow(plot_node_tmp) == freq_y_max){
      plot_node_tmp$y <- 1:nrow(plot_node_tmp)
    }else{
      interval_y <- (freq_y_max - 1) / (nrow(plot_node_tmp) + 1)
      y_new <- 1 + interval_y
      for (j in 1:nrow(plot_node_tmp)) {
        plot_node_tmp$y[j] <- y_new
        y_new <- y_new + interval_y
      }
    }
    plot_node_new <- rbind(plot_node_new, plot_node_tmp)

  }

  plot_node_new$celltype <- celltype_receiver
  plot_node_new[which(plot_node_new$dest == ligand),"celltype"] <- celltype_sender
  plot_node_new$celltype <- factor(plot_node_new$celltype,
                                   levels = c(celltype_sender,celltype_receiver))
  # Gene Classification
  plot_node_new$Function <- "Mediator"
  plot_node_new[which(plot_node_new$dest == ligand), "Function"] <- "Ligand"
  plot_node_new[which(plot_node_new$dest == receptor), "Function"] <- "Receptor"
  plot_node_new[which(plot_node_new$target == "Yes"), "Function"] <- "Target"
  plot_node_new[which(plot_node_new$TF == "Yes"), "Function"] <- "TF"

  plot_node_new$Function <- factor(plot_node_new$Function,
                                   levels = c("Ligand","Receptor","Mediator","TF","Target"))
  # coords of segment
  plot_segment <- plot_node_new
  plot_segment$dest_x <- plot_segment$x
  plot_segment$dest_y <- plot_segment$y
  plot_segment_tmp <- plot_segment[which(plot_segment$x > 1),]
  plot_segment_new <- data.frame()
  for (i in 1:nrow(plot_segment_tmp)) {
    plot_segment_1 <- plot_segment_tmp[i,]
    plot_segment_src <- plot_segment[plot_segment$x == plot_segment_1$x-1,]
    src_index <- which(plot_segment_1$src == plot_segment_src$dest)
    plot_segment_1$src_x <- plot_segment_src[src_index,"x"]
    plot_segment_1$src_y <- plot_segment_src[src_index,"y"]
    plot_segment_new <- rbind(plot_segment_new, plot_segment_1)

  }
  # color
  if(is.null(color_celltype)){
    color_celltype <- ggpubr::get_palette(palette = "Set3", k = 2)
  }
  if(is.null(color_function)){
    color_function <- ggpubr::get_palette(palette = "Set1", k = length(unique(plot_node_new$Function)))
  }
  # flip
  if(flip == TRUE){
    plot_node_new$x <- max(plot_node_new$x) - plot_node_new$x
    plot_segment_new$src_x <- max(plot_node_new$x) - plot_segment_new$src_x + 1
    plot_segment_new$dest_x <- max(plot_node_new$x) - plot_segment_new$dest_x + 1
  }
  # Plot ----------------------------------------
  require(ggplot2)
  p <- ggplot2::ggplot()+
    ggplot2::geom_segment(data = plot_segment_new,
                          ggplot2::aes(x=src_x,y=src_y,xend=dest_x,yend=dest_y),lwd=edge_width)+
    geom_point(data = plot_node_new,
               ggplot2::aes(x=x, y=y, color=celltype),size=size)+
    ggrepel::geom_label_repel(data = plot_node_new,
                              ggplot2::aes(x=x,y=y,label=dest,fill=Function))+
    ggplot2::scale_color_manual(values = color_celltype)+
    ggplot2::scale_fill_manual(values = color_function)+
    ggplot2::theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                   panel.background = element_rect(colour = "white", fill="white"),
                   plot.background = element_rect(colour = "white", fill="white"),
                   panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
                   axis.text =element_blank(),
                   axis.ticks =element_blank(),
                   axis.title =element_blank(),
                   legend.position = legend.position,
                   legend.title=element_text(size = 10,face="bold"),
                   legend.text=element_text(size = 8),
                   legend.key = element_rect(colour = "transparent", fill = "white"),
                   legend.key.size = unit(0.45, 'cm'),
                   strip.text = element_text(size = 10,face="bold"),
                   plot.title = element_text(hjust = 0.5,face = "bold"))+
    ggplot2::guides(color=guide_legend(title="Cell Type"),
                    fill=guide_legend(title="Gene Function"))+
    ggplot2::ggtitle(title)
  if(flip == TRUE){
    p <- p + coord_flip()
  }
  if(ReturnData){
    p <- lr_path
  }
  return(p)
}


#' @title plot_CCI_path2gene
#'
#' @description Draw a sankey plot containing pathways and genes.
#'
#' @param st_dec_data Deconvoluted Spatial matrix.
#' @param CCI_res The result of `find_CCI`, which contains data.frame including lrpair, tf, and cell_pair.
#' @param pathways Pathways with source and destination. Default is `SpaTalk::pathways`.
#' @param species Character. "Human" or "Mouse" or other specified. Default is "Human".
#' @param celltype_sender Character. The sender cell type.
#' @param celltype_receiver Character. The receiver cell type.
#' @param ligand Character. The expressed ligand of celltype_sender cells.
#' @param receptor Character. The expressed receptor of celltype_receiver cells.
#' @param max_hop Numeric. The maximum hop to plot. Default is `NULL`.
#' @param min_gene_num Numeric. The minimum number of genes in one pathways. Default is 5.
#' @param pvalue Numeric. The pvalue threshold of enrichment analysis. Default is 0.05.
#' @param geneRatio Numeric. The ratio threshold of enrichment analysis. Default is `NULL`.
#' @param color Character. The color to be used. Default is `NULL`.
#' @param size Numeric. The size of text. Default is 2.
#' @param title Character. The title of the plot. Default is `NULL`.
#' @param ReturnData Logical. Whether to return the enrichment result. Default is FALSE (sankey plot).
#'
#' @return A ggplot object (when `ReturnData` is FALSE) or a data.frame (when `ReturnData` is TRUE).
#'
#' @export
#'
#' @example plot_CCI_path2gene(st_dec_data = SpaTalk_dec_sc,CCI_res = CCI_res,
#'                   celltype_sender = "Fibroblasts", celltype_receiver = "Cancer_clone_A",
#'                   ligand = "DCN", receptor = "EGFR",max_hop = 1,geneRatio = 0.3,title = "DCN-EGFR")
#'
#' @author GM.W
#'
plot_CCI_path2gene <- function(st_dec_data, CCI_res,
                               pathways = NULL, species = "Human",
                               celltype_sender, celltype_receiver,
                               ligand, receptor, max_hop = NULL, min_gene_num = 5,
                               pvalue = 0.05, geneRatio = NULL, color = NULL, size = 2,
                               title = NULL, ReturnData = FALSE){
  # Check Input ---------------------------------
  lrpair <- CCI_res$lrpair
  TF <- CCI_res$tf
  if(is.null(pathways)){
    pathways <- SpaTalk::pathways
  }
  if(!species %in% pathways$species){
    stop("Please provide the correct name of species!")
  }
  pathways <- pathways[pathways$species == species,]
  if (!celltype_sender %in% lrpair$celltype_sender) {
    stop("Please provide the correct name of celltype_sender!")
  }
  if (!celltype_receiver %in% lrpair$celltype_receiver) {
    stop("Please provide the correct name of celltype_receiver!")
  }
  if (!ligand %in% lrpair$ligand) {
    stop("Please provide the correct name of ligand!")
  }
  if (!receptor %in% lrpair$receptor) {
    stop("Please provide the correct name of receptor!")
  }
  check <- lrpair[lrpair$ligand == ligand & lrpair$receptor == receptor,]
  if(nrow(check) == 0){
    stop("Please provide the correct ligand-receptor pair!")
  }
  if(!is.null(max_hop)){
    if(max_hop > max(TF$n_hop)){
      stop("Please provide the max_hop smaller than you used in `filter_lr_paths`!")
    }
  }

  # Data Organization ---------------------------
  TF <- TF[which(TF$receptor == receptor),]
  if(!is.null(max_hop)){
    TF <- TF[which(TF$n_hop <= max_hop),]
  }
  ### lr_path
  # L_R
  L_R <- data.frame(src = ligand, dest = receptor, hop = 0)
  # R_TF
  R_TF_tmp <- TF$mediator
  R_TF_tmp <- as.data.frame(stringr::str_split(R_TF_tmp, ",", simplify = T))
  R_TF <- NULL
  for (i in 1:nrow(R_TF_tmp)) {
    R_TF_1 <- as.character(R_TF_tmp[i,])
    R_TF_1 <- R_TF_1[which(R_TF_1 != "")]
    R_TF <- c(R_TF, R_TF_1)
  }
  R_TF <- as.data.frame(stringr::str_split(R_TF, "->|__",simplify = T))
  colnames(R_TF) <- c("src","dest","hop")
  R_TF <- unique(R_TF)
  # TF_Target
  TF_Target_tmp <- TF[, c("tf","n_hop","target")]
  TF_Target_genes <- as.data.frame(stringr::str_split(TF_Target_tmp$target, ",", simplify = T))
  TF_Target <- data.frame()
  for (i in 1:nrow(TF_Target_tmp)) {
    gene1 <- as.character(TF_Target_genes[i,])
    gene1 <- gene1[which(gene1 != "")]
    TF_Target_1 <- data.frame(src = TF_Target_tmp[i,"tf"], dest = gene1, hop = TF_Target_tmp[i,"n_hop"]+1)
    TF_Target <- rbind(TF_Target, TF_Target_1)
  }
  TF_Target <- unique(TF_Target)
  lr_path <- rbind(L_R, R_TF, TF_Target)
  # genes
  genes_interest <- unique(c(lr_path$src, lr_path$dest))
  num_genes_interest <- length(genes_interest)
  genes_all <- rownames(st_dec_data)
  num_genes_all <- length(genes_all)

  pathways_rec <- pathways[which(pathways$src == receptor | pathways$dest == receptor),]

  # For result data.frame
  pathway_names <- unique(pathways_rec$pathway)
  pathway_state <- rep("No", length(pathway_names))
  geneID <- rep("", length(pathway_names))
  pathway_pvalue <- as.double(rep(1, length(pathway_names)))
  gene_ratio <- rep(0, length(pathway_names))
  gene_ratio_str <- rep("", length(pathway_names))
  bg_ratio_str <- rep("", length(pathway_names))
  count_intersect_gene <- rep(0, length(pathway_names))

  # Hypergeometric distribution (fisher.test)
  for (i in 1:length(pathway_names)) {
    pathways_1 <- pathways[which(pathways$pathway == pathway_names[i]),]
    genes_pathway <- unique(c(pathways_1$src, pathways_1$dest))
    genes_pathway <- genes_pathway[genes_pathway %in% genes_all]
    num_genes_pathway <- length(genes_pathway)
    if(num_genes_pathway >= min_gene_num){
      genes_pathway_interest <- intersect(genes_interest, genes_pathway)
      num_genes_pathway_interest <- length(genes_pathway_interest)
      mtx <- matrix(c(num_genes_pathway_interest, num_genes_interest - num_genes_pathway_interest,
                      num_genes_pathway - num_genes_pathway_interest,
                      num_genes_all - num_genes_interest - num_genes_pathway + num_genes_pathway_interest),
                    nrow = 2,byrow = F)
      fisher_pvalue <- fisher.test(mtx)
      fisher_pvalue <- as.double(fisher_pvalue$p.value)

      # collect result
      pathway_state[i] <- "Yes"
      geneID[i] <- paste0(genes_pathway_interest,collapse = "/")
      pathway_pvalue[i] <- fisher_pvalue
      gene_ratio[i] <- num_genes_pathway_interest / num_genes_interest
      gene_ratio_str[i] <- paste0(num_genes_pathway_interest, "/", num_genes_interest, collapse = "")
      bg_ratio_str[i] <- paste0(num_genes_pathway, "/", num_genes_all, collapse = "")
      count_intersect_gene[i] <- num_genes_pathway_interest
    }
  }
  res <- data.frame(pathway = pathway_names, state = pathway_state, GeneRatio = gene_ratio,
                    GeneRatio_str = gene_ratio_str, BgRatio_str = bg_ratio_str,
                    pvalue = pathway_pvalue, geneID = geneID, Count = count_intersect_gene)
  res <- res[which(res$pvalue < pvalue & res$state == "Yes"),]
  if(!is.null(geneRatio)){
    res <- res[which(res$GeneRatio > geneRatio),]
  }

  # sankey_df
  sankey_df <- data.frame()
  for (i in 1:nrow(res)) {
    path_1 <- res$pathway[i]
    genes_1 <- res$geneID[i]
    genes_1 <- as.character(stringr::str_split(genes_1,"/",simplify = T))
    sankey_tmp <- data.frame(pathway = path_1, gene = genes_1)
    sankey_df <- rbind(sankey_df, sankey_tmp)
  }
  suppressWarnings(suppressMessages(library(dplyr)))
  suppressMessages(
    sankey_freq <- sankey_df %>%
      group_by(pathway, gene) %>%
      summarise(count = n())
  )
  sankey_long <- ggalluvial::to_lodes_form(data = sankey_freq, key = "x_label", axes = 1:2)

  # color
  if(is.null(color)){
    color <- ggpubr::get_palette(palette = "Set3", k=2)
  }

  # Plot ----------------------------------------
  require(ggplot2)
  require(ggalluvial)
  p <- ggplot2::ggplot(sankey_long,
                       ggplot2::aes(x=x_label,stratum=stratum,alluvium=alluvium,
                                    y=count,fill=x_label,label=stratum))+
    ggalluvial::geom_flow()+
    ggalluvial::geom_stratum()+
    ggplot2::geom_text(stat = "stratum",size = size)+
    ggplot2::labs(x="", y="")+
    ggplot2::theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
          panel.background = element_rect(colour = "white", fill="white"),
          plot.background = element_rect(colour = "white", fill="white"),
          panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          legend.title= element_text(size = 10,face="bold"),
          legend.text= element_text(size = 8),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          legend.key.size = unit(0.45, 'cm'),
          strip.text = element_text(size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,face = "bold"))+
    ggplot2::scale_fill_manual(values = color)+
    ggplot2::ggtitle(title)

  # Return --------------------------------------
  if(ReturnData){
    p <- res
  }
  return(p)
}
