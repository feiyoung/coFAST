
# Cluster assessment --------------------------------------------------

#' Add the spatial coordinates to the reduction slot
#' @description Calculate the adjacency matrix given a spatial coordinate matrix with 2-dimension or 3-dimension or more.
#' @param seu a SeuratObject with spatial coordinate information in the \code{meta.data} slot.
#' @param coord.name a character vector, specify the names of spatial coordinates in the \code{meta.data} slot. For example, \code{c("x", "y")}.
#' @param assay a string, specify the assay.
#' @return return a revised Seurat object with a slot `Spatial` in the reductions slot.
#'
#' @seealso None
#' @references None
#' @export
#' @importFrom Seurat CreateDimReducObject
#' @examples
#' data(CosMx_subset)
#' Addcoord2embed(CosMx_subset, coord.name = c("x", "y"))
#'
#'
Addcoord2embed <- function(seu, coord.name,  assay = "RNA"){

  # coord.name <- c("x", "y")
  coord <- seu@meta.data[,coord.name]
  row.names(coord) <- colnames(seu)
  colnames(coord) <- paste0(coord, 1:ncol(coord))
  seu@reductions[['Spatial']] <- CreateDimReducObject(embeddings = as.matrix(coord),
                                                       key = paste0(toupper("Spatial"), "_"), assay = assay)
  seu
}


#' Calculate the aggregation score for specific clusters
#' @description Calculate the adjacency matrix given a spatial coordinate matrix with 2-dimension or 3-dimension or more.
#' @param seu  a SeuratObject with reductions not \code{NULL}.
#' @param reduction.name an character, specify the reduction name for calculating the aggregation score.
#' @param random.seed a positive integer, specify the random seed for reproducibility.
#' @return return a data.frame with two columns: the first column is the number of spots in each category (cluster/cell type); the second column is the corresponding aggregation score.
#'
#' @seealso None
#' @references None
#' @export
#' @importFrom Seurat Idents
#' @examples
#' data(CosMx_subset)
#' library(Seurat)
#' CosMx_subset <- Addcoord2embed(CosMx_subset, coord.name = c("x", "y"))
#' Idents(CosMx_subset) <- 'cell_type'
#' dat.sp.score <- AggregationScore(CosMx_subset, reduction.name = 'Spatial')
#' print(dat.sp.score)
#'
AggregationScore <- function(seu, reduction.name='cofast', random.seed=1){

  subtype <- Idents(seu)
  idents_levels <- levels(subtype)
  n_celltypes <- length(idents_levels)
  cell.use <- colnames(seu)
  cell.set.list <- lapply(idents_levels, function(x) cell.use[which(subtype == x)])
  set.seed(random.seed)
  res.rgTest2 <- pathway.rgTest(
    Embeddings(seu, reduction = reduction.name), cell.use, cell.set.list, gene.set.cutoff = 3,
    test.type = list("ori", "gen", "wei", "max"), k = 5, wei.fun = "weiMax",
    perm.num = 1000, ncores = 20, parallel = TRUE, asy_threshold = NULL
  )
  spot_vec <- as.vector(table(subtype))
  res.rgTest2$num.spots <- spot_vec[spot_vec >= 3]
  row.names(res.rgTest2) <- idents_levels[spot_vec >= 3]

  ## save rgTest result for cell embedding and pos
  columns.use <- c("num.spots", "asy.max.statistic")
  rgtest_max_res <- res.rgTest2[,columns.use]
  rgtest_max_res <- as.data.frame(rgtest_max_res)
  colnames(rgtest_max_res) <- c("num.spots", 'AggScore')
  return(rgtest_max_res)
}


pathway.rgTest <- function (coembed, genes.use, gene.set.list, gene.set.cutoff = 3,
                            test.type = list("ori", "gen", "wei", "max"), k = 5, wei.fun = c("weiMax",
                                                                                             "weiGeo", "weiArith"), perm.num = 0, progress_bar = FALSE,
                            ncores = 20, eta = 1e-04, parallel = TRUE, asy_threshold = 100){
  wei.fun <- match.arg(wei.fun)
  coembed.scale <- scale(coembed[genes.use, ])
  gene.set.list.cut <- lapply(gene.set.list, function(gene.set) intersect(gene.set,
                                                                          genes.use))
  l.gene.set <- sapply(gene.set.list.cut, length)
  gene.set.list.cut <- gene.set.list.cut[l.gene.set >= gene.set.cutoff]
  E <- kmst(y = coembed.scale, k = k)
  if (is.null(asy_threshold)) {
    res <- pathway.rgTest.body.full(E = E, genes.use = genes.use,
                                    gene.set.list.cut = gene.set.list.cut, test.type = test.type,
                                    wei.fun = wei.fun, perm.num = perm.num, progress_bar = progress_bar,
                                    ncores = ncores, eta = eta, parallel = parallel)
  }
  else {
    res <- pathway.rgTest.body(E = E, genes.use = genes.use,
                               gene.set.list.cut = gene.set.list.cut, test.type = test.type,
                               wei.fun = wei.fun, perm.num = perm.num, progress_bar = progress_bar,
                               ncores = ncores, eta = eta, parallel = parallel,
                               asy_threshold = asy_threshold)
  }
  return(res)
}

#' @importFrom ade4 mstree
#' @importFrom stats as.dist
kmst <- function(y=NULL, dis=NULL, k=1){
  require(ade4)
  if (is.null(dis) && is.null(y)){
    cat("Please input data or the distance matrix!\n")
    return(0)
  }
  if (is.null(dis)) dis = pdistance.mat(y, y)
  mymst = ade4::mstree(stats::as.dist(dis),k)
  cbind(mymst[,1], mymst[,2])
}
#' @useDynLib coFAST, .registration = TRUE
pdistance.mat <- function (Ar, Br, eta = 1e-10){
  dis <- pdistance_cpp(Ar, Br, eta = eta)
  rownames(dis) <- rownames(Ar)
  colnames(dis) <- rownames(Br)
  return(dis)
}

#' @importFrom dplyr `%>%`
pathway.rgTest.body.full <- function (E, genes.use, gene.set.list.cut, test.type = list("ori",
                                                                                        "gen", "wei", "max"), wei.fun = c("weiMax", "weiGeo", "weiArith"),
                                      perm.num = 0, progress_bar = FALSE, ncores = 10, eta = 1e-04,
                                      parallel = TRUE){
  wei.fun <- match.arg(wei.fun)
  l.gene.set <- sapply(gene.set.list.cut, length)
  n_gene <- length(genes.use)
  n_pathway <- length(gene.set.list.cut)
  ii <- c(E[, 1], E[, 2])
  jj <- c(E[, 2], E[, 1])
  nodedeg <- as.integer(table(factor(ii, levels = 1:n_gene)))
  if (wei.fun == "weiMax") {
    wei_sp_mat <- Matrix::sparseMatrix(i = ii, j = jj, x = 1/pmax(nodedeg[ii],
                                                                  nodedeg[jj]))
  }
  else if (wei.fun == "weiGeo") {
    wei_sp_mat <- Matrix::sparseMatrix(i = ii, j = jj, x = 2/(nodedeg[ii] +
                                                                nodedeg[jj]))
  }
  else {
    wei_sp_mat <- Matrix::sparseMatrix(i = ii, j = jj, x = 1/sqrt(nodedeg[ii] *
                                                                    nodedeg[jj]))
  }
  match_idx <- 1:n_gene
  names(match_idx) <- genes.use
  idx.list <- lapply(gene.set.list.cut, function(gene.set) match_idx[gene.set])
  test_fun <- function(idx) {
    n1 <- length(idx)
    idx <- c(idx, setdiff(1:n_gene, idx))
    test_res <- fast.rg.test.full(wei_sp_mat = wei_sp_mat[idx,
                                                          idx], n1 = n1, test.type = test.type, perm.num = perm.num,
                                  progress_bar = progress_bar, eta = eta)
    unlist(test_res)
  }
  if (parallel) {
    library(furrr)
    library(future)
    plan("multicore", workers = ncores)
    pb <- progress::progress_bar$new(format = "[:bar] :percent ETA: :eta")
    res <- idx.list %>% future_map(test_fun, .progress = TRUE,
                                   .options = furrr_options(seed = TRUE))
    res <- t(Reduce(cbind, res))
  }
  else {
    res <- t(pbapply::pbsapply(idx.list, test_fun))
  }
  rownames(res) <- names(gene.set.list.cut)
  res <- cbind(res, gene.set.length = l.gene.set)
  p.adj <- apply(res[, paste0("asy.", unlist(test.type), ".pval"),
                     drop = FALSE], 2, p.adjust, method = "fdr")
  colnames(p.adj) <- paste0("asy.", unlist(test.type), ".pval.adj")
  res <- as.data.frame(cbind(res, p.adj))
  if (perm.num > 0) {
    p.adj <- apply(res[, paste0("perm.", unlist(test.type),
                                ".pval"), drop = FALSE], 2, p.adjust, method = "fdr")
    colnames(p.adj) <- paste0("perm.", unlist(test.type),
                              ".pval.adj")
    res <- as.data.frame(cbind(res, p.adj))
  }
  return(res)
}

pathway.rgTest.body <- function (E, genes.use, gene.set.list.cut, test.type = list("ori",
                                                                                   "gen", "wei", "max"), wei.fun = c("weiMax", "weiGeo", "weiArith"),
                                 perm.num = 0, progress_bar = FALSE, ncores = 10, eta = 1e-04,
                                 parallel = TRUE, asy_threshold = 100)
{
  wei.fun <- match.arg(wei.fun)
  l.gene.set <- sapply(gene.set.list.cut, length)
  n_gene <- length(genes.use)
  n_pathway <- length(gene.set.list.cut)
  ii <- c(E[, 1], E[, 2])
  jj <- c(E[, 2], E[, 1])
  nodedeg <- as.integer(table(factor(ii, levels = 1:n_gene)))
  if (wei.fun == "weiMax") {
    wei_sp_mat <- Matrix::sparseMatrix(i = ii, j = jj, x = 1/pmax(nodedeg[ii],
                                                                  nodedeg[jj]))
  }
  else if (wei.fun == "weiGeo") {
    wei_sp_mat <- Matrix::sparseMatrix(i = ii, j = jj, x = 2/(nodedeg[ii] +
                                                                nodedeg[jj]))
  }
  else {
    wei_sp_mat <- Matrix::sparseMatrix(i = ii, j = jj, x = 1/sqrt(nodedeg[ii] *
                                                                    nodedeg[jj]))
  }
  match_idx <- 1:n_gene
  names(match_idx) <- genes.use
  idx.list <- lapply(gene.set.list.cut, function(gene.set) match_idx[gene.set])
  test_fun <- function(idx) {
    n1 <- length(idx)
    idx <- c(idx, setdiff(1:n_gene, idx))
    test_res <- fast.rg.test(wei_sp_mat = wei_sp_mat[idx,
                                                     idx], n1 = n1, test.type = test.type, perm.num = perm.num,
                             progress_bar = progress_bar, eta = eta, asy_threshold = asy_threshold)
    unlist(test_res)
  }
  if (parallel) {
    library(furrr)
    library(future)
    plan("multicore", workers = ncores)
    pb <- progress::progress_bar$new(format = "[:bar] :percent ETA: :eta")
    res <- idx.list %>% future_map(test_fun, .progress = TRUE,
                                   .options = furrr_options(seed = TRUE))
    res <- t(Reduce(cbind, res))
  }
  else {
    res <- t(pbapply::pbsapply(idx.list, test_fun))
  }
  rownames(res) <- names(gene.set.list.cut)
  res <- cbind(res, gene.set.length = l.gene.set)
  p.adj <- apply(res[, paste0(unlist(test.type), ".pval"),
                     drop = FALSE], 2, p.adjust, method = "fdr")
  colnames(p.adj) <- paste0(unlist(test.type), ".pval.adj")
  res <- as.data.frame(cbind(res, p.adj))
  return(res)
}
#' @importFrom stats pnorm
fast.rg.test.full <- function(
    wei_sp_mat, n1, test.type = list("ori", "gen", "wei", "max"),
    perm.num = 0, progress_bar = FALSE, eta = 1e-4){
  n = nrow(wei_sp_mat)
  n2 <- n - n1

  if(perm.num != 0 && !is.numeric(perm.num)){
    stop('incorrect perm_num')
  }

  weighted_edge_count <- fast_weighted_R1R2(wei_sp_mat, 1:n1)
  # browser()
  asy_res <- fast_theo_mu_sig(wei_sp_mat, n1, n2)

  res_list <- list()
  if("ori" %in% test.type){
    #original edge count test
    ori_test <- (weighted_edge_count$R - asy_res$mu)/sqrt(asy_res$sig)
    ori_asy_p <- pnorm(ori_test)
    res_list <- c(res_list, asy.ori.statistic = ori_test, asy.ori.pval = ori_asy_p)
  }
  if("gen" %in% test.type){
    #generalized edge count test
    gen_asy_p <- asy_gen(asy_res, weighted_edge_count$R1, weighted_edge_count$R2)
    res_list <- c(res_list, asy.gen.statistic = gen_asy_p$test_statistic, asy.gen.pval = gen_asy_p$p_value)
  }
  if("wei" %in% test.type){
    #weighted edge count test
    wei_asy_p <- asy_wei(asy_res, weighted_edge_count$R1, weighted_edge_count$R2, n1, n2)
    res_list <- c(res_list, asy.wei.statistic = wei_asy_p$test_statistic, asy.wei.pval = wei_asy_p$p_value)
  }
  if("max" %in% test.type){
    #max-type edge count test
    max_asy_p <- asy_max(asy_res, weighted_edge_count$R1, weighted_edge_count$R2, n1, n2)
    res_list <- c(res_list, asy.max.statistic = max_asy_p$test_statistic, asy.max.pval = max_asy_p$p_value)
  }

  if (perm.num != 0 && is.numeric(perm.num)) {
    if (perm.num > 10000) {
      message("Doing more than 10,000 permutations could be time consuming.\n")
    }
    n_per = perm.num
    # temp_list = permu_edge(n_per, weis, progress_bar)
    temp_list = fast.permu_edge(n_per, wei_sp_mat, n1, n2, progress_bar)
    if("ori" %in% test.type){
      #original edge count test
      ori = sum(temp_list$R < weighted_edge_count$R)
      ori = ori/n_per
      res_list = c(res_list, perm.ori.pval = ori)
    }
    if("gen" %in% test.type){
      #generalized edge count test
      gen_p_res <- fast.permu_gen(temp_list$R1, temp_list$R2, weighted_edge_count$R1, weighted_edge_count$R2, n_per, eta)
      res_list <- c(res_list, perm.gen.pval = gen_p_res$pval)
    }
    if("wei" %in% test.type){
      #weighted edge count test
      wei_p_res <- fast.permu_wei(temp_list$R1, temp_list$R2, weighted_edge_count$R1, weighted_edge_count$R2, n1, n2, n_per)
      res_list <- c(res_list, perm.wei.pval = wei_p_res$pval)
    }
    if("max" %in% test.type){
      #max-type edge count test
      max_p_res <- fast.permu_max(temp_list$R1, temp_list$R2, weighted_edge_count$R1, weighted_edge_count$R2, n1, n2, n_per)
      res_list <- c(res_list, perm.max.pval = max_p_res$pval)
    }
  }

  return(res_list)
}

fast.rg.test <- function(
    wei_sp_mat, n1, test.type = list("ori", "gen", "wei", "max"),
    perm.num = 0, progress_bar = FALSE, eta = 1e-4, asy_threshold = 100){
  n = nrow(wei_sp_mat)
  n2 <- n - n1

  if(perm.num != 0 && !is.numeric(perm.num)){
    stop('incorrect perm_num')
  }

  weighted_edge_count <- fast_weighted_R1R2(wei_sp_mat, 1:n1)
  # browser()
  asy_res <- fast_theo_mu_sig(wei_sp_mat, n1, n2)

  res_list <- list()
  if("ori" %in% test.type){
    #original edge count test
    ori_test <- (weighted_edge_count$R - asy_res$mu)/sqrt(asy_res$sig)
    ori_asy_p <- pnorm(ori_test)
    res_list <- c(res_list, asy.ori.statistic = ori_test, asy.ori.pval = ori_asy_p)
  }
  if("gen" %in% test.type){
    #generalized edge count test
    gen_asy_p <- asy_gen(asy_res, weighted_edge_count$R1, weighted_edge_count$R2)
    res_list <- c(res_list, asy.gen.statistic = gen_asy_p$test_statistic, asy.gen.pval = gen_asy_p$p_value)
  }
  if("wei" %in% test.type){
    #weighted edge count test
    wei_asy_p <- asy_wei(asy_res, weighted_edge_count$R1, weighted_edge_count$R2, n1, n2)
    res_list <- c(res_list, asy.wei.statistic = wei_asy_p$test_statistic, asy.wei.pval = wei_asy_p$p_value)
  }
  if("max" %in% test.type){
    #max-type edge count test
    max_asy_p <- asy_max(asy_res, weighted_edge_count$R1, weighted_edge_count$R2, n1, n2)
    res_list <- c(res_list, asy.max.statistic = max_asy_p$test_statistic, asy.max.pval = max_asy_p$p_value)
  }

  if (perm.num != 0 && is.numeric(perm.num) && n1 >= asy_threshold) {
    if (perm.num > 10000) {
      message("Doing more than 10,000 permutations could be time consuming.\n")
    }
    n_per = perm.num
    # temp_list = permu_edge(n_per, weis, progress_bar)
    temp_list = fast.permu_edge(n_per, wei_sp_mat, n1, n2, progress_bar)
    if("ori" %in% test.type){
      #original edge count test
      ori = sum(temp_list$R < weighted_edge_count$R)
      ori = ori/n_per
      res_list = c(res_list, perm.ori.pval = ori)
    }
    if("gen" %in% test.type){
      #generalized edge count test
      gen_p_res <- fast.permu_gen(temp_list$R1, temp_list$R2, weighted_edge_count$R1, weighted_edge_count$R2, n_per, eta)
      res_list <- c(res_list, perm.gen.pval = gen_p_res$pval)
    }
    if("wei" %in% test.type){
      #weighted edge count test
      wei_p_res <- fast.permu_wei(temp_list$R1, temp_list$R2, weighted_edge_count$R1, weighted_edge_count$R2, n1, n2, n_per)
      res_list <- c(res_list, perm.wei.pval = wei_p_res$pval)
    }
    if("max" %in% test.type){
      #max-type edge count test
      max_p_res <- fast.permu_max(temp_list$R1, temp_list$R2, weighted_edge_count$R1, weighted_edge_count$R2, n1, n2, n_per)
      res_list <- c(res_list, perm.max.pval = max_p_res$pval)
    }

    output <- c(
      res_list[paste0("asy.", test.type, ".statistic")],
      res_list[paste0("perm.", test.type, ".pval")])
  } else {
    output <- c(
      res_list[paste0("asy.", test.type, ".statistic")],
      res_list[paste0("asy.", test.type, ".pval")])
  }

  names(output) <- c(paste0(test.type, ".statistic"), paste0(test.type, ".pval"))

  return(output)
}


fast_weighted_R1R2 <- function(wei_sp_mat, G1){
  R1 = sum(wei_sp_mat[G1, G1]) * 0.5
  R2 = sum(wei_sp_mat[-G1, -G1]) * 0.5
  R = sum(wei_sp_mat[G1, -G1])

  return(list(R1=R1, R2=R2, R=R))
}

fast_theo_mu_sig <- function(wei_sp_mat, n1, n2){
  N = n1+n2

  sum_wei_2 <- sum(wei_sp_mat@x^2) * 0.5
  sum_wei <- sum(wei_sp_mat@x) * 0.5

  mu_1 = sum_wei*n1*(n1-1)/N/(N-1)
  mu_2 = sum_wei*n2*(n2-1)/N/(N-1)
  mu = 2*sum_wei*n1*n2/N/(N-1)

  part11 = sum_wei_2*n1*(n1-1)/N/(N-1)
  part21 = sum_wei_2*n2*(n2-1)/N/(N-1)
  part1 = 2*sum_wei_2*n1*n2/N/(N-1)
  #edge pair

  ss = sum(Matrix::rowSums(wei_sp_mat)^2) * 0.5 - sum_wei_2
  part12 = 2*ss*n1*(n1-1)*(n1-2)/N/(N-1)/(N-2)
  part22 = 2*ss*n2*(n2-1)*(n2-2)/N/(N-1)/(N-2)
  part2 = 2*ss*(n1*n2*(n2-1)/N/(N-1)/(N-2) + n1*n2*(n1-1)/N/(N-1)/(N-2))
  #different edges
  s_4 = sum_wei^2
  s_4 = s_4 - 2*ss - sum_wei_2

  part13 = s_4*n1*(n1-1)*(n1-2)*(n1-3)/N/(N-1)/(N-2)/(N-3)
  part23 = s_4*n2*(n2-1)*(n2-2)*(n2-3)/N/(N-1)/(N-2)/(N-3)
  part3 = s_4*4*n1*n2*(n1-1)*(n2-1)/N/(N-1)/(N-2)/(N-3)

  sig1 = part11+part12+part13 - mu_1^2
  sig2 = part21+part22+part23 - mu_2^2
  sig = part1+part2+part3 - mu^2

  part_cov = s_4*n1*(n1-1)*n2*(n2-1)/N/(N-1)/(N-2)/(N-3)
  sig_12 = part_cov - mu_1*mu_2
  return(list(mu = mu, mu1 = mu_1, mu2 = mu_2, sig = sig, sig11 = sig1, sig22 = sig2, sig12 = sig_12))
}

fast.permu_edge <- function(n_per, wei_sp_mat, n1, n2, progress_bar = FALSE){
  n = n1+n2
  if(progress_bar){
    pb = txtProgressBar(min = 0, max = n_per, initial = 0)
    temp = sapply(1:n_per, function(peri){
      setTxtProgressBar(pb,peri)
      edgeinfo = fast_weighted_R1R2(wei_sp_mat, sample(1:n, n1))
      return(unlist(edgeinfo))
    })
    close(pb)
    distri_R1 = temp['R1', ]
    distri_R2 = temp['R2', ]
    distri_R = temp['R', ]
  }else{
    temp = sapply(1:n_per, function(peri){
      edgeinfo = fast_weighted_R1R2(wei_sp_mat, sample(1:n, n1))
      unlist(edgeinfo)
    })

    distri_R1 = temp['R1', ]
    distri_R2 = temp['R2', ]
    distri_R = temp['R', ]
  }
  return(list(R1=distri_R1, R2=distri_R2, R=distri_R))
}

#' @importFrom stats cov
fast.permu_gen <- function(R1_list, R2_list, R1_test, R2_test, n_per, eta = 1e-4){
  mu1 = mean(R1_list)
  mu2 = mean(R2_list)
  sigma = cov(cbind(R1_list, R2_list))
  if (sigma[1,1]*sigma[2,2] - sigma[1,2]*sigma[2,1] < eta) {
    sigma_inverse <- solve(sigma + diag(2)*eta)
  } else {
    sigma_inverse <- solve(sigma)
  }
  temp = c(R1_test - mu1, R2_test - mu2)
  gen_t = sum((temp %*% sigma_inverse) * temp)

  null_gen = cbind(R1_list - mu1, R2_list - mu2)
  null_gen = rowSums((null_gen %*% sigma_inverse) * null_gen)

  pval <- mean(null_gen > gen_t)
  # tmp <- log(null_gen)
  # an.pval <- pnorm(log(gen_t), mean(tmp), sd(tmp), lower.tail = FALSE)
  return(
    list(
      # an.pval = an.pval,
      pval = pval
    )
  )
}

#' @importFrom stats cov
fast.permu_wei <- function(R1_list, R2_list, R1_test, R2_test, n1, n2, n_per){
  mu1 = mean(R1_list)
  mu2 = mean(R2_list)
  sigma = cov(cbind(R1_list, R2_list))
  p_w = (n1-1)/(n1+n2-2)
  q_w = 1 - p_w

  null_wei = rep(0, n_per)
  deno2 = sqrt(q_w^2*sigma[1, 1] + p_w^2*sigma[2, 2] + 2*p_w*q_w*sigma[1, 2])
  num2 <- q_w*mu1 + p_w*mu2
  w_t = ((q_w*R1_test + p_w*R2_test) - num2) / deno2
  null_wei = ((q_w*R1_list + p_w*R2_list) - num2)/deno2

  pval <- mean(null_wei > w_t)
  # tmp <- log(q_w*R1_list + p_w*R2_list)
  # an.pval <- pnorm(log(q_w*R1_test + p_w*R2_test), mean(tmp), sd(tmp), lower.tail = FALSE)
  return(
    list(
      # an.pval = an.pval,
      pval = pval
    )
  )
}

#' @importFrom stats cov
fast.permu_max <- function(R1_list, R2_list, R1_test, R2_test, n1, n2, n_per){
  mu1 = mean(R1_list)
  mu2 = mean(R2_list)
  sigma = cov(cbind(R1_list, R2_list))
  p_w = (n1-1)/(n1+n2-2)
  q_w = 1 - p_w

  deno1 = sqrt(sigma[1, 1] + sigma[2, 2] - 2*sigma[1, 2])
  deno2 = sqrt(q_w^2*sigma[1, 1] + p_w^2*sigma[2, 2] + 2*p_w*q_w*sigma[1, 2])
  max_diff_t = ((R1_test - R2_test) - (mu1 - mu2)) / deno1
  max_w_t = ((q_w*R1_test + p_w*R2_test) - (q_w*mu1 + p_w*mu2)) / deno2
  max_t = max(abs(max_diff_t), max_w_t)
  null_max = pmax(abs(((R1_list - R2_list) - (mu1 - mu2))/deno1), ((q_w*R1_list + p_w*R2_list) - (q_w*mu1 + p_w*mu2))/deno2)

  pval <- mean(null_max > max_t)
  # tmp <- log(null_max)
  # an.pval <- pnorm(log(max_t), mean(tmp), sd(tmp), lower.tail = FALSE)
  return(
    list(
      # an.pval = an.pval,
      pval = pval
    )
  )
}
#' @importFrom stats pnorm
asy_max <- function(asy_res, R1_test, R2_test, n1, n2){
  max_diff_t = ((R1_test - R2_test) - (asy_res$mu1 - asy_res$mu2)) /
    (sqrt(asy_res$sig11 + asy_res$sig22 - 2*asy_res$sig12))
  p_w = (n1-1)/(n1+n2-2)
  q_w = 1 - p_w
  max_w_t = ((q_w*R1_test + p_w*R2_test) - (q_w*asy_res$mu1 + p_w*asy_res$mu2)) /
    (sqrt(q_w^2*asy_res$sig11 + p_w^2*asy_res$sig22 + 2*p_w*q_w*asy_res$sig12))
  max_t = max(abs(max_diff_t), max_w_t)
  p_max = 1-2*(pnorm(max_t)-0.5)*pnorm(max_t)
  test_statistic = max_t
  return(list(test_statistic = test_statistic, p_value = p_max))
}

#' @importFrom stats pnorm
asy_wei <- function(asy_res, R1_test, R2_test, n1, n2){
  p_w = (n1-1)/(n1+n2-2)
  q_w = 1 - p_w
  max_w_t = ((q_w*R1_test + p_w*R2_test) - (q_w*asy_res$mu1 + p_w*asy_res$mu2)) /
    (sqrt(q_w^2*asy_res$sig11 + p_w^2*asy_res$sig22 + 2*p_w*q_w*asy_res$sig12))
  test_statistic = max_w_t
  return(list(test_statistic = test_statistic, p_value = 1-pnorm(max_w_t)))
}
#' @importFrom stats pchisq
asy_gen <- function(asy_res, R1_test, R2_test){
  sigma = matrix(c(asy_res$sig11, asy_res$sig12, asy_res$sig12, asy_res$sig22), nrow = 2)
  temp_asy = c(R1_test - asy_res$mu1, R2_test - asy_res$mu2) %*% solve(sigma)
  z_gen = temp_asy %*% c(R1_test - asy_res$mu1, R2_test - asy_res$mu2)
  test_statistic = z_gen[1, 1]
  return(list(test_statistic = test_statistic, p_value = 1-pchisq(test_statistic, df = 2)))
}
