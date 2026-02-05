# functions to run robust PIC regression with permulations and early stopping

#' Run permulations given a phenotype and a tree
#'
#' @param pheno Named numeric vector giving phenotype values per species.
#' @param tree A tree object.
#' @param n Number of permulations to generate.
#' @param p_model Optionally, a model of evolution other than Brownian motion (used to transform the tree; see `phylolm::transf.branch.lengths()`).
#' @param ub Upper bound on parameter in p_model.
#' @param lb Lower bound on parameter in p_model.
#' @param pheno_sd Optionally, a vector of standard deviations around each measurement in `pheno`.
#' @export permulate
permulate <- function(
  pheno,
  tree,
  n = 1,
  p_model = "BM",
  ub = 0.99,
  lb = 0.01,
  pheno_sd = NULL,
  use_density = FALSE
) {
  rbt_params <- get_best_model_fit(
    pheno,
    tree,
    p_model,
    ub = ub,
    lb = lb,
    yield = "params"
  )
  random_traits <- phylolm::rTrait(
    n,
    tree,
    model = p_model,
    parameters = rbt_params
  )
  if (n == 1) {
    random_traits <- data.frame(random_traits)
  }
  if (use_density) {
    bandwidth <- bw.SJ(pheno)
  }
  rt <- apply(random_traits, 2, \(x) {
    if (use_density) {
      resampled_pheno <- sample(pheno, replace = TRUE) +
      rnorm(length(pheno), mean = 0, sd = bandwidth)
      tr <- sort(resampled_pheno)[rank(x, ties.method = "random")]
    } else {
      tr <- sort(pheno)[rank(x, ties.method = "random")]
    }
    names(tr) <- names(x)
    if (!is.null(pheno_sd)) {
      tr_sd <- pheno_sd[names(sort(pheno))[rank(x, ties.method = "random")]]
      names(tr_sd) <- names(x)
      return(list(tr = tr, s = tr_sd))
    } else {
      return(tr)
    }
  })
  if (!is.null(pheno_sd)) {
    return(list(
      tr = Reduce(cbind, lapply(rt, \(x) x$tr)),
      s = Reduce(cbind, lapply(rt, \(x) x$s))
    ))
  } else {
    return(rt)
  }
}

#' Run rate permutations given a phenotype and a tree
#' @param pheno Named numeric vector giving phenotype values per species.
#' @param tree A tree object.
#' @param n Number of permulations to generate.
#' @param n_bins Number of bins to divide edges into.
#' @param p_model Optionally, a model of evolution other than Brownian motion (used to transform the tree; see `phylolm::transf.branch.lengths()`).
#' @param ub Upper bound on parameter in p_model.
#' @param lb Lower bound on parameter in p_model.
#' @param cut_by How to bin edges (default: "parent_depth"; can also be "edge_length". Any other value bins by child depth)
#' @param rank Transform the permutrated values back to the original data points, as in permulation (default: FALSE).
#' @export permutrate
permutrate <- function(pheno, tree, n = 1, nbins = 10, p_model = "BM", ub=0.99, lb=0.01, cut_by="parent_depth", rank=FALSE) {
  # p_model is currently ignored since it seems to mess up more than it helps
  # if (!p_model == "BM") {
  #   tree <- get_best_model_fit(
  #     pheno, tree, p_model, ub = ub, lb = lb, yield = "tree"
  #   )
  # }
  tree_r <- castor::reorder_tree_edges(
    tree,
    root_to_tips = TRUE,
    depth_first_search = TRUE
  )
  perm_mtx <- matrix(nr = length(pheno), nc = n)
  node_tip_vals <- c(
    pheno,
    castor::asr_independent_contrasts(tree_r, pheno)$ancestral_states
  )
  tip_node_depths <- c(
    rep(0, length(pheno)),
    castor::get_all_node_depths(tree_r)
  )
  max_depth <- max(tip_node_depths)
  # get edge table with normalized rates
  edge_tbl <- tibble::tibble(V1 = tree_r$edge[, 1], V2 = tree_r$edge[, 2]) |>
    dplyr::mutate(len = tree_r$edge.length) |>
    dplyr::mutate(top_node_depth = map_dbl(V1, ~ tip_node_depths[.x])) |>
    dplyr::mutate(btm_node_depth = map_dbl(V2, ~ tip_node_depths[.x])) |>
    dplyr::mutate(rootlen = sqrt(len)) |>
    dplyr::mutate(diff = node_tip_vals[V2] - node_tip_vals[V1]) |>
    dplyr::mutate(V3 = diff / rootlen)
  if (cut_by=="edge_length") {
    edge_tbl <- dplyr::mutate(edge_tbl, cut_by = rootlen)
  } else if (cut_by=="parent_depth") {
    edge_tbl <- dplyr::mutate(edge_tbl, cut_by = top_node_depth)
  } else {
    edge_tbl <- dplyr::mutate(edge_tbl, cut_by = btm_node_depth)
  }
  for (i in 1:n) {
    # permute rates within bins of equal size based on length, parent node depth, or child node depth
    binned_edge_tbl <- edge_tbl |>
      dplyr::mutate(binned_lengths = cut(
        top_node_depth,
        quantile(top_node_depth, seq(0, 1, 1/nbins))
      )) |>
      dplyr::group_by(binned_lengths) |>
      dplyr::mutate(scrambled_V3 = sample(V3)) |>
      dplyr::mutate(scrambled_rescale = scrambled_V3 * rootlen)
    new_nt_vals <- node_tip_vals
    # since these are already sorted in depth-first root to tip traversal order:
    for (x in 1:nrow(binned_edge_tbl)) {
      new_nt_vals[binned_edge_tbl[["V2"]][x]] <-
        new_nt_vals[binned_edge_tbl[["V1"]][x]] +
        binned_edge_tbl[["scrambled_rescale"]][x]
    }
    perm_mtx[, i] <- new_nt_vals[tree_r$tip.label]
  }
  if (rank) {
    new_pheno_vals <- sort(pheno)[order(perm_mtx[,i])]
    perm_mtx[, i] <- new_pheno_vals
  }
  rownames(perm_mtx) <- tree_r$tip.label
  perm_mtx
}


#' Helper to z-score a vector
#' @export zsc
zsc <- function(x) (x - mean(na.omit(x))) / sd(na.omit(x))

#' Helper to estimate p-values using an empirical Z test
#' @export normal_estimate_pval
normal_estimate_pval <- function(bg_stats, fg_stat) {
  center <- mean(bg_stats)
  sd <- sd(bg_stats)
  2 * (1 - pnorm(abs(fg_stat - center)/sd))
}

#' Run robust test given a real phenotype, permulated/permutrated phenotypes, and genes (all in PIC space)
#' @param realphenoPIC Real PIC vector for a phenotype.
#' @param permPICs A matrix of PIC vectors for permulated phenotypes.
#' @param genePIC PIC vector for a gene (predictor).
#' @param early_a The `a` parameter for early stopping. Controls how fast the threshold decays before it is possible to stop early.
#' @param early_p The `p` parameter for early stopping. All p-values below this threshold should be reliable.
#' @param mean_center_nulls Should the null distribution be mean-centered?
#' @param add_pseudocount Should a pseudocount be added to the top and bottom of the empirical p-value fraction?
#' @param approx_pvals Should the p-values be approximated using a continuous function? (Default: TRUE)
#' @param approx_method Can be either "normal" or "density." Controls how p-values are estimated if `approx_pvals==TRUE`.
#' @param regression_method Controls the function used to fit the data (either "rlm" or "lm" are recognized).
#' @param return_bgs Return the full set of null statistics?
#' @export repermulize_test
repermulize_test <- function(
  realphenoPIC,
  permPICs,
  genePIC,
  early_a = 10,
  early_p = 0.1,
  mean_center_nulls = FALSE,
  add_pseudocount = FALSE,
  approx_pvals = TRUE,
  approx_method = "normal",
  regression_method = "rlm",
  return_bgs=FALSE
) {
  early_c <- early_p / (1 + early_p)
  if (regression_method=="rlm") regress <- MASS::rlm else regress <- stats::lm
  tryCatch({
    model_fit <- regress(realphenoPIC ~ genePIC - 1)
    model_coef <- summary(model_fit)$coefficients
  }, error = \(e) {
    warning(paste0(e))
    return(c(Estimate = NA, p.value = NA))
  })
  model_fit <- regress(realphenoPIC ~ genePIC - 1) # nb, leave out intercept in PIC
  model_coef <- summary(model_fit)$coefficients
  fg_tv <- model_coef[1, "t value"]
  fg_est <- model_coef[1, 1]
  fg_df <- summary(model_fit)$df
  bg_tvs <- numeric(ncol(permPICs)) # empty vector to start
  permPICs <- permPICs[, sample(1:ncol(permPICs))] # "shuffle the deck" since we are using early stopping
  nperms = ncol(permPICs)
  for (ix in 1:nperms) {
    tryCatch(
      bg_tvs[ix] <- summary(regress(permPICs[, ix] ~ genePIC - 1))$coefficients[
        1,
        "t value"
      ],
      error = \(e) {
        warning(paste0(e))
        bg_tvs[ix] <- NA
      })
    n_actual_tests <- sum(!is.na(bg_tvs[1:ix]))
    # Use this if the null distribution is skewed away from 0 towards positives or negatives,
    # but not before there is a meaningful mean
    if (mean_center_nulls && (n_actual_tests >= 5)) {
      null_tvs <- abs(bg_tvs[1:ix] - mean(na.omit(bg_tvs[1:ix])))
      test_tv <- abs(fg_tv - mean(na.omit(bg_tvs[1:ix])))
    } else {
      null_tvs <- abs(na.omit(bg_tvs[1:ix]))
      test_tv <- abs(fg_tv)
    }
    running_pv <- mean(na.omit(null_tvs) >= test_tv)
    if (running_pv > ((early_a / n_actual_tests) + early_c) / (1 + early_c)) {
      break # early stopping rule
    }
  }
  if (approx_pvals) {
    if (approx_method == "density") true_pv <- density_estimate_pval(na.omit(bg_tvs[1:ix]), fg_tv)
    if (approx_method == "normal") true_pv <- normal_estimate_pval(na.omit(bg_tvs[1:ix]), fg_tv)
  } else {
    # exact 0 p-values can be problematic, so we can optionally add the equivalent of one pseudocount
    if (add_pseudocount) {
      true_pv <- ((running_pv * nperms) + 1 + 0) / (nperms + 2)
    } else {
      true_pv <- running_pv
    }
  }
  to_return <- c(Estimate = fg_est, p.value = true_pv)
  if (return_bgs) attr(to_return, "bg_vals") <- bg_tvs[1:ix]
  if (return_bgs) attr(to_return, "fg_val") <- fg_tv
  to_return
}

#' Helper to put PICs in node order, leaving NAs if any missing (shouldn't be)
#' @export order_pic_wrapper
order_pic_wrapper <- function(a, b) {
  ic <- castor::get_independent_contrasts(a, b)
  p <- rep(NA, a$Nnode)
  for (i in 1:a$Nnode) {
    p[ic$nodes[i]] <- ic$PICs[i]
  }
  p
}

#' Wrap the entire thing including making PICs and running genes in parallel
#' 
#' Run robust test given a real phenotype, permulated/permutrated phenotypes, and genes (all in PIC space)
#' @param real_pheno Real phenotype vector (named values).
#' @param real_genes Real matrix of gene presence/absence (or pangenome prevalence).
#' @param real_tree Tree to use for phylogenetic regressions.
#' @param perm_pheno Optionally, can directly provide permulated or otherwise randomized phenotypes.
#' @param real_pheno_sd Optionally, a vector giving the standard deviation of the mean for each phenotype value.
#' @param save_everything Should intermediate results be returned? Useful for debugging. (Default: FALSE) 
#' @param p_model Optionally, a model of evolution other than Brownian motion (used to transform the tree; see `phylolm::transf.branch.lengths()`).
#' @param ub Upper bound on parameter in p_model.
#' @param lb Lower bound on parameter in p_model.
#' @param verbose Talk more?
#' @param perm_method Should be "permulate" or "permutrate".
#' @param rank For permutration specifically, should the values be replaced with the original values based on their rank, like in permulation? (Default: FALSE)
#' @export repermulize_wrapper
repermulize_wrapper <- function(
  real_pheno,
  real_genes,
  real_tree,
  perm_pheno = NULL,
  real_pheno_sd = NULL,
  save_everything = FALSE,
  genes_are_PICs = FALSE,
  real_PICs = NULL,
  just_scramble_PICs = FALSE,
  perm_model = "BM",
  n = 1000,
  ub = 0.99,
  lb = 0.01,
  verbose = TRUE,
  perm_method = "permulate",
  rank = FALSE,
  ...
) {
  # make sure same species represented in tree, genes, and phenotype
  tips <- intersect(real_tree$tip.label,
    intersect(colnames(real_genes), names(real_pheno)))
  if (length(tips) < 3) {
    stop(
      "Too few tips remaining after matching gene matrix, phenotype, and tree"
    )
  }
  real_genes <- real_genes[, tips]
  reduced_tree <- ape::keep.tip(real_tree, tips)
  real_pheno <- real_pheno[reduced_tree$tip.label] # put in same order as tree
  if (!is.null(real_pheno_sd)) {
    real_pheno_sd <- real_pheno_sd[reduced_tree$tip.label]
  }

  if (verbose) pz.message("Calculating real PICs...")
  # calculate real PICs
  if (!is.null(real_pheno_sd)) {
    # use the actual SDs to best augment the tree
    real_rescale <- add_uncertainty_to_tree(
      reduced_tree,
      real_pheno,
      real_pheno_sd
    )
    real_rescaled_tree <- real_rescale$aug_tree
  } else {
    # still allow using an alternative model to capture, e.g., general measurement error
    real_rescaled_tree <- get_best_model_fit(
      real_pheno,
      reduced_tree,
      perm_model,
      yield = "tree",
      ub = ub,
      lb = lb
    )
  }
  if (is.null(real_PICs)) {
      real_PICs <- order_pic_wrapper(real_rescaled_tree, real_pheno)
  }
  if (verbose) pz.message("Calculating fake PICs...")
  if (is.null(real_pheno_sd)) perm_pheno_sd = NULL
  if (is.null(perm_pheno)) {
    if (just_scramble_PICs) {
      perm_pheno <- NULL
      perm_PICs <- Reduce(cbind, purrr::map(1:n, ~ sample(real_PICs)))
    } else {
      if (perm_method=="permulate") {
        perm_pheno <- permulate(
          real_pheno,
          real_rescaled_tree,
          n,
          pheno_sd = NULL, # don't re-fit based on the standard deviation
          p_model = "BM" # don't re-scale the tree since we're already doing that
        )
        
      } else if (perm_method=="permutrate") {
        perm_pheno <- permutrate(real_pheno, real_rescaled_tree, n, rank=rank)
      } else {
        stop(sprintf("unknown perm_method: %s", perm_method))
      }
      perm_pheno <- perm_pheno[reduced_tree$tip.label, ]
      # Don't use parallel for this -- doesn't seem to be worth the overhead
      perm_PICs <- pbapply::pbapply(perm_pheno, 2, \(x) {
					    order_pic_wrapper(real_rescaled_tree, x)
      })
    }
  }
  
  # loop across genes, using futures to do parallel computing without requiring lots of memory
  # note, we will calculate PICs inside the gene loop to avoid casting to a dense matrix
  named_indices <- 1:nrow(real_genes)
  names(named_indices) <- rownames(real_genes)

  if (verbose) message("Getting empirical p-values...")
  pbapply::pboptions(type = "timer")
  res <- pbapply::pblapply(
    named_indices,
    \(i) {
      if (genes_are_PICs) {
        g_pic <- real_genes[i, ]
      } else {
        g_pic <- order_pic_wrapper(real_rescaled_tree, 1 * real_genes[i, ])
      }
      tryCatch(
        repermulize_test(real_PICs, perm_PICs, g_pic, ...),
        error = function(e) { warning(paste(e)); c(Estimate=NA,p.value=NA) }
      )
    },
    cl = "future",
    future.seed=TRUE
  )

  # Collect and convert to a more familiar format
  r_df <- res |> as.data.frame()
  
  # Return PICs, etc. if desired (can't return the gene PICs since those are inside the parallel loop)
  if (save_everything) {
    return(list(
      r_df = r_df,
      perm_pheno = perm_pheno,
      perm_PICs = perm_PICs,
      real_PICs = real_PICs
    ))
  } else {
    return(r_df)
  }
}

#' Helper to get best fit to a more complex model and either return model parameters or a rescaled tree
#' @export get_best_model_fit
get_best_model_fit <- function(pheno, tree, p_model="lambda", ub=0.99, lb=0.01, yield="tree") {
  model_fit <- phylolm::phylolm(
    pheno ~ 1,
    phy = tree,
    model = p_model,
    upper.bound = ub,
    lower.bound = lb
  )
  optpar <- model_fit$optpar
  rbt_params <- list(sigma2 = model_fit$sigma2)
  optparam_names <- c(
    lambda = "lambda",
    delta = "delta",
    kappa = "kappa",
    EB = "rate"
  )
  if (p_model %in% names(optparam_names)) {
    rbt_params[optparam_names[p_model]] <- optpar
  }
  if (yield == "tree") {
    phylolm::transf.branch.lengths(tree, p_model, parameters = rbt_params)$tree
  } else if (yield == "params") {
    return(rbt_params)
  } else if (yield == "both") {
    return(list(
      tree=phylolm::transf.branch.lengths(tree, p_model, parameters = rbt_params)$tree,
      params=rbt_params
    ))
  } else {
    warning(paste0("unrecognized value: ", yield))
    return(rbt_params)
  }
}

#' Helper to estimate p-values using a Gaussian kernel density estimator instead of using the empirical distribution
#' @export density_estimate_pval
density_estimate_pval <- function(bg_stats, fg_stat, bw=bw.SJ, ...) {
  bg_bw <- bw(bg_stats)
  tail_prob <- mean(pnorm(fg_stat, mean = bg_stats, sd = bg_bw))
  return(min(2 * tail_prob, 2 * (1-tail_prob)))
}

#' Function to transform a tree by lengthening the tips to match measurement error for each individual taxon.
#' @param phy Tree object to transform.
#' @param pheno Named phenotype vector (numeric).
#' @param sd Named vector of standard deviations around the means represented in `pheno`.
#' @param lb In log10-space, the lower bound to search for the scale parameter.
#' @param ub In log10-space, the upper bound to search for the scale parameter.
#' @export add_uncertainty_to_tree
add_uncertainty_to_tree <- function(phy, pheno, pheno_sd, lb = -3, ub = 3) {
  # make sure all aligned
  common_taxa <- na.omit(intersect(
    phy$tip.label,
    intersect(names(pheno), names(pheno_sd))
  ))
  phy <- ape::keep.tip(phy, common_taxa)
  pheno <- pheno[phy$tip.label]
  pheno_sd <- pheno_sd[phy$tip.label]

  pheno_var_std <- (pheno_sd) / mean(pheno_sd)

  best_scale <- optimize(
    f = function(s) {
      AIC_scale_and_augment(phy, pheno, pheno_var_std, s)
    },
    interval = 10**seq(lb, ub, 0.01)
  )
  if (abs(best_scale$minimum - lb) <= 1e-5) {
    warning("Scale parameter was close to the lower bound")
  }
  if (abs(best_scale$minimum - ub) <= 1e-5) {
    warning("Scale parameter was close to the upper bound")
  }
  aug_tree <- augment_tip_branch_lengths(
    phy,
    pheno_var_std * best_scale$minimum
  )

  list(best_scale = best_scale, aug_tree = aug_tree)

}

#' Helper function to compute best AIC for a given scaling of per-tip stdevs
#' @export AIC_scale_and_augment
AIC_scale_and_augment <- function(phy, pheno, v, scale, meas_err=FALSE) {
  phy_aug <- augment_tip_branch_lengths(phy, v * scale)
  AIC(phylolm::phylolm(pheno ~ 1, phy = phy_aug, lower.bound=0, measurement_error = meas_err))
}

#' Helper function to add a vector of "extra" tip-adjacent branch lengths to a tree
#' @export augment_tip_branch_lengths
augment_tip_branch_lengths <- function(phy, tip_amounts) {
  tip_amounts <- tip_amounts[phy$tip.label]
  edge_orders <- purrr::map_dbl(1:length(phy$tip.label),
    \(i) which(phy$edge[, 2] == i))
  phy$edge.length[edge_orders] <- phy$edge.length[edge_orders] + tip_amounts
  phy
}
