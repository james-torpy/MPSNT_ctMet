DMR_analysis <- function(
  sample1_beta,
  sample1_name,
  sample2_beta,
  sample2_name,
  row_ranges,
  min_sample_per_probe,
  plot_dir
) {

  print(
    paste0(
      "Detecting differentially methylated probes between ",
      sample1_name, " and ", sample2_name, "..."
    )
  )

	# fetch beta dfs:
  unfilt_mtx <- list(
    sample1_beta = as.matrix(sample1_beta),
    sample2_beta = as.matrix(sample2_beta)
  )
  
  # remove rows with < min_sample_per_probe samples not NA:
  filt_mtx <- lapply(unfilt_mtx, function(x) {
    
    return(
      x[
        apply(x, 1, function(y) {
          length(which(!is.na(y))) >= min_sample_per_probe
        }),
      ]
    )
  
  })
  
  # keep only rows in both, and merge:
  filt_mtx$sample1_beta <- filt_mtx$sample1_beta[
    rownames(filt_mtx$sample1_beta) %in% rownames(filt_mtx$sample2_beta),
  ]
  filt_mtx$sample2_beta <- filt_mtx$sample2_beta[
    rownames(filt_mtx$sample2_beta) %in% rownames(filt_mtx$sample1_beta),
  ]
  beta_mtx <- do.call("cbind", filt_mtx)
  
  # adjust other objects for new se:
  row_ranges <- row_ranges[names(row_ranges) %in% rownames(beta_mtx)]
  
  final_se <- SummarizedExperiment(
    assays = list(
      counts = beta_mtx
    ),
    rowRanges = row_ranges,
    colData = DataFrame(
      tissue = c(
        rep(sample1_name, ncol(filt_mtx$sample1_beta)),
        rep(sample2_name, ncol(filt_mtx$sample2_beta))
      )
    )
  )
  
  res <- list(
    DMR = TCGAanalyze_DMC(
      data = final_se,
      groupCol = "tissue",
      group1 = sample1_name,
      group2 = sample2_name,
      plot.filename = paste0(
      	plot_dir, 
      	sample1_name, 
      	"_vs_", 
      	sample2_name,
      	"_DMR_volcano.png"
      )
    ),
    filtered_probes = nrow(beta_mtx)
  )

}