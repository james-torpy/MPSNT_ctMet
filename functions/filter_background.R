filter_background <- function(
  dmr,
  background_cutoff,
  hypo_beta,
  min_sample_per_probe
) {

  # fetch beta values of DM probes and remove those represented 
  # by < min_sample_per_probe:
  dm_beta <- as.data.frame(hypo_beta[rownames(dmr),])
  dm_beta <- dm_beta[
    rowSums(
      apply(dm_beta, 2, function(y) !is.na(y))
    ) > min_sample_per_probe,
  ]

  if (length(dm_beta) > 0) {
    # calculate median beta values of DM probes:
    median_dm_beta <- apply(dm_beta, 1, median, na.rm = TRUE)
  
    # remove probes with median above background cutoff:
    filt_dm <- names(median_dm_beta)[median_dm_beta <= background_cutoff]
  
    return(dmr[filt_dm,])
  } else {
    return(NA)
  }
  
}