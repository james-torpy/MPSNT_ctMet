blood_filter <- function(
    dmr, 
    blood_mtx,
    min_sample_per_probe
  ) {

  if (!is.na(dmr)) {
    # fetch blood values for DM probes and ensure data for enough samples
    # for each probe:
    blood_beta <- blood_mtx[rownames(dmr),]
  
    blood_beta <- blood_beta[
        rowSums(!is.na(blood_beta)) > min_sample_per_probe,
    ]

    # calculate blood median values:
    blood_medians <- apply(
        blood_beta, 
        1, 
        median,
        na.rm = TRUE
      )

      # filter out values above/below blood_cutoffs:
      if (length(grep("Hypermethylated", dmr$status[1])) > 0) {
        filt_blood <- names(blood_medians)[
          blood_medians < blood_cutoffs$MPNST_hypermethylated
        ]
      } else {
        filt_blood <- names(blood_medians)[
          blood_medians > blood_cutoffs$MPNST_hypomethylated
        ]
      }
      return(dmr[filt_blood,])
  } else {
    return(dmr)
  }
}