fetch_DM_beta <- function(
  DMR_obj,
  all_mean_beta,
  probe_coords
) {

  library(reshape)

  # fetch normal mean beta vals:
  for (j in 1:length(all_mean_beta)) {

    # isolate probes DM in all non malig:
    DM_beta <- all_mean_beta[[j]][
      names(all_mean_beta[[j]]) %in% 
      rownames(DMR_obj$all_non_malig)
    ]
    
    # return as df:
    if (j==1) {

      mean_beta_list <- list(
        data.frame(
          probe = names(DM_beta),
          beta = DM_beta,
          tissue = names(all_mean_beta)[j]
        )
      )

    } else {

      mean_beta_list[[j]] <- data.frame(
        probe = names(DM_beta),
        beta = DM_beta,
        tissue = names(all_mean_beta)[j]
      )

    }
      
  }
  names(mean_beta_list) <- names(all_mean_beta)

  # merge to dataframe:
  if (nrow(mean_beta_list[[1]]) > 0) {

  	mean_beta <- do.call("rbind", mean_beta_list)

    # format tissue names:
    mean_beta$tissue <- gsub(
      "TCGA-|_450K_hg19", "", 
      mean_beta$tissue
    )
  
    # fetch probe locations and add:
    specific_coords <- probe_coords[
      names(probe_coords) %in% mean_beta$probe
    ]
    mean_beta <- merge(
      mean_beta,
      data.frame(
        probe = names(specific_coords),
        chr = seqnames(specific_coords),
        coord = start(specific_coords)
      ),
      by = "probe"
    )

    # melt for plotting:
    return(mean_beta)
  } else {
  	return(NULL)
  }

}