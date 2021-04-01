fetch_DM_beta <- function(
  DMR_obj,
  beta_stats,
  probe_coords,
  tissue_key
) {

  library(reshape)
  library(dplyr)
  library(naturalsort)

  # fetch normal beta stats:
  for (j in 1:length(beta_stats)) {

    # isolate probes DM in all non malig:
    DM_stats <- as.data.frame(
      beta_stats[[j]][
        rownames(beta_stats[[j]]) %in% 
        DMR_obj$probe,
      ]
    )
    DM_stats$probe <- rownames(DM_stats)
    DM_stats$tissue <- names(beta_stats)[j]

    # return as df:
    if (j==1) {
      DM_stats_list <- list(DM_stats)
    } else {
      DM_stats_list[[j]] <- DM_stats
    }
      
  }
  names(DM_stats_list) <- names(beta_stats)

  # merge to dataframe:
  if (nrow(DM_stats_list[[1]]) > 0) {

  	DM_stats <- do.call("rbind", DM_stats_list)
  
    # fetch probe locations and add:
    specific_coords <- probe_coords[
      names(probe_coords) %in% DM_stats$probe
    ]
    DM_stats <- merge(
      DM_stats,
      data.frame(
        probe = names(specific_coords),
        chr = seqnames(specific_coords),
        coord = start(specific_coords)
      ),
      by = "probe"
    )

    # remove all_non_malig:
    DM_stats <- DM_stats[
      DM_stats$tissue != "all_non_malig",
    ]

    # substitute tissue codes for names:
    m <- match(DM_stats$tissue, tissue_key$code)
    DM_stats$tissue <- tissue_key$name[m]

    # order by chromosome, then position:
    DM_stats$chr <- factor(
      DM_stats$chr,
      levels = naturalsort(
        unique(
          as.character(DM_stats$chr)
        )
      )
    )

    # add NF scores:
    colnames(DM_stats)[3] <- "quart_0.1"
    colnames(DM_stats)[4] <- "quart_0.9"
    DM_stats <- arrange(DM_stats, chr, coord)
    temp_df <- DM_stats[!duplicated(DM_stats$probe),]
    DM_stats <- rbind(
      DM_stats,
      data.frame(
        probe = DMR_obj$probe,
        mean = DMR_obj$score,
        quart_0.1 = NA,
        quart_0.9 = NA,
        sd = NA,
        tissue = rep("NF", nrow(temp_df)),
        chr = temp_df$chr,
        coord = temp_df$coord
      )
    )
    DM_stats <- arrange(DM_stats, chr, coord)

    # add column for each x position:
    DM_stats$x_pos <- factor(DM_stats$probe)
    levels(DM_stats$x_pos) <- 1:length(unique(DM_stats$x_pos))

    # add column to make NF points different shape:
    DM_stats$shape <- "non_malignant"
    DM_stats$shape[DM_stats$tissue == "NF"] <- "NF"
    DM_stats$shape[DM_stats$tissue == "MPNST"] <- "MPNST"
    DM_stats$shape[DM_stats$tissue == "blood"] <- "blood"

    # order shapes:
    DM_stats$shape <- factor(DM_stats$shape)
    DM_stats$shape <- relevel(DM_stats$shape, "non_malignant")

    return(DM_stats)

  } else {
  	return(NULL)
  }

}