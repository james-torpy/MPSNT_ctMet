fetch_DM_beta <- function(
  dmr,
  beta_val_list,
  probe_coords,
  tissue_key,
  n_vals,
  plot_min_sample_no
) {

  library(reshape)
  library(dplyr)
  library(naturalsort)

  # fetch beta values for DM probes:
  DM_combined <- do.call("rbind", dmr)
  DM_probes <- DM_combined$probe
  DM_beta <- lapply(beta_val_list, function(x) {
    return(x[rownames(x) %in% DM_probes,])
  })

  # calculate stats on beta values:
  for (j in 1:length(DM_beta)) {

    if (j==1) {

      beta_stats <- list(
        as.data.frame(
          t(
            apply(DM_beta[[j]], 1, function(y) {
              # calculate mean, 10% and 90% quartiles:
              quantiles <- quantile(
                y, c(0.1, 0.25, 0.75, 0.9), na.rm = TRUE
              )
              return(
                c(
                  tissue = names(DM_beta)[j],
                  median = median(y, na.rm = TRUE),
                  mean = mean(y, na.rm = TRUE),
                  sd = sd(y, na.rm = TRUE),
                  quart_0.1 = quantiles[1],
                  quart_0.25 = quantiles[2],
                  quart_0.75 = quantiles[3],
                  quart_0.9 = quantiles[4]
                )
              )
            })
          )
        )
      )

    } else {

      beta_stats[[j]] <- as.data.frame(
        t(
          apply(DM_beta[[j]], 1, function(y) {
            # calculate mean, 10% and 90% quartiles:
            quantiles <- quantile(
              y, c(0.1, 0.25, 0.75, 0.9), na.rm = TRUE
            )
            return(
              c(
                tissue = names(DM_beta)[j],
                median = median(y, na.rm = TRUE),
                mean = mean(y, na.rm = TRUE),
                sd = sd(y, na.rm = TRUE),
                quart_0.1 = quantiles[1],
                quart_0.25 = quantiles[2],
                quart_0.75 = quantiles[3],
                quart_0.9 = quantiles[4]
              )
            )
          })
        )
      )

    }
    # make probe column:
    beta_stats[[j]]$probe <- rownames(beta_stats[[j]])

  }
  names(beta_stats) <- names(DM_beta)

  # filter out tissues with < 10 samples from beta stats:
  rm_tissues <- names(n_vals)[as.numeric(n_vals) < plot_min_sample_no]
  beta_stats <- beta_stats[
    !(names(beta_stats) %in% rm_tissues)
  ]

  # merge to dataframe:
  if (nrow(beta_stats[[1]]) > 0) {

  	DM_stats <- do.call("rbind", beta_stats)
  
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

    # fix up column names:
    colnames(DM_stats) <- sub(
      '^([^.]+.[^.]+).*', 
      '\\1', 
      colnames(DM_stats)
    )
    # add NF scores:
    temp_df <- DM_stats[!duplicated(DM_stats$probe),]
    DM_stats <- rbind(
      DM_stats,
      data.frame(
        probe = temp_df$probe,
        tissue = rep("NF", nrow(temp_df)),
        median = DM_combined$score,
        mean = DM_combined$score,
        sd = NA,
        quart_0.1 = NA,
        quart_0.25 = NA,
        quart_0.75 = NA,
        quart_0.9 = NA,
        chr = temp_df$chr,
        coord = temp_df$coord
      )
    )
    DM_stats <- arrange(DM_stats, chr, coord)

    # add column to make NF points different shape:
    DM_stats$shape <- "non_malignant"
    DM_stats$shape[DM_stats$tissue == "NF"] <- "NF"
    DM_stats$shape[DM_stats$tissue == "MPNST"] <- "MPNST"
    DM_stats$shape[DM_stats$tissue == "blood"] <- "blood"

    # order shapes:
    DM_stats$shape <- factor(DM_stats$shape)
    DM_stats$shape <- relevel(DM_stats$shape, "non_malignant")

    # convert stats into numeric values:
    for (
      k in which(
        colnames(DM_stats) == "median"
      ):which(
        colnames(DM_stats) == "quart_0.9"
      )
    ) {
      DM_stats[ ,k] <- round(as.numeric(DM_stats[ ,k]), 3)
    }
    
    # resplit into hyper and hypo candidates:
    split_stats <- list(
      hyper = DM_stats[DM_stats$probe %in% dmr$hyper$probe,],
      hypo = DM_stats[DM_stats$probe %in% dmr$hypo$probe,]
    )

    return(split_stats)

  } else {
  	return(NULL)
  }

}