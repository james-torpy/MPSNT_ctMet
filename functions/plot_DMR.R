### This function arranges label data and plots differentially 
# methylated candidates ###

plot_DMR <- function(
  plot_df,
  dot_vals = "median",
  bar_vals = "quartile",
  plot_type = "all",
  meth_type = "hyper",
  probe_annot = TRUE
) {

  # make tissue column factor and set MPNST as level 1:
  nm_df <- plot_df[
    !(plot_df$tissue %in% c("MPNST", "NF", "blood")),
  ]
  nm_df <- nm_df[order(nm_df$tissue),]

  # reorder so MPNST, NF and blood are top levels
  nm_blood <- rbind(
	plot_df[plot_df$tissue == "blood",],
	nm_df
  )
  nm_nf <- rbind(
	plot_df[plot_df$tissue == "NF",],
	nm_blood
  )
  plot_df <- rbind(
	plot_df[plot_df$tissue == "MPNST",],
	nm_nf
  )

  plot_df$tissue <- factor(plot_df$tissue,
  	levels = unique(plot_df$tissue))

  # add no samples to legend:
  legend_n_vals <- paste0(
  	levels(plot_df$tissue),
  	" (",
  	n_vals[levels(plot_df$tissue)],
  	")"
  )

  # set x positions and make character factor:
  plot_df <- arrange(
    plot_df, chr, coord
  )
  plot_df$x_pos <- as.character(
    rleid(plot_df$probe)
  )

  # set up temp df for labels:
  temp_df <- plot_df[
    !duplicated(plot_df$probe),
  ]

  # plot data:
  p <- ggplot(
    plot_df, 
    aes_string(
      x="x_pos", 
      y=dot_vals, 
      color=factor(
      	plot_df$tissue, labels = legend_n_vals
      ),
      group="shape"
    )
  )
  if (bar_vals == "quartile") {
  	p <- p + geom_pointrange(
    	aes(ymin = quart_0.25, ymax = quart_0.75, shape=shape), 
    	size = 0.5,
      position=position_jitter(width=0.3, seed = 754)
    )
  } else if (bar_vals == "wide_quartile") {
    p <- p + geom_pointrange(
    	aes(ymin = quart_0.1, ymax = quart_0.9, shape=shape), 
    	size = 0.5,
      position=position_jitter(width=0.3, seed = 754)
    )
  } else if (bar_vals == "sd") {
  	p <- p + geom_pointrange(
    	aes(ymin = mean-sd, ymax = mean+sd, shape=shape), 
    	size = 0.5,
      position=position_jitter(width=0.3, seed = 754)
    )
  }
  if (probe_annot) {
    p <- p + scale_x_discrete(
      labels = paste0(
        temp_df$probe, 
        "\n", 
        temp_df$chr,
        ":",
        temp_df$coord
      )
    )
  } else {
  	p <- p + scale_x_discrete(
  	  labels = c(
        "candidate 1", 
        seq(2, nrow(DMR_plot_df_top[[i]]), 1)
      )
    )
  }
  p <- p + labs(color = "Tissue")
  p <- p + guides(shape = FALSE)
  p <- p + xlab(
  	paste0(
  	  "CpG probes ", meth_type, "methylated ",
  	  "in MPNST"
  	)
  )
  p <- p + ylab(
  	paste0(
  	  gsub("m", "M", dot_vals), " methylation level"
  	)
  )
  p <- p + theme_cowplot(12)
  p <- p + theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
  p <- p + scale_color_manual(
  	values = col_pal[1:length(unique(plot_df$tissue)),1]
  )
  
  print(
    paste0(
      "Plotting candidates, filename = ", plot_dir, plot_type, 
      "_MPNST_", meth_type, "methylated_marker_probes_", dot_vals, 
      "_", bar_vals, "_bars_probes_annotated.png"
    )
  )

  # define plot names
  if (probe_annot) {
  	pprefix <- paste0(
      plot_dir, plot_type, "_MPNST_", meth_type, 
      "methylated_marker_probes_", dot_vals, "_", bar_vals, 
      "_bars_probes_annotated"
    )
  } else {
  	pprefix <- paste0(
      plot_dir, plot_type, "_MPNST_", meth_type, 
      "methylated_marker_probes_", dot_vals, "_", bar_vals, 
      "_bars"
    )
  }

  # write plots:
  png(
    paste0(pprefix, ".png"),
    res = 300,
    unit = "in",
    height = 7,
    width = 15
  )
    print(p)
  dev.off()

  ggsave(
    paste0(pprefix, ".pdf"),
    plot = p,
    useDingbats = F,
    width = 10
  )

}