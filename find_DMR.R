
# define parameters
min_sample_per_probe <- 5
min_beta_val <- 0.2
subset_mtx <- F

# fetch other variables:
args = commandArgs(trailingOnly=TRUE)

#beta_cutoff <- args[1]
#print(beta_cutoff)
#pval_cutoff <- as.character(args[2])
#background_cutoff <- args[3]
#blood_cutoffs <- list(
#  MPNST_hypermethylated = args[4],
#  MPNST_hypomethylated = args[5]
#)

beta_cutoff <- 0.35
pval_cutoff <- as.character("1e-4")
background_cutoff <- 0.2
blood_cutoffs <- list(
  MPNST_hypermethylated = 0.2,
  MPNST_hypomethylated = 0.8
)

# define and create directories:
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/MPNST_ctMet/")
raw_dir <- paste0(project_dir, "raw_files/")
script_dir <- paste0(project_dir, "scripts/")
func_dir <- paste0(script_dir, "functions/")

normal_in <- paste0(raw_dir, "TCGA_normal/")
Robject_in <- paste0(raw_dir, "Rdata/")
result_dir <- paste0(project_dir, "results/")
out_dir <- paste0(
  result_dir, "/",
  "beta_", beta_cutoff, "/",
  "pval_", pval_cutoff, "/",
  "background_", background_cutoff, "/",
  "blood_hypo_", blood_cutoffs$MPNST_hypermethylated, "/",
  "blood_hyper_", blood_cutoffs$MPNST_hypomethylated, "/"
)

if (subset_mtx) {
  result_dir <- paste0(result_dir, "subset/")
  out_dir <- paste0(out_dir, "subset/")
}

Robject_dir <- paste0(result_dir, "Rdata/")
plot_dir <- paste0(out_dir, "plots/")
table_dir <- paste0(out_dir, "tables/")

system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", plot_dir))
system(paste0("mkdir -p ", table_dir))

print(paste0("Output candidate counts in ", table_dir, "DMR_counts.txt"))


####################################################################################
### 0. Load packages and functions ###
####################################################################################

library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)

DMR_analysis <- dget(paste0(func_dir, "DMR_analysis.R"))
DMR_filter <- dget(paste0(func_dir, "DMR_filter.R"))
filter_background <- dget(paste0(func_dir, "filter_background.R"))


####################################################################################
### 1. Identify MPNST vs normal DMRs ###
####################################################################################

print("Loading non-malignant beta values...")

if (!file.exists(paste0(Robject_dir, "normal_se_object.Rdata"))) {

  # load normal se objects:
  normal_se_files <- list.files(normal_in, pattern = "rda")
  
  for (i in 1:length(normal_se_files)) {
  
    print(paste0("Loading ", normal_se_files[i], "..."))
  
    load(paste0(normal_in, normal_se_files[i]))
  
    if (i==1) {
      normal_se <- list(data)
    } else {
      normal_se[[i]] <- data
    }
  
    names(normal_se)[i] <- gsub(".rda", "", normal_se_files[i])
  
  }
  
  saveRDS(normal_se, paste0(Robject_dir, "normal_se_object.Rdata"))

} else {

  normal_se <- readRDS(paste0(Robject_dir, "normal_se_object.Rdata"))

}

print("Formatting non-malignant beta values...")

# isolate beta matrices from se objects:
normal_mtx <- lapply(normal_se, function(x) assays(x)[[1]])

# merge those with less than 5 samples:
if (exists("mixed_mtx")) {
  rm(mixed_mtx)
}

for (i in 1:length(normal_mtx)) {

  if (ncol(normal_mtx[[i]]) < 5) {

  	if (!exists("mixed_mtx")) {
  	  mixed_mtx <- normal_mtx[[i]]
  	  to_rm <- names(normal_mtx)[i]
  	} else {
  	  mixed_mtx <- cbind(mixed_mtx, normal_mtx[[i]])
  	  to_rm <- c(to_rm, names(normal_mtx)[i])
  	}

  }
}
normal_mtx <- normal_mtx[!(names(normal_mtx) %in% to_rm)]
normal_mtx$mixed <- mixed_mtx

# add all normal matrix:
normal_mtx$all <- do.call("cbind", normal_mtx)

# rename elements:
names(normal_mtx) <- paste0(
  gsub(
    "_.*$", "", 
    gsub(
    	"^.*-", "", 
    	names(normal_mtx)
    )
  ),
  "_non_malig"
)

print("Loading MPNST beta values...")

# Load MPNST beta values:
MPNST_beta_df <- readRDS(paste0(Robject_in, "MPNST_beta_df.Rdata"))

# label columns:
colnames(MPNST_beta_df) <- paste0(colnames(MPNST_beta_df), "-MPNST")

if (subset_mtx) {

  normal_mtx <- lapply(normal_mtx, function(x) x[1:100,])

  MPNST_beta_df <- MPNST_beta_df[1:100,]

}

# find DMRs between each group and MPNST groups:
if (!file.exists(paste0(Robject_dir, "MPNST_vs_normal_groups.Rdata"))) {
	
  for (i in 1:length(normal_mtx)) {

  	print(
  	  paste0(
  	  	"Identifying DM probes, MPNST vs ", names(normal_mtx)[i], "..."
  	  )
  	)

    DMR <- DMR_analysis(
      sample1_beta = MPNST_beta_df, 
      sample1_name = "MPNST", 
      sample2_beta = normal_mtx[[i]], 
      sample2_name = names(normal_mtx)[i],
      row_ranges = rowRanges(normal_se[[1]]),
      plot_dir
    )
  
    if (i==1) {
    	MPNST_vs_normal <- list(DMR)
    } else {
    	MPNST_vs_normal[[i]] <- DMR
    }
  
  }
  names(MPNST_vs_normal) <- names(normal_mtx)

  saveRDS(
  	MPNST_vs_normal, 
  	paste0(Robject_dir, "MPNST_vs_normal_groups.Rdata")
  )

} else {

  MPNST_vs_normal <- readRDS(
  	paste0(Robject_dir, "MPNST_vs_normal_groups.Rdata")
  )

}

# record values:
unfiltered_hyper <- unlist(
  lapply(MPNST_vs_normal, function(x) {
    nrow(
    	x[x$status == "Hypermethylated in MPNST",]
    )
  })
)

unfiltered_hypo <- unlist(
  lapply(MPNST_vs_normal, function(x) {
    nrow(
    	x[x$status == "Hypomethylated in MPNST",]
    )
  })
)

DMR_record <- rbind(unfiltered_hypo, unfiltered_hyper)


####################################################################################
### 2. Filter DM probes ###
####################################################################################

print(
  paste0(
    "Filtering for canditates with at least ", 
    beta_cutoff, 
    "beta value difference and < ", pval_cutoff, "p-value..."
  )
)

DMR_initial <- lapply(
  MPNST_vs_normal,
  DMR_filter,
  beta_cutoff,
  as.numeric(pval_cutoff)
)

# record values:
beta_filt_hyper <- unlist(
  lapply(DMR_initial, function(x) {
    nrow(
      x$hyper
    )
  })
)

beta_filt_hypo <- unlist(
  lapply(DMR_initial, function(x) {
    nrow(
      x$hypo
    )
  })
)

DMR_record <- rbind(
  DMR_record,
  rbind(beta_filt_hyper, beta_filt_hypo)
)

print(
  paste0(
  	beta_filt_hyper[16], " hypermethylated and ", beta_filt_hypo[16],
  	" hypomethylated beta value/p-value filtered candiates identified..."
  )
)

print(
  paste0(
    "Filtering for canditates hypermethylated in MPNST with normal ", 
    "background level values < ", 
    background_cutoff, 
    " and hypomethylated in MPNST with MPNST background values > ", 
    background_cutoff, "..."
  )
)


for (i in 1:length(DMR_initial)) {

  initial_filt <- lapply(DMR_initial[[i]], function(x) {

  	if (length(grep("Hypermethylated", x$status[1])) > 0) {

  	  return(
  	  	filter_background(
          dmr = x,
          background_cutoff,
          hypo_beta = normal_mtx[[i]],
          min_sample_per_probe
        )
      )

  	} else {

  	  return(
    	filter_background(
          dmr = x,
          background_cutoff,
          hypo_beta = MPNST_beta_df,
          min_sample_per_probe
        )
      )

  	}
  	
  })

  if (i==1) {
  	DMR_bg_filt <- list(initial_filt)
  } else {
  	DMR_bg_filt[[i]] <- initial_filt
  }

}
names(DMR_bg_filt) <- names(DMR_initial)

# record values:
bg_filt_hyper <- unlist(
  lapply(DMR_bg_filt, function(x) {
    nrow(
      x$hyper
    )
  })
)

bg_filt_hypo <- unlist(
  lapply(DMR_bg_filt, function(x) {
    nrow(
      x$hypo
    )
  })
)

DMR_record <- rbind(
  DMR_record,
  rbind(bg_filt_hyper, bg_filt_hypo)
)

print(
  paste0(
  	bg_filt_hyper[16], " hypermethylated and ", bg_filt_hypo[16],
  	" hypomethylated high background filtered candiates identified..."
  )
)


####################################################################################
### 3. Remove hypo/hypermethylated blood regions ###
####################################################################################

print("Loading blood beta values...")

# load GEO blood data and fetch beta matrix:
blood_se <- readRDS(paste0(raw_dir, "GEO_blood/SE_GEO_Meth_Data.rds"))
blood_mtx <- assays(blood_se)[[1]][-1,]

print(
  paste0(
    "Filtering for canditates hypermethylated in MPNST with blood values < ", 
    blood_cutoffs[1], 
    " and hypomethylated in MPNST with blood values > ", 
    blood_cutoffs[2], "..."
  )
)

# identify and remove hypermethylated probes with median beta values > 
# blood_cutoffs$MPNST_hypermethylated in blood, and hypomethylated probes with 
# median beta values < blood_cutoffs$MPNST_hypomethylated in blood:
DMR_blood_filt <- lapply(DMR_bg_filt, function(x) {

	res <- lapply(x, function(y) {

	  # fetch blood values for DM probes and ensure data for enough samples
	  # for each probe:
	  blood_beta <- blood_mtx[rownames(y),]

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
  	  if (length(grep("Hypermethylated", y$status[1])) > 0) {

  	  	filt_blood <- names(blood_medians)[
  	  	  blood_medians < blood_cutoffs$MPNST_hypermethylated
  	  	]

  	  } else {

  	  	filt_blood <- names(blood_medians)[
  	  	  blood_medians > blood_cutoffs$MPNST_hypomethylated
  	  	]

  	  }

  	  return(y[filt_blood,])

	})
	
	return(res)
  
})

# record values:
blood_filt_hyper <- unlist(
  lapply(DMR_blood_filt, function(x) {
    nrow(
      x$hyper
    )
  })
)

blood_filt_hypo <- unlist(
  lapply(DMR_blood_filt, function(x) {
    nrow(
      x$hypo
    )
  })
)

DMR_record <- rbind(
  DMR_record,
  rbind(blood_filt_hyper, blood_filt_hypo)
)

print(
  paste0(
  	blood_filt_hyper[16], " hypermethylated and ", blood_filt_hypo[16],
  	" hypomethylated blood filtered candiates identified..."
  )
)


####################################################################################
### 4. Filter for probes DM in MPNST vs each normal ###
####################################################################################

# separate hyper and hypo dfs:
DMR_list <- list(
  hyper = lapply(DMR_blood_filt, function(x) x$hyper),
  hypo = lapply(DMR_blood_filt, function(x) x$hypo)
)

# find DMR probes common to all normal tissue:
DMR_probes <- lapply(DMR_list, function(x) {
  lapply(x, rownames)
})
DMR_final_probes <- lapply(DMR_probes, function(x) {

  if (exists("common_probes")) {
    rm(common_probes)
  }

  for (i in 2:length(x)) {

  	if (!exists("common_probes")) {
  	  common_probes <- intersect(x[[i-1]], x[[i]])
  	} else {
  	  common_probes <- intersect(common_probes, x[[i]])
  	}

  }
  return(common_probes)

})

# fetch final DM values:
if (exists("DMR_final")) {
  rm(DMR_final)
}
for (i in 1:length(DMR_list)) {

  if (length(DMR_list[[i]]) > 0) {

  	if (!exists("DMR_final")) {

  	  DMR_final <- list(
  	  	lapply(DMR_list[[i]], function(y) {
          y[rownames(y) %in% DMR_final_probes[[i]],]
        })
      )
      names(DMR_final) <- names(DMR_list)[i]

  	} else {

  	  DMR_final[[i]] <- lapply(DMR_list[[i]], function(y) {
        y[rownames(y) %in% DMR_final_probes[[i]],]
      })
      names(DMR_final)[i] <- names(DMR_list)[i]

  	}
  	
  }
  
}

# record values:
final_hyper <- unlist(
  lapply(DMR_final$hyper, function(x) {
    nrow(x)
  })
)

final_hypo <- unlist(
  lapply(DMR_final$hypo, function(x) {
    nrow(x)
  })
)

DMR_record <- rbind(
  DMR_record,
  rbind(final_hyper, final_hypo)
)

print(
  paste0(
  	final_hyper[16], " hypermethylated and ", final_hypo[16],
  	" hypomethylated final candiates identified!"
  )
)

# save as table:
write.table(
  DMR_record,
  paste0(table_dir, "DMR_counts.txt"),
  row.names = T,
  col.names = T,
  quote = F
)


####################################################################################
### 6. Plot candidates ###
####################################################################################
# write function to fetch mean beta values of all tissues:
fetch_DM_beta <- function(
  DMR_obj
) {

  # fetch normal mean beta vals:
  mean_beta_list <- lapply(DMR_obj, function(y) {

    y_sub <- y[,grep("mean.*non.malig", colnames(y))]
    return(
      data.frame(
        row.names = rownames(y),
        non_malig = y_sub[, grep("minus", colnames(y_sub), invert = T)]
      )
    )

  })

  # fetch MPNST mean beta vals:
  mean_beta_list$MPNST <- data.frame(
  	row.names = rownames(DMR_obj[[1]]),
  	MPNST = DMR_obj[[1]][
      , colnames(DMR_obj[[1]]) == "mean.MPNST"
    ]
  )

  # remove elements with no DM probes:
  mean_beta_list <- mean_beta_list[lapply(mean_beta_list, nrow) > 0]

  # merge to dataframe:
  if (length(mean_beta_list) > 0) {
  	mean_beta <- do.call("cbind", mean_beta_list)
    colnames(mean_beta) <- names(mean_beta_list)
  
    # add probe column and melt for plotting:
    mean_beta$probe <- rownames(mean_beta)
    return(melt(mean_beta))
  } else {
  	return(NULL)
  }

}

# fetch mean beta values of all tissues:
DMR_plot_dfs <- lapply(DMR_final, fetch_DM_beta)

# create line plot:
p <- ggplot(DMR_plot_dfs$hypo, aes(x=probe, y=value, group=variable, color=variable))
p <- p + geom_line()

ggsave(
  paste0(plot_dir, "MPNST_hypomethylated_marker_probes.pdf"),
  plot = p,
  useDingbats = F
)


####################################################################################
### 5. Identify stringent DMRs ###
####################################################################################

# filter matrices for common probes:
normal_mtx <- lapply(normal_mtx, function(x) {
  return(x[rownames(MPNST_beta_df), ])
})

MPNST_beta_df <- MPNST_beta_df[rownames(normal_mtx[[1]]), ]

# cbind MPNST and normal beta values:
all_beta <- c(
  normal_mtx,
  list(MPNST_beta_df)
)
names(all_beta)[length(all_beta)] <- "MPNST"

# identify methylation probes with mean beta values < 0.05 and others > 0.8:
fringe_beta <- lapply(all_beta, function(x) { 

  mean_beta <- apply(x, 1, mean, na.rm = TRUE)

  low_beta <- names(mean_beta[mean_beta < 0.2])
  print(
  	paste0(
  	  "No. low beta vals = ", length(low_beta), " out of ",
  	  length(mean_beta)
  	)
  )

  high_beta <- names(mean_beta[mean_beta > 0.8])
  print(
  	paste0(
  	  "No. high beta vals = ", length(high_beta), " out of ",
  	  length(mean_beta)
  	)
  )

  return(
  	list(
  	  low_beta = low_beta,
  	  high_beta = high_beta
  	)
  )

})

# fetch list of probes with mean beta values < 0.05 in one tissue and > 0.8 
# in other tissue:
dm_ind <- list(
  MPNST_high = which(
    fringe_beta$MPNST$high_beta %in% fringe_beta$BRCA$low_beta
  ),
  MPNST_low = which(
    fringe_beta$MPNST$low_beta %in% fringe_beta$BRCA$high_beta
  )
)

dm_probes <- list(
  MPNST_high = fringe_beta$MPNST$high_beta[dm_ind$MPNST_high],
  MPNST_low = fringe_beta$MPNST$low_beta[dm_ind$MPNST_low]
)

# isolate only candidate regions of > 1 DM probes:
DMR <- lapply(dm_ind, function(x) {

  splt <- split(
    x, 
    cumsum(
    	seq_along(x) %in% (which(diff(x)>1)+1)
    )
  )

  dm_regions <- lapply(splt, function(x) {

    if (length(x) > 1) {
    	return(x)
    } else {
    	return(NULL)
    }

  })
  return(
  	dm_regions[
      !sapply(dm_regions, is.null)
    ]
  )

})

