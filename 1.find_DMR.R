
# define parameters
min_sample_per_probe <- 5
min_beta_val <- 0.2
subset_mtx <- F

# fetch other variables:
args = commandArgs(trailingOnly=TRUE)

beta_cutoff <- args[1]
print(beta_cutoff)
pval_cutoff <- as.character(args[2])
background_cutoff <- args[3]
blood_cutoffs <- list(
  MPNST_hypermethylated = args[4],
  MPNST_hypomethylated = args[5]
)
NF_cutoffs <- list(
  MPNST_hypermethylated = args[6],
  MPNST_hypomethylated = args[7]
)
check_indiv_tissue <- as.logical(args[8])
plot_min_sample_no <- as.numeric(args[9])

#beta_cutoff <- 0.35
#pval_cutoff <- as.character("1e-5")
#background_cutoff <- 0.2
#blood_cutoffs <- list(
#  MPNST_hypermethylated = 0.2,
#  MPNST_hypomethylated = 0.8
#)
#NF_cutoffs <- list(
#  MPNST_hypermethylated = 0.2,
#  MPNST_hypomethylated = 0.8
#)
#check_indiv_tissue <- as.logical("FALSE")
#plot_min_sample_no <- 10

print(paste0("beta_cutoff = ", beta_cutoff))
print(paste0("pval_cutoff = ", pval_cutoff))
print(paste0("background_cutoff = ", background_cutoff))
print(paste0("blood_cutoffs = ", blood_cutoffs$MPNST_hypermethylated, ", ", blood_cutoffs$MPNST_hypomethylated))
print(paste0("NF_cutoffs = ", NF_cutoffs$MPNST_hypermethylated, ", ", NF_cutoffs$MPNST_hypermethylated))

# define tissues:
tissue_key <- data.frame(
  code = c(
    "BLCA_non_malig", "blood", "BRCA_non_malig", 
    "CESC_non_malig", "CHOL_non_malig", "COAD_non_malig", 
    "ESCA_non_malig", "HNSC_non_malig", "KIRP_non_malig", 
    "LIHC_non_malig", "LUAD_non_malig", "LUSC_non_malig", 
    "mixed_non_malig", "MPNST", "NF", 
    "PAAD_non_malig", "PRAD_non_malig", "READ_non_malig", 
    "SARC_non_malig", "STAD_non_malig", "THCA_non_malig", 
    "THYM_non_malig"
  ),
  name = c(
    "bladder", "blood", "breast",
    "cervical", "bile duct", "colon",
    "esophagus", "head/neck", "kidney",
    "liver", "lung", "lung",
    "mixed", "MPNST", "NF", 
    "pancreas", "prostate", "rectum", 
    "sarcoma", "stomach", "thyroid", 
    "thymus"
  )
)

# define and create directories:
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/MPNST_ctMet/")
raw_dir <- paste0(project_dir, "raw_files/")
script_dir <- paste0(project_dir, "scripts/")
func_dir <- paste0(script_dir, "functions/")
col_dir <- paste0(home_dir, "R/colour_palettes/")
genome_dir <- paste0(project_dir, "genome/")

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

if (!is.na(NF_cutoffs[1])) {
  out_dir <- paste0(
  	out_dir, "NF_hypo_", NF_cutoffs$MPNST_hypermethylated,
  	"/NF_hyper_", NF_cutoffs$MPNST_hypomethylated, "/"
  )
}

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

print(
  paste0(
    "Output candidate plots are ", 
    plot_dir, "MPNST_hypomethylated_marker_probes.png & ",
    plot_dir, "MPNST_hypermethylated_marker_probes.png"
  )
)


####################################################################################
### 0. Load packages, functions and colours ###
####################################################################################

library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)
library(cowplot)
library(rtracklayer)

DMR_analysis <- dget(paste0(func_dir, "DMR_analysis.R"))
DMR_filter <- dget(paste0(func_dir, "DMR_filter.R"))
filter_background <- dget(paste0(func_dir, "filter_background.R"))
fetch_DM_beta <- dget(paste0(func_dir, "fetch_DM_beta.R"))

col_pal <- read.table(
  paste0(col_dir, "labelled_colour_palette_2.txt"),
  sep = "\t",
  header = F,
  comment.char = ""
)


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
orig_normal_mtx <- lapply(normal_se, function(x) assays(x)[[1]])
# rename elements:
names(orig_normal_mtx) <- paste0(
  gsub(
    "_.*$", "", 
    gsub(
    	"^.*-", "", 
    	names(orig_normal_mtx)
    )
  ),
  "_non_malig"
)

# merge those with less than 5 samples:
if (exists("mixed_mtx")) {
  rm(mixed_mtx)
}

for (i in 1:length(orig_normal_mtx)) {

  if (ncol(orig_normal_mtx[[i]]) < 5) {

  	if (!exists("mixed_mtx")) {
  	  mixed_mtx <- orig_normal_mtx[[i]]
  	  to_rm <- names(orig_normal_mtx)[i]
  	} else {
  	  mixed_mtx <- cbind(mixed_mtx, orig_normal_mtx[[i]])
  	  to_rm <- c(to_rm, names(orig_normal_mtx)[i])
  	}

  }
}
normal_mtx <- orig_normal_mtx[!(names(orig_normal_mtx) %in% to_rm)]
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

print(
  paste0(
    unfiltered_hyper[16], " hypermethylated and ", unfiltered_hypo[16],
    " hypomethylated beta value/p-value filtered candiates identified..."
  )
)


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
    "Filtering for canditates hypermethylated in MPNST with median normal ", 
    "background level values < ", 
    background_cutoff, 
    " and hypomethylated in MPNST with median MPNST background values > ", 
    background_cutoff, "..."
  )
)

for (i in 1:length(DMR_initial)) {

  initial_filt <- lapply(DMR_initial[[i]], function(x) {

    if (nrow(x) > 0) {

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

    } else {

      initial_filt <- NA

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
  	bg_filt_hyper[["all_non_malig"]], " hypermethylated and ", 
    bg_filt_hypo[["all_non_malig"]], " hypomethylated high ", 
    "background filtered candiates identified..."
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

	  if (!is.na(y)) {

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

	  } else {

	  	return(y)

	  }

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
  	blood_filt_hyper[["all_non_malig"]], " hypermethylated and ", 
    blood_filt_hypo[["all_non_malig"]]," hypomethylated blood ",
    "filtered candiates identified..."
  )
)


####################################################################################
### 4. Filter for probes DM in MPNST vs each normal ###
####################################################################################

if (check_indiv_tissue) {

  # separate hyper and hypo dfs:
  DMR_list <- list(
    hyper = lapply(DMR_blood_filt, function(x) x$hyper),
    hypo = lapply(DMR_blood_filt, function(x) x$hypo)
  )
  
  # find DMR probes common to all normal tissue:
  DMR_probes <- lapply(DMR_list, function(x) {
    lapply(x, rownames)
  })
  DMR_indiv_probes <- lapply(DMR_probes, function(x) {
  
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
  
  # fetch indiv DM values:
  if (exists("DMR_indiv")) {
    rm(DMR_indiv)
  }
  for (i in 1:length(DMR_list)) {
  
    if (length(DMR_list[[i]]) > 0) {
  
     if (!exists("DMR_indiv")) {
  
       DMR_indiv <- list(
         lapply(DMR_list[[i]], function(y) {
            y[rownames(y) %in% DMR_indiv_probes[[i]],]
          })
        )
        names(DMR_indiv) <- names(DMR_list)[i]
  
     } else {
  
       DMR_indiv[[i]] <- lapply(DMR_list[[i]], function(y) {
          y[rownames(y) %in% DMR_indiv_probes[[i]],]
        })
        names(DMR_indiv)[i] <- names(DMR_list)[i]
  
     }
     
    }
    
  }
  
  # record values:
  indiv_hyper <- unlist(
    lapply(DMR_indiv$hyper, function(x) {
      nrow(x)
    })
  )
  
  indiv_hypo <- unlist(
    lapply(DMR_indiv$hypo, function(x) {
      nrow(x)
    })
  )
  
  DMR_record <- rbind(
    DMR_record,
    rbind(indiv_hyper, indiv_hypo)
  )
  
  print(
    paste0(
     indiv_hyper[16], " hypermethylated and ", indiv_hypo[16],
     " hypomethylated final candiates identified"
    )
  )

  prefinal_DMR <- DMR_indiv$all_non_malig

} else {

  prefinal_DMR <- DMR_blood_filt$all_non_malig

}


####################################################################################
### 5. Filter by NF methylation  ###
####################################################################################

# load NF methylation bed file:
if (!file.exists(paste0(Robject_dir, "NF_methylation_ranges.Rdata"))) {

  NF_met <- import(paste0(raw_dir, "feber_2011/GSM541730_benign_batman.gff"))
  NF_gr <- GRanges(
    seqnames <- paste0("chr", seqnames(NF_met)),
    ranges <- ranges(NF_met),
    strand <- strand(NF_met)
  )
  values(NF_gr) <- values(NF_met)

  saveRDS(NF_gr, paste0(Robject_dir, "NF_methylation_ranges.Rdata"))

} else {

  NF_gr <- readRDS(paste0(Robject_dir, "NF_methylation_ranges.Rdata"))

}

# fetch positions of each probe:
probe_coords <- normal_se[[1]]@rowRanges

# annotate DMR dfs with probe info:
prefinal_DMR <- lapply(prefinal_DMR, function(x) {

  if (!is.na(x)) {

    # make probe column:
    x$probe <- rownames(x)
    # fetch probe locations and add:
    specific_coords <- probe_coords[
      names(probe_coords) %in% x$probe
    ]
    return(
      merge(
        x,
        data.frame(
          probe = names(specific_coords),
          chr = seqnames(specific_coords),
          coord = start(specific_coords)
        ),
        by = "probe"
      )
    )

  } else {
    return(x)
  }

})

# import gencode annot:
gencode <- import(
  paste0(genome_dir, "gencode.v19.annotation.gtf")
)

# filter hypermethylated candidates for those hypomethylated in NF, and vice versa:
for (i in 1:length(prefinal_DMR)) {

  if (nrow(prefinal_DMR[[i]]) > 0) {

    # convert DM probe df to granges object:
    dm_gr <- GRanges(
      seqnames = prefinal_DMR[[i]]$chr,
      ranges = IRanges(
        start = prefinal_DMR[[i]]$coord,
        end = prefinal_DMR[[i]]$coord
      ),
      strand = "*"
    )
    values(dm_gr) <- subset(
      prefinal_DMR[[i]], 
      select = -c(chr, coord)
    )

    # find overlaps with genes and keep in separate gr:
    gene_olaps <- findOverlaps(gencode, dm_gr)

    if (i==1) {
      DMR_genes <- list(
      	hyper = gencode[queryHits(gene_olaps)]
      )
      DMR_genes$hyper$probe <- dm_gr$probe[subjectHits(gene_olaps)]
    } else {
      DMR_genes$hypo <- gencode[queryHits(gene_olaps)]
      DMR_genes$hypo$probe <- dm_gr$probe[subjectHits(gene_olaps)]
    }
  
    if (length(grep("Hyper", prefinal_DMR[[i]]$status[i])) > 0) {
  
      # isolate hypomethylated NF regions:
      NF_spec <- NF_gr[NF_gr$score < NF_cutoffs$MPNST_hypermethylated,]
  
      # find overlaps with gr:
      olaps <- findOverlaps(dm_gr, NF_spec)
  
    } else {
  
      # isolate hypermethylated NF regions:
      NF_spec <- NF_gr[NF_gr$score > NF_cutoffs$MPNST_hypomethylated,]
  
      # find overlaps with gr:
      olaps <- findOverlaps(dm_gr, NF_spec)
  
    }

    # combine array and MeDIP data:
    DMR_res <- as.data.frame(dm_gr[queryHits(olaps),])
    DMR_res <- cbind(
      DMR_res,
      subset(values(NF_spec[subjectHits(olaps),]), select = score)
    )
  
    # return filtered candidates:
    if (i==1) {
      DMR_final <- list(DMR_res)
    } else {
      DMR_final[[i]] <- DMR_res
    }

  } else {

    if (i==1) {
      DMR_final <- list(NULL)
    } else {
      DMR_final[[i]] <- NULL
    }

  }
  
}
names(DMR_final) <- names(prefinal_DMR)

print(
  paste0(
   nrow(DMR_final$hyper), " hypermethylated and ", nrow(DMR_final$hypo),
   " hypomethylated final candiates identified"
  )
)


####################################################################################
### 6. Prepare plot data ###
####################################################################################

# calculate mean values for all probes and tissues:
if (!file.exists(paste0(Robject_dir, "all_beta_stats.Rdata"))) {

  all_beta_df <- orig_normal_mtx
  all_beta_df$MPNST <- MPNST_beta_df
  all_beta_df$blood <- blood_mtx

  all_beta_stats <- lapply(all_beta_df, function(x) {

    return(
      t(
      	apply(x, 1, function(y) {
          # calculate mean, 10% and 90% quartiles:
          quantiles <- quantile(
          	y, c(0.1, 0.25, 0.75, 0.9), na.rm = TRUE
          )
          return(
            c(
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

  })

#  ######
#
#  temp_beta_stats <- lapply(all_beta_df, function(x) {
#
#    return(
#      t(
#      	apply(x, 1, function(y) {
#          # calculate mean, 10% and 90% quartiles:
#          quantiles <- quantile(y, c(0.25, 0.75), na.rm = TRUE)
#          return(
#            c(
#              quart_0.25 = quantiles[1],
#              quart_0.75 = quantiles[2],
#              median = median(y, na.rm = TRUE)
#         	)
#          )
#        })
#      )
#    )
#
#  })
#
#  for (i in 1:length(all_beta_stats)) {
#  	if (i==1) {
#  	  all_beta_stats_temp <- list(
#  	  	cbind(all_beta_stats[[i]], temp_beta_stats[[i]])
#  	  )
#  	} else {
#  	  all_beta_stats_temp[[i]] <- cbind(all_beta_stats[[i]], temp_beta_stats[[i]])
#  	}
#  }
#
#  all_beta_stats <- all_beta_stats_temp
#  names(all_beta_stats) <- names(all_beta_df)
#  ######
  
  saveRDS(all_beta_stats, paste0(Robject_dir, "all_beta_stats.Rdata"))

} else {

  all_beta_stats <- readRDS(paste0(Robject_dir, "all_beta_stats.Rdata"))

}

# fetch sample numbers for tissues:
n_vals <- c(
  unlist(lapply(orig_normal_mtx, ncol)),
  ncol(blood_mtx),
  as.character(10),
  ncol(MPNST_beta_df)
)
names(n_vals)[length(n_vals)-2] <- "blood"
names(n_vals)[length(n_vals)-1] <- "NF"
names(n_vals)[length(n_vals)] <- "MPNST"

# filter out tissues with < 10 samples from beta stats:
rm_tissues <- names(n_vals)[as.numeric(n_vals) < plot_min_sample_no]
plot_beta_stats <- all_beta_stats[
  !(names(all_beta_stats) %in% rm_tissues)
]

# update names for n_vals to plot labels:
m <- match(names(n_vals), tissue_key$code)
names(n_vals) <- tissue_key$name[m]

# fetch mean beta values, coords and overlapping genes for 
# DM probes and plot:
DMR_plot_df <- lapply(
  DMR_final,
  fetch_DM_beta,
  plot_beta_stats,
  probe_coords,
  tissue_key
)

# save final nos:
final_no <- lapply(DMR_final, nrow)

DMR_record <- rbind(
  DMR_record,
  t(
    data.frame(
      final_hyper = c(
        rep(NA, (ncol(DMR_record)-1)),
        final_no$hyper
      ),
      final_hypo = c(
        rep(NA, (ncol(DMR_record)-1)),
        final_no$hypo
      )
    )
  )
)

write.table(
  DMR_record,
  paste0(table_dir, "DMR_counts.txt"),
  row.names = T,
  col.names = T,
  quote = F
)


####################################################################################
### 7. Plot top 10 DM probes ###
####################################################################################

# find and take top 10 candidates by beta difference for plot:
DMR_top <- lapply(DMR_final, function(x) {
  ordered_x <- x[
    order(abs(x$mean.MPNST.minus.mean.all.non.malig), decreasing = T), 
  ]
  return(head(ordered_x, 10))
})

DMR_plot_df_top <- DMR_plot_df

for (i in 1:length(DMR_plot_df_top)) {

  DMR_plot_df_top[[i]] <- DMR_plot_df_top[[i]][
  	DMR_plot_df_top[[i]]$probe %in% DMR_top[[i]]$probe,
  ]

  # determine whether any final candidates overlap genes:
  top_DMR_genes <- DMR_genes[[i]][
    DMR_genes[[i]]$probe %in% DMR_plot_df_top[[i]]$probe
  ]
  
}

# plot data:
for (i in 1:length(DMR_plot_df_top)) {

  # make tissue column factor and set MPNST as level 1:
  nm_df <- DMR_plot_df_top[[i]][
    !(DMR_plot_df_top[[i]]$tissue %in% c("MPNST", "NF", "blood")),
  ]
  nm_df <- nm_df[order(nm_df$tissue),]

  # reorder so MPNST, NF and blood are top levels
  nm_blood <- rbind(
	DMR_plot_df_top[[i]][DMR_plot_df_top[[i]]$tissue == "blood",],
	nm_df
  )
  nm_nf <- rbind(
	DMR_plot_df_top[[i]][DMR_plot_df_top[[i]]$tissue == "NF",],
	nm_blood
  )
  DMR_plot_df_top[[i]] <- rbind(
	DMR_plot_df_top[[i]][DMR_plot_df_top[[i]]$tissue == "MPNST",],
	nm_nf
  )

  DMR_plot_df_top[[i]]$tissue <- factor(DMR_plot_df_top[[i]]$tissue,
  	levels = unique(DMR_plot_df_top[[i]]$tissue))

  # add no samples to legend:
  legend_n_vals <- paste0(
  	levels(DMR_plot_df_top[[i]]$tissue),
  	" (",
  	n_vals[levels(DMR_plot_df_top[[i]]$tissue)],
  	")"
  )

  # update x_pos and make character factor:
  DMR_plot_df_top[[i]] <- arrange(
    DMR_plot_df_top[[i]], chr, coord
  )
  DMR_plot_df_top[[i]]$x_pos <- as.character(
    rleid(DMR_plot_df_top[[i]]$probe)
  )

  # set up temp df for labels:
  temp_df <- DMR_plot_df_top[[i]][
    !duplicated(DMR_plot_df_top[[i]]$probe),
  ]

  # plot data:
  p <- ggplot(
    DMR_plot_df_top[[i]], 
    aes(
      x=x_pos, 
      y=median, 
      color=factor(
      	tissue, labels = legend_n_vals
      ),
      group=shape
    )
  )
#  p <- p + geom_jitter(
#  	aes(shape=shape),
#  	size = 2, 
#  	width = 0.1, 
#  	height = 0
#  )
  p <- p + geom_pointrange(
  	aes(ymin = quart_0.25, ymax = quart_0.75, shape=shape), 
  	size = 0.5,
    position=position_jitter(width=0.3, seed = 754)
  )
  p <- p + labs(color = "Tissue")
#  p <- p + scale_x_discrete(
#  	labels = c(
#  	  "candidate 1", 
#  	  seq(2, nrow(DMR_plot_df_top[[i]]), 1)
#  	)
#    labels = paste0(
#      temp_df$probe, 
#      "\n", 
#      temp_df$chr,
#      ":",
#      temp_df$coord
#    )
#  )
  p <- p + guides(shape = FALSE)
  p <- p + xlab(
  	paste0(
  	  "CpG probes ", names(DMR_plot_df_top)[i], "methylated ",
  	  "in MPNST"
  	)
  )
  p <- p + ylab("Median methylation level")
  p <- p + theme_cowplot(12)
  p <- p + theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
  p <- p + scale_color_manual(
  	values = col_pal[1:length(unique(DMR_plot_df_top[[i]]$tissue)),1]
  )
  
  print(
    paste0(
      "Plotting candidates, filename = ", plot_dir, "MPNST_", names(DMR_plot_df_top)[i], "methylated_marker_probes.png"
    )
  )

  png(
  	paste0(
  	  plot_dir, "top_10_MPNST_", names(DMR_plot_df_top)[i], 
  	  "methylated_marker_probes_median_quartile_bars.png"
  	),
  	res = 300,
  	unit = "in",
  	height = 7,
  	width = 15
  )
    print(p)
  dev.off()
  
  ggsave(
    paste0(plot_dir, "top_10_MPNST_", names(DMR_plot_df_top)[i], 
    	"methylated_marker_probes_median_quartile_bars.pdf"),
    plot = p,
    useDingbats = F,
    width = 10
  )

}


####################################################################################
### 6. Plot DM probes ###
####################################################################################

for (i in 1:length(DMR_plot_df)) {

  # make tissue column factor and set MPNST as level 1:
  nm_df <- DMR_plot_df[[i]][
    !(DMR_plot_df[[i]]$tissue %in% c("MPNST", "NF")),
  ]
  nm_df <- nm_df[order(nm_df$tissue),]

  nm_mpnst <- rbind(
	DMR_plot_df[[i]][DMR_plot_df[[i]]$tissue == "MPNST",],
	nm_df
  )
  DMR_plot_df[[i]] <- rbind(
  	nm_mpnst,
  	DMR_plot_df[[i]][DMR_plot_df[[i]]$tissue == "NF",]
  )

  DMR_plot_df[[i]]$tissue <- factor(DMR_plot_df[[i]]$tissue,
  	levels = unique(DMR_plot_df[[i]]$tissue))

  # add no samples to legend:
  legend_n_vals <- paste0(
  	levels(DMR_plot_df[[i]]$tissue),
  	" (",
  	n_vals[levels(DMR_plot_df[[i]]$tissue)],
  	")"
  )

  # update x_pos and make character factor:
  DMR_plot_df[[i]] <- arrange(
    DMR_plot_df[[i]], chr, coord
  )
  DMR_plot_df[[i]]$x_pos <- as.character(
    rleid(DMR_plot_df[[i]]$probe)
  )

  # set up temp df for labels:
  temp_df <- DMR_plot_df[[i]][
    !duplicated(DMR_plot_df[[i]]$probe),
  ]

  p <- ggplot(
    DMR_plot_df[[i]], 
    aes(x=x_pos, y=beta, color=tissue, group=shape)
  )
  p <- p + geom_jitter(
  	aes(shape=shape),
  	size = 2, 
  	width = 0.1, 
  	height = 0
  )
  p <- p + scale_x_discrete(
    labels = paste0(
      unique(DMR_plot_df[[i]]$probe), 
      "\n", 
      unique(DMR_plot_df[[i]]$chr),
      ":",
      unique(DMR_plot_df[[i]]$coord)
    )
  )
  p <- p + xlab(
  	paste0(
  	  "CpG probes ", names(DMR_plot_df)[i], "methylated ",
  	  "in MPNST"
  	)
  )
  p <- p + ylab("Methylation level")
  p <- p + theme_cowplot(12)
  p <- p + theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
  p <- p + scale_color_manual(
  	values = col_pal[1:length(unique(DMR_plot_df[[i]]$tissue)),1]
  )
  
  print(
    paste0(
      "Plotting candidates, filename = ", plot_dir, "MPNST_", names(DMR_plot_df)[i], "methylated_marker_probes.png"
    )
  )

  png(
  	paste0(
  	  plot_dir, "MPNST_", names(DMR_plot_df)[i], 
  	  "methylated_marker_probes.png"
  	),
  	res = 300,
  	unit = "in",
  	height = 7,
  	width = 10
  )
    print(p)
  dev.off()
  
  ggsave(
    paste0(plot_dir, "MPNST_", names(DMR_plot_df)[i], 
    	"methylated_marker_probes.pdf"),
    plot = p,
    useDingbats = F
  )

}
