
min_sample_per_probe <- 5
min_beta_val <- 0.2
subset_mtx = F

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/MPSNT_ctMet/")
raw_dir <- paste0(project_dir, "raw_files/")
script_dir <- paste0(project_dir, "scripts/")
func_dir <- paste0(script_dir, "functions/")

normal_in <- paste0(raw_dir, "TCGA_normal/")
Robject_in <- paste0(raw_dir, "Rdata/")
result_dir <- paste0(project_dir, "results/")

if (sub) {
  result_dir <- paste0(result_dir, "subset/")
}

Robject_dir <- paste0(result_dir, "Rdata/")
plot_dir <- paste0(result_dir, "plots/")
table_dir <- paste0(result_dir, "tables/")

system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", plot_dir))
system(paste0("mkdir -p ", table_dir))


####################################################################################
### 0. Load packages and functions ###
####################################################################################

library("TCGAbiolinks")
library("SummarizedExperiment")

DMR_analysis <- dget(paste0(func_dir, "DMR_analysis.R"))


####################################################################################
### 1. Identify MPNST vs normal DMRs ###
####################################################################################

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

# Load MPNST beta values:
MPNST_beta_df <- readRDS(paste0(Robject_in, "MPNST_beta_df.Rdata"))

# label columns:
colnames(MPNST_beta_df) <- paste0(colnames(MPNST_beta_df), "-MPNST")

if (subset_mtx) {

  normal_mtx <- lapply(normal_mtx, function(x) x[1:100,])

  MPNST_beta_df <- MPNST_beta_df[1:100,]

}

# find DMRs between each group and MPNST groups:
for (i in 1:length(normal_mtx)) {

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

# filter for probes with at least 0.35 and p-value < 10e-5:
MPNST_vs_all_normal_filt <- MPNST_vs_normal$all_non_malig[
  MPNST_vs_normal$all_non_malig$p.value.adj.MPNST.all.non.malig < 10e-5 &
  MPNST_vs_normal$all_non_malig$mean.MPNST.minus.mean.all.non.malig >= 0.35,
]

MPNST_vs_all_normal_hyper <- MPNST_vs_all_normal_filt[
  MPNST_vs_all_normal_filt$mean.MPNST.minus.mean.all.non.malig > 0,
]

MPNST_vs_all_normal_hypo <- MPNST_vs_all_normal_filt[
  MPNST_vs_all_normal_filt$mean.MPNST.minus.mean.all.non.malig < 0,
]


####################################################################################
### 2. Remove hypo/hypermethylated blood regions ###
####################################################################################





####################################################################################
### 4. Identify stringent DMRs ###
####################################################################################

# identify methylation probes with 0.2 < mean beta values < 0.3 and others > 0.8:
both_beta <- lapply(filt_mtx, function(x) { 

  mean_beta <- apply(x, 1, mean, na.rm = TRUE)

  low_beta <- names(mean_beta[mean_beta > 0.2 & mean_beta < 0.3])
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

# fetch list of probes with 0.2 < mean beta values < 0.3 in one tissue and > 0.8 
# in other tissue:
dm_ind <- list(
  MPNST_high = which(
    both_beta$MPNST$high_beta %in% both_beta$BRCA$low_beta
  ),
  MPNST_low = which(
    both_beta$MPNST$low_beta %in% both_beta$BRCA$high_beta
  )
)

dm_probes <- list(
  MPNST_high = both_beta$MPNST$high_beta[dm_ind$MPNST_high],
  MPNST_low = both_beta$MPNST$low_beta[dm_ind$MPNST_low]
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








