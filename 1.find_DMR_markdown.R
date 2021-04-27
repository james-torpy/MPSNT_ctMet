
merge_tissues_min <- 5
min_sample_per_probe <- 5

# define and create directories:
home_dir <- "/Users/torpor/clusterHome/"
#home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/MPNST_ctMet/")
raw_dir <- paste0(project_dir, "raw_files/")
script_dir <- paste0(project_dir, "scripts/")
func_dir <- paste0(script_dir, "functions/")
normal_in <- paste0(raw_dir, "TCGA_normal/")
Robject_in <- paste0(raw_dir, "Rdata/")
result_dir <- paste0(project_dir, "results2/DMR/")

Robject_dir <- paste0(result_dir, "Rdata/")
plot_dir <- paste0(result_dir, "plots/")
system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", plot_dir))


####################################################################################
### 0. Load packages, functions and colours ###
####################################################################################

library(TCGAbiolinks)
library(SummarizedExperiment)

DMR_analysis <- dget(paste0(func_dir, "DMR_analysis.R"))


####################################################################################
### 1. Identify MPNST vs normal DMRs ###
####################################################################################

if (
  !file.exists(
    paste0(Robject_dir, "MPNST_vs_normal_groups.Rdata")
  )
) {

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
  
    if (ncol(orig_normal_mtx[[i]]) < merge_tissues_min) {
  
      if (!exists("mixed_mtx")) {
        mixed_mtx <- orig_normal_mtx[[i]]
        to_rm <- names(orig_normal_mtx)[i]
      } else {
        mixed_mtx <- cbind(mixed_mtx, orig_normal_mtx[[i]])
        to_rm <- c(to_rm, names(orig_normal_mtx)[i])
      }
  
    }
  }
  
  if (exists("to_rm")) {
    normal_mtx <- orig_normal_mtx[
      !(names(orig_normal_mtx) %in% to_rm)
    ]
    normal_mtx$mixed <- mixed_mtx
  } else {
    normal_mtx <- orig_normal_mtx
  }
  
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
  MPNST_df <- readRDS(paste0(Robject_in, "MPNST_beta_df.Rdata"))
  
  # label columns:
  colnames(MPNST_df) <- paste0(colnames(MPNST_df), "-MPNST")

  # find no common probes:
  common_probes <- rownames(normal_mtx$all_non_malig)[
    rownames(normal_mtx$all_non_malig) %in% rownames(MPNST_df)
  ]
  
  # find DMRs between each group and MPNST groups: 
  for (i in 1:length(normal_mtx)) {

    print(
      paste0(
        "Identifying DM probes, MPNST vs ", 
         names(normal_mtx)[i], "..."
      )
    )

    DMR <- DMR_analysis(
      sample1_beta = MPNST_df, 
      sample1_name = "MPNST", 
      sample2_beta = normal_mtx[[i]], 
      sample2_name = names(normal_mtx)[i],
      row_ranges = rowRanges(normal_se[[1]]),
      min_sample_per_probe,
      plot_dir
    )
  
    if (i==1) {
     MPNST_vs_normal <- list(DMR$DMR)
     filtered_probes <- list(DMR$filtered_probes)
    } else {
     MPNST_vs_normal[[i]] <- DMR$DMR
     filtered_probes[[i]] <- DMR$filtered_probes
    }
  
  }
  names(MPNST_vs_normal) <- names(normal_mtx)
  names(filtered_probes) <- names(normal_mtx)

  # save DMR:
  DMR_list <- list(
    MPNST_vs_normal = MPNST_vs_normal,
    non_malig_probes = nrow(normal_mtx$all_non_malig),
    MPNST_probes = nrow(MPNST_df),
    common_probes = length(common_probes),
    filtered_probes = filtered_probes$all_non_malig
  )

  saveRDS(
   DMR_list, 
   paste0(Robject_dir, "MPNST_vs_normal_groups.Rdata")
  )

} else {

  print("DMR previously performed, exiting...")

}




