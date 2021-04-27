
min_sample_per_probe <- 5
subset_data <- TRUE

## fetch other variables:
#args = commandArgs(trailingOnly=TRUE)
#
#merge_tissues_min <- args[1]
#beta_cutoff <- args[2]
#print(beta_cutoff)
#pval_cutoff <- as.character(args[3])
#background_cutoff <- args[4]
#blood_cutoffs <- list(
#  MPNST_hypermethylated = args[5],
#  MPNST_hypomethylated = args[6]
#)
#NF_cutoffs <- list(
#  MPNST_hypermethylated = args[7],
#  MPNST_hypomethylated = args[8]
#)
#check_indiv_tissue <- as.logical(args[9])
#plot_min_sample_no <- as.numeric(args[10])

## @knitr merge_tissues_min

merge_tissues_min <- 5
beta_cutoff <- 0.35
pval_cutoff <- as.character("1e-5")
background_cutoff <- 0.2
blood_cutoffs <- list(
  MPNST_hypermethylated = 0.2,
  MPNST_hypomethylated = 0.8
)
NF_cutoffs <- list(
  MPNST_hypermethylated = 0.2,
  MPNST_hypomethylated = 0.8
)
check_indiv_tissue <- as.logical("FALSE")
plot_min_sample_no <- 10

print(paste0("beta_cutoff = ", beta_cutoff))
print(paste0("pval_cutoff = ", pval_cutoff))
print(paste0("background_cutoff = ", background_cutoff))
print(
  paste0(
    "blood_cutoffs = ", 
    blood_cutoffs$MPNST_hypermethylated, ", ", 
    blood_cutoffs$MPNST_hypomethylated
  )
)
print(
  paste0(
    "NF_cutoffs = ", 
    NF_cutoffs$MPNST_hypermethylated, ", ", 
    NF_cutoffs$MPNST_hypermethylated
  )
)

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
genome_dir <- paste0(project_dir, "genome/")
ref_dir <- paste0(project_dir, "refs/")

normal_in <- paste0(raw_dir, "TCGA_normal/")
Robject_in <- paste0(raw_dir, "Rdata/")
result_dir <- paste0(project_dir, "results2/")
out_path <- paste0(
  result_dir,
  "beta_", beta_cutoff, "_",
  "pval_", pval_cutoff, "/",
  "background_", background_cutoff, "/",
  "blood_hypo_", blood_cutoffs$MPNST_hypermethylated, "_",
  "blood_hyper_", blood_cutoffs$MPNST_hypomethylated, "/"
)

if (!is.na(NF_cutoffs[1])) {
  out_path <- paste0(
  	out_path, "NF_hypo_", NF_cutoffs$MPNST_hypermethylated,
  	"_NF_hyper_", NF_cutoffs$MPNST_hypomethylated, "/"
  )
}

if (subset_data) {
  result_dir <- paste0(result_dir, "subset/")
  out_path <- paste0(out_path, "subset/")
}

Robject_dir <- paste0(result_dir, "Rdata/")
plot_dir <- paste0(out_path, "plots/")
table_dir <- paste0(out_path, "tables/")

system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", plot_dir))
system(paste0("mkdir -p ", table_dir))

print(
  paste0(
    "Output candidate plots are in ", 
    plot_dir
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
library(data.table)

DMR_analysis <- dget(paste0(func_dir, "DMR_analysis.R"))
DMR_filter <- dget(paste0(func_dir, "DMR_filter.R"))
background_filter <- dget(paste0(func_dir, "background_filter.R"))
record_vals <- dget(paste0(func_dir, "record_vals.R"))
blood_filter <- dget(paste0(func_dir, "blood_filter.R"))
annotate_probes <- dget(paste0(func_dir, "annotate_probes.R"))
filter_by_NF <- dget(paste0(func_dir, "filter_by_NF.R"))
fetch_DM_beta <- dget(paste0(func_dir, "fetch_DM_beta.R"))
plot_DMR <- dget(paste0(func_dir, "plot_DMR.R"))

col_pal <- read.table(
  paste0(ref_dir, "custom_colour_palette.txt"),
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

  if (subset_data) {
    normal_se_files <- normal_se_files[1]
  }
  
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

if (subset_data) {

  normal_mtx <- lapply(normal_mtx, function(x) x[1:40000,])

  MPNST_df <- MPNST_df[1:40000,]

}

# find DMRs between each group and MPNST groups:
if (
  !file.exists(
    paste0(Robject_dir, "MPNST_vs_normal_groups.Rdata")
  )
) {
	
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
DMR_record <- as.data.frame(
  rbind(
    sapply(MPNST_vs_normal, function(x) {
      nrow(
        x[x$status == "Hypermethylated in MPNST",]
      )
    }), 
    sapply(MPNST_vs_normal, function(x) {
      nrow(
        x[x$status == "Hypomethylated in MPNST",]
      )
    })
  )
)
rownames(DMR_record) <- c("unfiltered_hyper", "unfiltered_hypo")

print(
  paste0(
    DMR_record["unfiltered_hyper",]$all_non_malig, 
    " hypermethylated and ",
    DMR_record["unfiltered_hypo",]$all_non_malig,
    " hypomethylated unfiltered candiates ", 
    "identified..."
  )
)


####################################################################################
### 2. Filter DM probes by beta/pvals ###
####################################################################################

print(
  paste0(
    "Filtering for canditates with at least ", 
    beta_cutoff, 
    " beta value difference and < ", pval_cutoff, "p-value..."
  )
)

DMR_initial <- lapply(
  MPNST_vs_normal,
  DMR_filter,
  beta_cutoff,
  as.numeric(pval_cutoff)
)

# record values:
DMR_record <- rbind(
  DMR_record,
  record_vals(DMR_initial, "beta_filt")
)


####################################################################################
### 2. Filter DM probes for low background candidates ###
####################################################################################

print(
  paste0(
    "Filtering for canditates hypermethylated in MPNST with median ",
    "normal background level values < ", 
    background_cutoff, 
    " and hypomethylated in MPNST with median MPNST background ", 
    "values > ", background_cutoff, "..."
  )
)

for (i in 1:length(DMR_initial)) {

  initial_filt <- lapply(DMR_initial[[i]], function(x) {

    if (nrow(x) > 0) {

      if (length(grep("Hypermethylated", x$status[1])) > 0) {

        return(
          background_filter(
            dmr = x,
            background_cutoff,
            hypo_beta = normal_mtx[[i]],
            min_sample_per_probe
          )
        )
  
      } else {
  
        return(
        background_filter(
            dmr = x,
            background_cutoff,
            hypo_beta = MPNST_df,
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
DMR_record <- rbind(
  DMR_record,
  record_vals(DMR_bg_filt, "bg_filt")
)


####################################################################################
### 3. Remove hypo/hypermethylated blood regions ###
####################################################################################

print("Loading blood beta values...")

# load GEO blood data and fetch beta matrix:
if (
  !file.exists(
    paste0(Robject_dir, "blood_mtx.Rdata")
  )
) {

  blood_se <- readRDS(paste0(raw_dir, "GEO_blood/SE_GEO_Meth_Data.rds"))
  blood_mtx <- assays(blood_se)[[1]][-1,]
  
  if (subset_data) {
    blood_mtx <- blood_mtx[1:40000,]
  }

  saveRDS(blood_mtx, paste0(Robject_dir, "blood_mtx.Rdata"))

} else {
  blood_mtx <- readRDS(paste0(Robject_dir, "blood_mtx.Rdata"))
}

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
  return(lapply(x, blood_filter, blood_mtx, min_sample_per_probe))
})


# record values:
DMR_record <- rbind(
  DMR_record,
  record_vals(DMR_blood_filt, "blood_filt")
)


####################################################################################
### 4. Filter for probes DM in MPNST vs each normal ###
####################################################################################

if (check_indiv_tissue) {

  DMR_indiv <- indiv_tissue_filter(DMR_blood_filt)

  # record values:
  DMR_record <- rbind(
    DMR_record,
    record_vals(DMR_indiv, "indiv_filt")
  )

  prefinal_DMR <- DMR_indiv$all_non_malig

} else {

  prefinal_DMR <- DMR_blood_filt$all_non_malig

}


####################################################################################
### 5. Filter by NF methylation  ###
####################################################################################

# load NF methylation bed file:
if (
  !file.exists(
    paste0(Robject_dir, "NF_methylation_ranges.Rdata")
  )
) {

  NF_met <- import(
    paste0(raw_dir, "feber_2011/GSM541730_benign_batman.gff")
  )
  NF_gr <- GRanges(
    seqnames <- paste0("chr", seqnames(NF_met)),
    ranges <- ranges(NF_met),
    strand <- strand(NF_met)
  )
  values(NF_gr) <- values(NF_met)

  saveRDS(NF_gr, paste0(Robject_dir, "NF_methylation_ranges.Rdata"))

} else {
  NF_gr <- readRDS(
    paste0(Robject_dir, "NF_methylation_ranges.Rdata")
  )
}

# fetch positions of each probe:
probe_coords <- normal_se[[1]]@rowRanges

# annotate DMR dfs with probe info:
prefinal_DMR <- lapply(prefinal_DMR, annotate_probes)

# import gencode annot:
if (
  !file.exists(
    paste0(genome_dir, "gencode.v19.annotation.Rdata")
  )
) {

  gencode <- import(
    paste0(genome_dir, "gencode.v19.annotation.gtf")
  )

  saveRDS(
    gencode, paste0(genome_dir, "gencode.v19.annotation.Rdata")
  )

} else {
  gencode <- readRDS(
    paste0(genome_dir, "gencode.v19.annotation.Rdata")
  )
}

# filter hypermethylated candidates for those hypomethylated in NF, 
# and vice versa:
final_DMR_and_genes <- filter_by_NF(
  dmr = prefinal_DMR,
  NF_vals = NF_gr, 
  gene_annot = gencode
)
final_DMR <- final_DMR_and_genes$final_DMR
DMR_genes <- final_DMR_and_genes$DMR_genes

# record and write values to file:
DMR_record <- rbind(
  DMR_record,
  record_vals(
  	dmr = final_DMR, 
  	label = "final",
  	col_names = colnames(DMR_record)
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
### 6. Prepare plot data ###
####################################################################################

# add all beta values to list:
all_beta_list <- orig_normal_mtx
all_beta_list$MPNST <- MPNST_df
all_beta_list$blood <- blood_mtx

# determine number of samples of each tissue:
n_vals <- c(
  unlist(lapply(orig_normal_mtx, ncol)),
  ncol(blood_mtx),
  as.character(10),
  ncol(MPNST_df)
)
names(n_vals)[length(n_vals)-2] <- "blood"
names(n_vals)[length(n_vals)-1] <- "NF"
names(n_vals)[length(n_vals)] <- "MPNST"

# update names for n_vals to plot labels:
m <- match(names(n_vals), tissue_key$code)
names(n_vals) <- tissue_key$name[m]

# fetch beta value stats for DM probes:
if (!file.exists(paste0(Robject_dir, "DMR_plot_df.Rdata"))) {

  DMR_plot_df <- fetch_DM_beta(
    dmr = final_DMR,
    beta_val_list = all_beta_list,
    probe_coords,
    tissue_key,
    n_vals,
    plot_min_sample_no
  )

  saveRDS(DMR_plot_df, paste0(Robject_dir, "DMR_plot_df.Rdata"))

} else {
  DMR_plot_df <- readRDS(paste0(Robject_dir, "DMR_plot_df.Rdata"))
}


####################################################################################
### 7. Plot top 10 and all DM probes ###
####################################################################################

# combine mean beta difference and adj p-values to order best 
# candidates and take top 10 candidates for plot:
DMR_top <- lapply(final_DMR, function(x) {

  # calculate log fold change between groups:
  log_fc <- log(x$mean.MPNST/x$mean.all.non.malig)

  # use Xiao, Bioinformatics 2014 to combine log FC and p-values:
  x$combined_score <- log_fc*(
  	-log10(x$p.value.adj.MPNST.all.non.malig)
  )

  # order by combined score:
  ordered_x <- x[
    order(abs(x$combined_score), decreasing = T), 
  ]
  return(head(ordered_x, 10))
})

DMR_plot_df_top <- DMR_plot_df

for (i in 1:length(DMR_plot_df_top)) {
  DMR_plot_df_top[[i]] <- DMR_plot_df_top[[i]][
  	DMR_plot_df_top[[i]]$probe %in% DMR_top[[i]]$probe,
  ]
}

# plot both all and top 10 DM probes:
plot_list <- list(
  top_10 = DMR_plot_df_top,
  all = DMR_plot_df
)

for (i in 1:length(plot_list)) {

  for (j in 1:length(plot_list[[i]])) {

    plot_DMR(
      plot_df = plot_list[[i]][[j]],
      dot_vals = "median",
      bar_vals = "quartile",
      plot_type = names(plot_list)[i],
      meth_type = names(plot_list[[i]])[j],
      probe_annot = FALSE
    )
  
    plot_DMR(
      plot_df = plot_list[[i]][[j]],
      dot_vals = "median",
      bar_vals = "wide_quartile",
      plot_type = names(plot_list)[i],
      meth_type = names(plot_list[[i]])[j],
      probe_annot = FALSE
    )
  
    plot_DMR(
      plot_df = plot_list[[i]][[j]],
      dot_vals = "mean",
      bar_vals = "sd",
      plot_type = names(plot_list)[i],
      meth_type = names(plot_list[[i]])[j],
      probe_annot = FALSE
    )
  
  }

}


