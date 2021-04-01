

# define and create directories:
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/MPNST_ctMet/")
script_dir <- paste0(project_dir, "scripts/")
func_dir <- paste0(script_dir, "functions/")

result_dir <- paste0(project_dir, "results/")


####################################################################################
### 0. Load packages and functions ###
####################################################################################


####################################################################################
### 1. Load counts ###
####################################################################################

print("Loading candidate counts for all parameter sets...")

filenames <- list.files(
  result_dir, 
  pattern = "DMR_counts.txt",
  recursive = T,
  full.names = T
)

for (i in 1:length(filenames)) {

  params <- strsplit(
  	gsub(
  	  "^.*results/", "",
  	  filenames[i]
  	), 
    "/"
  )[[1]]
  params <- params[grep(".txt|tables", params, invert = T)]
  params <- params[params != ""]

  df <- read.table(
  	paste0(filenames[i]),
  	header = T,
  	stringsAsFactors = F
  )

  param_vals <- gsub("^.*_", "", params)

  if (i==1) {

  	DMR_vs_non_malig_counts <- data.frame(
  	  row.names = c(
  	  	gsub("_[0-9].*$", "", params),
  	  	"hyper_candidates", "hypo_candidates"
  	  ),
  	  c(
  	  	param_vals,
  	  	df[rownames(df) == "final_hyper",]$all_non_malig,
  	  	df[rownames(df) == "final_hypo",]$all_non_malig
  	  )
  	)

  } else {

    DMR_vs_non_malig_counts <- cbind(
      DMR_vs_non_malig_counts,
      data.frame(
        row.names = c(
          gsub("_[0-9].*$", "", params),
          "hyper_candidates", "hypo_candidates"
        ),
        c(
          param_vals,
          df[rownames(df) == "final_hyper",]$all_non_malig,
          df[rownames(df) == "final_hypo",]$all_non_malig
        )
      )
    )

  }

}
colnames(DMR_vs_non_malig_counts) <- 1:ncol(DMR_vs_non_malig_counts)

