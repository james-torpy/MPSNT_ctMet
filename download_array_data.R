library(GEOquery)
library(dplyr)
library(GenomicRanges)
library(TCGAbiolinks)
library(SummarizedExperiment)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/MPNST_ctMet/")
raw_dir <- paste0(project_dir, "raw_files/")
Robject_dir <- paste0(raw_dir, "Rdata/")

system(paste0("mkdir -p ", Robject_dir))


####################################################################################
### 1. Download Jour et. al. 2019 data
####################################################################################

if (!(file.exists(paste0(Robject_dir, "MPNST_beta_1.Rdata")))) {

  jour_soft <- getGEO("GSE112308", GSEMatrix=T)

  # fetch phenotype data:
  jour_pheno <- list(
    phenotype = subset(
      pData(jour_soft[[1]]), 
      select = c(
        "title", "geo_accession", "type", "characteristics_ch1", "label_ch1", 
        "data_row_count"
      )
    ),
    protocols = subset(
      pData(jour_soft[[1]]), 
      select = c(
        "title", "geo_accession", "extract_protocol_ch1", "label_protocol_ch1", 
        "hyb_protocol", "scan_protocol", "data_processing"
      )
    )
  )
  
  
  # fetch beta value matrix:
  jour_beta <- jour_soft$GSE112308_series_matrix.txt.gz@assayData$exprs
  
  # isolate MPNST samples only:
  jour_pheno <- lapply(jour_pheno, function(x) {
    x <- x[grep("MPNST", x$title),]
    x$title <- gsub("^.* ", "", x$title)
    return(x)
  })
  
  jour_beta <- jour_beta[,colnames(jour_beta) %in% jour_pheno$phenotype$geo_accession]
  
  # initiate list for MPNST data:
  MPNST_beta <- list(
    jour = jour_beta
  )
  saveRDS(MPNST_beta, paste0(Robject_dir, "MPNST_beta_1.Rdata"))
  
  MPNST_pheno <- list(
    jour_pheno$phenotype
  )
  saveRDS(MPNST_pheno, paste0(Robject_dir, "MPNST_pheno_1.Rdata"))
  
  MPNST_protocols <- list(
    jour_pheno$protocols
  )
  saveRDS(MPNST_protocols, paste0(Robject_dir, "MPNST_protocols_1.Rdata"))

} else {

  MPNST_beta <- readRDS(paste0(Robject_dir, "MPNST_beta_1.Rdata"))
  
  MPNST_pheno <- readRDS(paste0(Robject_dir, "MPNST_pheno_1.Rdata"))
  
  MPNST_protocols <- readRDS(paste0(Robject_dir, "MPNST_protocols_1.Rdata"))

}


####################################################################################
### 2. Download TCGA MPNST data
####################################################################################

# fetch barcodes for MPNST samples:
sarc_meta <- read.table(
  paste0(raw_dir, "TCGA_MPNST/combined_study_clinical_data.tsv"),
  sep = "\t",
  header = T,
  stringsAsFactors = F
)
sarc_ids <- subset(
  sarc_meta, 
  select = c("Study.ID", "Sample.ID", "Cancer.Type.Detailed")
)
MPNST_ids <- sarc_ids[
  sarc_ids$Cancer.Type.Detailed == "Malignant Peripheral Nerve Sheath Tumor",
]
MPNST_ids <- MPNST_ids[!duplicated(MPNST_ids$Sample.ID),]

if (!file.exists(paste0(raw_dir, "TCGA_MPNST/TCGA_SARC_450K_hg19.rda"))) {

  # download all sarcoma data from TGCA:
  TCGAbiolinks:::getGDCprojects()$project_id
  
  cohort_name <- "TCGA-SARC"
  
  query.met <- GDCquery(
    project = cohort_name,
    data.category = "DNA methylation", 
    platform = "Illumina Human Methylation 450", 
    sample.type = "Primary Tumor",
    legacy = TRUE
  )
  
  GDCdownload(query.met)
    
  data <- GDCprepare(
    query = query.met,
    save = TRUE,
    save.filename = paste0(raw_dir, "TCGA_MPNST/TCGA_SARC_450K_hg19.rda"),
    summarizedExperiment = TRUE
  )

} else {

  load(paste0(raw_dir, "TCGA_MPNST/TCGA_SARC_450K_hg19.rda"))

}

# isolate beta values:
sarc_beta <- data@assays@data@listData[[1]]

# isolate MPNST samples:
for (i in 1:nrow(MPNST_ids)) {

  print(paste0("Fetching ", MPNST_ids$Sample.ID[i], "..."))

  if (i==1) {

    TCGA_MPNST_beta <- data.frame(
      sarc_beta[
        , grep(MPNST_ids$Sample.ID[i], colnames(sarc_beta))
      ]
    )
    colnames(TCGA_MPNST_beta)[i] <- MPNST_ids$Sample.ID[i]

  } else {

    TCGA_MPNST_beta <- cbind(
      TCGA_MPNST_beta,
      data.frame(
        sarc_beta[
          , grep(MPNST_ids$Sample.ID[i], colnames(sarc_beta))
        ]
      )
    )
    colnames(TCGA_MPNST_beta)[i] <- MPNST_ids$Sample.ID[i]

  }

}

MPNST_beta$TGCA <- TCGA_MPNST_beta

# isolate common probes from both datasets:
jour_450k <- MPNST_beta$jour[
  rownames(MPNST_beta$jour) %in% rownames(MPNST_beta$TGCA),
]
TGCA_450k <- MPNST_beta$TGCA[
  rownames(MPNST_beta$TGCA) %in% rownames(MPNST_beta$jour),
]

MPNST_beta_df <- cbind(jour_450k, TGCA_450k)

saveRDS(MPNST_beta_df, paste0(Robject_dir, "MPNST_beta_df.Rdata"))


####################################################################################
### 3. Download non-malignant TGCA data
####################################################################################

normTGCA_files <- list.files(paste0(raw_dir, "TCGA_normal"))

for (i in 1:length(normTGCA_files)) {

  print(paste0("Loading file ", i, "; ", normTGCA_files[i]))

  load(paste0(raw_dir, "TCGA_normal/", normTGCA_files[i]))

  if (i==1) {
    normal_beta <- list(data@assays@data@listData[[1]])
    names(normal_beta) <- gsub("_450K.*$", "", normTGCA_files[i])
  } else {
    normal_beta[[i]] <- data@assays@data@listData[[1]]
    names(normal_beta) <- gsub("_450K.*$", "", normTGCA_files[i])
  }

}

saveRDS(normal_beta, paste0(Robject_dir, "normal_beta.Rdata"))


####################################################################################
### 4. Download non-malignant Moss 2018 data
####################################################################################

if (!(file.exists(paste0(Robject_dir, "all_beta_3.Rdata")))) {

  moss_soft <- getGEO("GSE122126", GSEMatrix=T)

  # fetch phenotype data:
  moss_pheno <- list(
    inf450k = list(
      phenotype = subset(
        pData(moss_soft[["GSE122126-GPL13534_series_matrix.txt.gz"]]), 
        select = c(
          "title", "geo_accession", "type", "characteristics_ch1", "label_ch1", 
          "data_row_count"
        )
      ),
      protocols = subset(
        pData(moss_soft[["GSE122126-GPL13534_series_matrix.txt.gz"]]), 
        select = c(
          "title", "geo_accession", "extract_protocol_ch1", "label_protocol_ch1", 
          "hyb_protocol", "scan_protocol", "data_processing"
        )
      )
    ),
    epic = list(
      phenotype = subset(
        pData(moss_soft[["GSE122126-GPL21145_series_matrix.txt.gz"]]), 
        select = c(
          "title", "geo_accession", "type", "characteristics_ch1", "label_ch1", 
          "data_row_count"
        )
      ),
      protocols = subset(
        pData(moss_soft[["GSE122126-GPL21145_series_matrix.txt.gz"]]), 
        select = c(
          "title", "geo_accession", "extract_protocol_ch1", "label_protocol_ch1", 
          "hyb_protocol", "scan_protocol", "data_processing"
        )
      )
    )
  )
  
  
  # fetch beta value matrix:
  moss_beta <- list(
    inf450k = moss_soft[["GSE122126-GPL13534_series_matrix.txt.gz"]]@assayData$exprs,
    epic = moss_soft[["GSE122126-GPL21145_series_matrix.txt.gz"]]@assayData$exprs
  )

  # keep only 450k probes:
  moss_beta_450k <- cbind(
    moss_beta$inf450k,
    moss_beta$epic[rownames(moss_beta$epic) %in% rownames(moss_beta$inf450k),]
  )

  # isolate MPNST samples only:
  moss_pheno <- lapply(moss_pheno, function(x) {
    x <- x[grep("MPNST", x$title),]
    x$title <- gsub("^.* ", "", x$title)
    return(x)
  })
  
  moss_beta <- moss_beta[,colnames(moss_beta) %in% moss_pheno$phenotype$geo_accession]
  
  # initiate list for all data:
  all_beta <- list(
    moss = moss_beta
  )
  saveRDS(all_beta, paste0(Robject_dir, "all_beta_3.Rdata"))
  
  all_pheno <- list(
    moss_pheno$phenotype
  )
  saveRDS(all_pheno, paste0(Robject_dir, "all_pheno_3.Rdata"))
  
  all_protocols <- list(
    moss_pheno$protocols
  )
  saveRDS(all_protocols, paste0(Robject_dir, "all_protocols_3.Rdata"))

} else {

  all_beta <- readRDS(paste0(Robject_dir, "all_beta_3.Rdata"))
  
  all_pheno <- readRDS(paste0(Robject_dir, "all_pheno_3.Rdata"))
  
  all_protocols <- readRDS(paste0(Robject_dir, "all_protocols_3.Rdata"))

}


####################################################################################
### 4. Download TCGA non-malignant data
####################################################################################

# check available datasets:
checkTCGA('Dates')

# check available cohorts:
(cohorts <- infoTCGA() %>% 
   rownames() %>% 
   sub("-counts", "", x=.))

# check available datasets:
avail_ds <- checkTCGA("DataSets", "SARC", date = "2016-01-28")
met_ds <- avail_ds[grep("methyl", avail_ds$Name, ignore.case = T),]

# specify dataset details:
releaseDate <- "2016-01-28"
cohort <- "SARC"
datasets <- met_ds$Name

# download data:
all_met_files <- list.files(
  raw_dir, pattern = "methyl", recursive = TRUE, full.names = TRUE
)

if (length(all_met_files) == 0) {
  for (dataset in datasets) {
    try(downloadTCGA(
      cancerTypes = cohort, 
      destDir = raw_dir, 
      date = releaseDate, 
      dataSet = dataset
      ),
      silent=TRUE
    )
  }
}

# read and format data:
all_met_files <- list.files(
  raw_dir, pattern = "met.*txt", recursive = TRUE, full.names = TRUE
)


for (i in 1:length(all_met_files)) {

  print(paste0("Adding methylation values from ", all_met_files[i]))

  # load met_vals:
  met_values <- read.table(all_met_files[i], fill = T, header = T)

  ######
  met_sub <- met_values[1:100,1:100]

  # merge 2 colnames rows:
  colnames(met_sub) <- paste0(colnames(met_sub), met_sub[1,])
  met_sub <- met_sub[-1,]

  # change chromosome 23 and 24 to X and Y:
  met_sub$Chromosome[met_sub$Chromosome == "23"] <- "X"

  # convert to granges object:
  met_gr <- GRanges(
    seqnames = Rle(paste0("chr", met_sub$Chromosome)),
    ranges = IRanges(
      start = met_sub$Start, 
      end = met_sub$End
    ),
    strand = Rle("*")
  )

  ######

  # merge 2 colnames rows:
  colnames(met_values) <- paste0(colnames(met_values), met_values[1,])
  met_values <- met_values[-1,]

  # change chromosome 23 and 24 to X and Y:
  met_values$Chromosome[met_values$Chromosome == "23"] <- "X"

  # convert to granges object:
  met_gr <- GRanges(
    seqnames = Rle(paste0("chr", met_values$Chromosome)),
    ranges = IRanges(
      start = met_values$Start, 
      end = met_values$End
    ),
    strand = Rle("*")
  )

  print(paste0("Length of granges object is: ", length(CNA_gr)))

  if ( length(grep("hg18", all_cna_files[i])) > 0 ) {

    new_filename <- gsub("(.*)hg18.*", "\\1hg38.Rdata", all_cna_files[i])
    
    if (!file.exists(new_filename)) {

      print("Co-ordinates mapped to hg18, converting to hg38...")

      # liftover to hg38:
      CNA_gr_hg38 <- do_liftover(CNA_gr, "hg18", "hg38")
  
      # for each element, keep only min and max coordinates of main chromosome:
      CNA_gr_hg38 <- lapply(CNA_gr_hg38, collapse_lifted_CNA)
      CNA_gr_hg38 <- do.call("c", CNA_gr_hg38)
  
      saveRDS(CNA_gr_hg38, new_filename)

    } else {

      print("Loading hg38 version and adding...")
      CNA_gr_hg38 <- readRDS(new_filename)

    }

    CNA_gr <- CNA_gr_hg38

  } else if ( length(grep("hg19", all_cna_files[i])) > 0 ) {

    new_filename <- gsub("(.*)hg19.*", "\\1hg38.Rdata", all_cna_files[i])
    
    if (!file.exists(new_filename)) {

      print("Co-ordinates mapped to hg18, converting to hg38...")
  
      # liftover to hg38:
      CNA_gr_hg38 <- do_liftover(CNA_gr, "hg19", "hg38")
  
      # for each element, keep only min and max coordinates of main chromosome:
      CNA_gr_hg38 <- lapply(CNA_gr_hg38, collapse_lifted_CNA)
      CNA_gr_hg38 <- do.call("c", CNA_gr_hg38)
  
      saveRDS(CNA_gr_hg38, new_filename)

    } else {

      print("Loading hg38 version...")
      CNA_gr_hg38 <- readRDS(new_filename)

    }

    CNA_gr <- CNA_gr_hg38

  }

  if (c==1) {
    all_CNA <- CNA_gr
  } else {
    all_CNA <- c(all_CNA, CNA_gr)
  }

}


