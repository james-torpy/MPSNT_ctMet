filter_by_NF <- function(
  dmr, 
  NF_vals, 
  gene_annot
) {

  for (i in 1:length(dmr)) {

    if (nrow(dmr[[i]]) > 0) {
  
      # convert DM probe df to granges object:
      dm_gr <- GRanges(
        seqnames = dmr[[i]]$chr,
        ranges = IRanges(
          start = dmr[[i]]$coord,
          end = dmr[[i]]$coord
        ),
        strand = "*"
      )
      values(dm_gr) <- subset(
        dmr[[i]], 
        select = -c(chr, coord)
      )
  
      # find overlaps with genes and keep in separate gr:
      gene_olaps <- findOverlaps(gene_annot, dm_gr)
  
      if (i==1) {
        DMR_genes <- list(
          hyper = gene_annot[queryHits(gene_olaps)]
        )
        DMR_genes$hyper$probe <- dm_gr$probe[subjectHits(gene_olaps)]
      } else {
        DMR_genes$hypo <- gene_annot[queryHits(gene_olaps)]
        DMR_genes$hypo$probe <- dm_gr$probe[subjectHits(gene_olaps)]
      }
    
      if (length(grep("Hyper", dmr[[i]]$status[i])) > 0) {
    
        # isolate hypomethylated NF regions:
        NF_spec <- NF_vals[
          NF_vals$score < NF_cutoffs$MPNST_hypermethylated,
        ]
    
        # find overlaps with gr:
        olaps <- findOverlaps(dm_gr, NF_spec)
    
      } else {
    
        # isolate hypermethylated NF regions:
        NF_spec <- NF_vals[
          NF_vals$score > NF_cutoffs$MPNST_hypomethylated,
        ]
    
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
        final_DMR <- list(DMR_res)
      } else {
        final_DMR[[i]] <- DMR_res
      }
  
    } else {
  
      if (i==1) {
        final_DMR <- list(NULL)
      } else {
        final_DMR[[i]] <- NULL
      }
  
    }
    
  }
  names(final_DMR) <- names(dmr)

  return(
    list(
      final_DMR = final_DMR, 
      DMR_genes = DMR_genes
    )
  )

}
