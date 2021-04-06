indiv_tissue_filter <- function(dmr) {

  # separate hyper and hypo dfs:
  DMR_list <- list(
    hyper = lapply(dmr, function(x) x$hyper),
    hypo = lapply(dmr, function(x) x$hypo)
  )
  
  # find DMR probes common to all normal tissue:
  DMR_probes <- lapply(DMR_list, function(x) {
    lapply(x, rownames)
  })
  DMR_temp_probes <- lapply(DMR_probes, function(x) {
  
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
  for (i in 1:length(DMR_list)) {
  
    if (length(DMR_list[[i]]) > 0) {
  
     if (!exists("DMR_temp")) {
  
       DMR_temp <- list(
         lapply(DMR_list[[i]], function(y) {
            y[rownames(y) %in% DMR_temp_probes[[i]],]
          })
        )
        names(DMR_temp) <- names(DMR_list)[i]
  
     } else {
  
       DMR_temp[[i]] <- lapply(DMR_list[[i]], function(y) {
          y[rownames(y) %in% DMR_temp_probes[[i]],]
        })
        names(DMR_temp)[i] <- names(DMR_list)[i]
  
     }
     
    }
    
  }

  # turn list inside out:
  for (i in 1:length(DMR_temp)) {

    if (i==1) {
      DMR_indiv <- list(
        list(
          hyper = DMR_temp$hyper[[i]],
          hypo = DMR_temp$hypo[[i]]
        )
      )
    } else {
      DMR_indiv[[i]] <- list(
        hyper = DMR_temp$hyper[[i]],
        hypo = DMR_temp$hypo[[i]]
      )
    }

  }
  names(DMR_indiv) <- names(dmr)

  return(DMR_indiv)

}

