record_vals <- function(dmr, label) {

  if ("hyper" %in% names(dmr)) {

    record_obj <- rbind(
      nrow(dmr$hyper$all_non_malig),
      nrow(dmr$hypo$all_non_malig)
    )

    # fill in other tissue values with NA:
    record_obj <- cbind(
      data.frame(matrix(NA, nrow = 2, ncol = length(dmr$hyper)-1)),
      record_obj
    )
    colnames(record_obj) <- names(dmr$hyper)

  } else {

    record_obj <- as.data.frame(
      rbind(
        sapply(dmr, function(x) {
          nrow(
            x$hyper
          )
        }), 
        sapply(dmr, function(x) {
          nrow(
            x$hypo
          )
        })
      )
    )

  }
  
  rownames(record_obj) <- c(
    paste0(label, "_hyper"), 
    paste0(label, "_hypo")
  )

  print(
    paste0(
      record_obj[1,]$all_non_malig, 
      " hypermethylated and ",
      record_obj[2,]$all_non_malig,
      " hypomethylated ", label, " candiates ", 
      "identified..."
    )
  )

  return(record_obj)

}

