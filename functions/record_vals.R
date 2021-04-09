record_vals <- function(
  dmr, 
  label, 
  col_names
) {

  # bind hyper and hypo candidates together:
  if ("hyper" %in% names(dmr)) {

    if ("all_non_malig" %in% names(dmr$hyper)) {
      record_obj <- rbind(
        nrow(dmr$hyper$all_non_malig),
        nrow(dmr$hypo$all_non_malig)
      )
    } else {
      record_obj <- rbind(
        nrow(dmr$hyper),
        nrow(dmr$hypo)
      )
    }

    # fill in other tissue values with NA:
    record_obj <- cbind(
      data.frame(matrix(NA, nrow = 2, ncol = length(col_names)-1)),
      record_obj
    )
    colnames(record_obj) <- col_names

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

