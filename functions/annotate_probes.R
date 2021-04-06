annotate_probes <- function(dmr) {

  if (!is.na(dmr)) {

    # make probe column:
    dmr$probe <- rownames(dmr)
    # fetch probe locations and add:
    specific_coords <- probe_coords[
      names(probe_coords) %in% dmr$probe
    ]
    return(
      merge(
        dmr,
        data.frame(
          probe = names(specific_coords),
          chr = seqnames(specific_coords),
          coord = start(specific_coords)
        ),
        by = "probe"
      )
    )

  } else {
    return(dmr)
  }

}