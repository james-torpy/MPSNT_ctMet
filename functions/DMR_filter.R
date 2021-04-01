DMR_filter <- function(
  dmr,
  beta_cutoff,
  pval_cutoff
) {

  # filter probes by methylation difference and p-values:
  filt <- dmr[
    dmr[grep("p.value.adj", colnames(dmr))] < pval_cutoff &
    (dmr[grep("minus", colnames(dmr))] >= beta_cutoff |
    dmr[grep("minus", colnames(dmr))] <= -beta_cutoff),
  ]
  
  # identify probes either hyper- or hypo-methylated in MPNST vs normal:
  return(
  	list(
      hyper = filt[
        grep("Hypermethylated", filt$status),
      ],
      hypo = filt[
        grep("Hypomethylated", filt$status),
      ]
    )
  )

}