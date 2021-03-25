DMR_filter <- function(
  dmr,
  beta_cutoff,
  pval_cutoff
) {

  # filter probes by methylation difference and p-values:
  filt <- dmr[
    dmr[grep("p.value.adj", colnames(dmr))] < 10e-5 &
    (dmr[grep("minus", colnames(dmr))] >= 0.35 |
    dmr[grep("minus", colnames(dmr))] <= -0.35),
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