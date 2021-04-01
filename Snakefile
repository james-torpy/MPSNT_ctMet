# Run command:
# snakemake --reason --cores 100 --cluster-config cluster.json --cluster 'qsub -pe smp {cluster.cores} -N DMR.smk -wd '/share/ScratchGeneral/jamtor/projects/MPNST_ctMet/logs' -b y -j y -V -P DSGClinicalGenomics' -j 10

# DAG command:
# snakemake --dag | dot -Tsvg > dag.svg

# define/create directories:
project_name = "MPNST_ctMet"
home_dir = '/share/ScratchGeneral/jamtor/'
project_dir = home_dir + 'projects/' + project_name + '/'
script_dir = project_dir + 'scripts/'

# define variables:
BETA_VALS = list(['0.35', '0.4'])
PVALS = list(['1e-5'])
BG_CUTOFFS = list(['0.2'])
BLOOD_HYPOS = list(['0.2'])
BLOOD_HYPERS = list(['0.8'])
NF_HYPOS = list(['0.2'])
NF_HYPERS = list(['0.8'])

rule all:
    input:
        expand(
            'results/beta_{beta_val}/pval_{pval}/background_{bg_cutoff}/' + 
                'blood_hypo_{blood_hypo}/blood_hyper_{blood_hyper}/' +
                'NF_hypo_{NF_hypo}/NF_hyper_{NF_hyper}/plots/' + 
                'MPNST_hypomethylated_marker_probes.png',
            beta_val = BETA_VALS, pval = PVALS,
            bg_cutoff = BG_CUTOFFS, blood_hypo = BLOOD_HYPOS, 
            blood_hyper = BLOOD_HYPERS, NF_hypo = NF_HYPOS, 
            NF_hyper = NF_HYPERS
        )

rule find_DMR:
    output:
        'results/beta_{beta_val}/pval_{pval}/background_{bg_cutoff}/' + 
                'blood_hypo_{blood_hypo}/blood_hyper_{blood_hyper}/' +
                'NF_hypo_{NF_hypo}/NF_hyper_{NF_hyper}/plots/' + 
                'MPNST_hypomethylated_marker_probes.png',
    shell:
        "mkdir -p logs/find_DMR/beta_{wildcards.beta_val}/pval_{wildcards.pval}/background_{wildcards.bg_cutoff}/" + 
            "blood_hypo_{wildcards.blood_hypo}/blood_hyper_{wildcards.blood_hyper}/NF_hypo_{wildcards.NF_hypo}/" + 
            "NF_hyper_{wildcards.NF_hyper}/; " + 
        "cd logs/find_DMR/beta_{wildcards.beta_val}/pval_{wildcards.pval}/background_{wildcards.bg_cutoff}/" + 
            "blood_hypo_{wildcards.blood_hypo}/blood_hyper_{wildcards.blood_hyper}/NF_hypo_{wildcards.NF_hypo}/" + 
            "NF_hyper_{wildcards.NF_hyper}/; " + 
        "R CMD BATCH  --no-save '--args " + 
        "{wildcards.beta_val} " +
        "{wildcards.pval} " +
        "{wildcards.bg_cutoff} " +
        "{wildcards.blood_hypo} " +
        "{wildcards.blood_hyper} " + 
        "{wildcards.NF_hypo} " +
        "{wildcards.NF_hyper}" + 
        "FALSE " + # check individual tissues for DMR
        "10 " + # min samples to include tissue in plot
        "' ../../../../../../../../../scripts/1.find_DMR.R 2> find_DMR.sort.errors"




