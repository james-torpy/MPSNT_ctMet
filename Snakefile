# Run command:
# snakemake --reason --cores 100 --cluster-config cluster.json --cluster 'qsub -pe smp {cluster.cores} -N DMR.smk -wd '/share/ScratchGeneral/jamtor/projects/MPNST_ctMet/logs' -b y -j y -V -P DSGClinicalGenomics' -j 10

# DAG command:
# snakemake --dag | dot -Tsvg > dag.svg

# define/create directories:
project_name = "MPNST_ctMet"
home_dir = '/share/ScratchGeneral/jamtor/'
project_dir = home_dir + 'projects/' + project_name + '/'
script_dir = project_dir + 'scripts/'

# deine variables:
#BETA_VALS = list(['0.35'])
#PVALS = list(['1e-4'])
#BG_CUTOFFS = list(['0.1'])
#BLOOD_HYPOS = list(['0.1'])
#BLOOD_HYPERS = list(['0.8'])
BETA_VALS = list(['0.3', '0.35', '0.4'])
PVALS = list(['1e-5', '1e-4'])
BG_CUTOFFS = list(['0.1', '0.2'])
BLOOD_HYPOS = list(['0.1', '0.2'])
BLOOD_HYPERS = list(['0.8', '0.9'])

#rule all:
#    input:
#        expand(
#            'results/beta_{beta_val}/DMR_counts.txt',
#            beta_val = BETA_VALS
#        )
#
#rule mkdir:
#    output:
#        'results/beta_{beta_val}/DMR_counts.txt'
#    shell:
#        "mkdir -p logs/find_DMR; " + 
#        "cd logs/find_DMR; " + 
##        "touch ../../{output}"
#        "R CMD BATCH  --no-save '--args " + 
#        "{wildcards.beta_val} " + 
#        "' ../../scripts/find_DMR.R 2> find_DMR.sort.errors"

rule all:
    input:
        expand(
            'results/beta_{beta_val}/pval_{pval}/background_{bg_cutoff}/' + 
                'blood_hypo_{blood_hypo}/blood_hyper_{blood_hyper}/tables/' + 
                'DMR_counts.txt',
            beta_val = BETA_VALS, pval = PVALS,
            bg_cutoff = BG_CUTOFFS, blood_hypo = BLOOD_HYPOS, 
            blood_hyper = BLOOD_HYPERS
        )

rule find_DMR:
    output:
        'results/beta_{beta_val}/pval_{pval}/background_{bg_cutoff}/' + 
            'blood_hypo_{blood_hypo}/blood_hyper_{blood_hyper}/tables/' + 
            'DMR_counts.txt'
    shell:
        "mkdir -p logs/find_DMR/beta_{wildcards.beta_val}/pval_{wildcards.pval}/background_{wildcards.bg_cutoff}/" + 
            "blood_hypo_{wildcards.blood_hypo}/blood_hyper_{wildcards.blood_hyper}/; " + 
        "cd logs/find_DMR/beta_{wildcards.beta_val}/pval_{wildcards.pval}/background_{wildcards.bg_cutoff}/" + 
            "blood_hypo_{wildcards.blood_hypo}/blood_hyper_{wildcards.blood_hyper}/; " + 
        "R CMD BATCH  --no-save '--args " + 
        "{wildcards.beta_val} " +
        "{wildcards.pval} " +
        "{wildcards.bg_cutoff} " +
        "{wildcards.blood_hypo} " +
        "{wildcards.blood_hyper}" + 
        "' ../../../../../../../scripts/find_DMR.R 2> find_DMR.sort.errors"




