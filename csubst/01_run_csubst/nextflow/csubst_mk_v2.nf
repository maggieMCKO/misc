#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

params.conda_env_path = "/usr/users/ksadana/anaconda3/envs/env_csubst39_14"

params.foreground = "/usr/users/ksadana/Genome/CSUBST/genome-wide/00_inputs/foreground.tsv"
params.tree_file="/usr/users/ksadana/Genome/CSUBST/genome-wide/00_inputs/subset_tree.nwk"
params.sp_ls="/usr/users/ksadana/Genome/CSUBST/genome-wide/00_inputs/sp_list.txt"
params.aln_input="/usr/users/ksadana/Genome/CSUBST/genome-wide/00_inputs/nt_align_renamed"
params.transcript = "/usr/users/ksadana/Genome/CSUBST/genome-wide/00_inputs/checked.transcripts.lst"
// params.transcript = "$HOME/Nectar/analysis/csubst/01_run_csubst_2024_proto/sunbrid_affected.txt"

params.prune_tree_script = "$projectDir/prune_ete4.py"
params.outdir = "results"

// to compile results
params.combineRscript = "$projectDir/concatenate_results.R"
params.cusbstdir = "results/csubst_out/" 
params.combinedoutdir = "results/combined/"


params.size = 1000

workflow {

    query_ch = Channel.fromPath( params.transcript ).splitText() { it.strip() }.buffer(size: params.size, remainder: true)
    
    // run 1
    ch1 = cp_aln( query_ch ) 
    ch1.flatten().set{ ch_aln }

    // run 2
    ch2 = prune_tree( ch1 ) 
    ch2.flatten().set{ ch_tre }

    // prep for 3
    ch_aln
        .map { [it.toString().split("/").last().split(".fa")[0], it] }
        .set { ch_aln }
    ch_tre
        .map { [it.toString().split("/").last().split("pruned_")[1].split(".nw")[0], it] }
        .set { ch_tre }
    // ch_aln.view()
    // ch_tre.view()
    
    // run 3
    ch_aln.combine(ch_tre, by:0) | runcsubst

    // after csubst finished running, compile results
    query_ch = Channel.of( 'csubst_b.tsv', 'csubst_cb_3.tsv', 'csubst_cb_4.tsv', 'csubst_cb_5.tsv', 'csubst_cb_6.tsv', 
    'csubst_cb_7.tsv', 'csubst_cb_8.tsv', 'csubst_cb_9.tsv', 'csubst_cb_10.tsv', 'csubst_cb_stats.tsv')
    combineData( query_ch ) 

    // 'csubst_cb_2.tsv'
    query_ch = Channel.fromPath( params.transcript ).splitText() { it.strip() }.buffer(size: params.size, remainder: true)
    combineData_for2( query_ch ).flatten() \
    | collectFile(name:'results/combined/combined_csubst_cb_2.tsv', keepHeader: true, storeDir: params.combinedoutdir) \
    | view{it}

}

process cp_aln {
    tag "cp_aln $transcript"
    publishDir params.outdir, mode : 'copy' 

    conda "$conda_env_path"
    scratch true
    errorStrategy 'ignore'

    input:
        val id
    
    output:
        // tuple val ("${INDEX}"), path("aln/rm_HLcinPul1_${INDEX}.fa"), emit: aln // doesnt work
        path("aln/*.fa"), emit: aln

    """
    for INDEX in ${id.join(' ')}; do
        echo "task \${INDEX}"
        mkdir -p aln
        cp "${params.aln_input}/\${INDEX}.fasta" "aln/\${INDEX}.fa"
    done
    """
}

process prune_tree {
    tag "prune_tree $transcript"
    publishDir params.outdir, mode : 'copy' 

    conda "$conda_env_path"
    scratch true
    errorStrategy 'ignore'

    input:
        val alns
    
    output:
        // path "tmp.txt", emit: pruned_tree // for test
        path ("tree/pruned_*.nw"), emit: pruned_tree
        // tuple val (stem), path ("tree/pruned_*.nw"), emit:pruned_tree // doesnt work

    """    
    for aln in ${alns.join(' ')}; do
        stem=\$(basename \${aln} )
        stem=\${stem%%.fa}
        mkdir -p tree
        python ${params.prune_tree_script} ${params.tree_file} \$aln "tree/pruned_\${stem}.nw"
    done
    """
}

process runcsubst {
    tag "csubst $transcript"
    publishDir params.outdir, mode : 'copy',
    saveAs: {
        filename -> "csubst_out/$id/$filename"
    }

    conda "$conda_env_path"
    scratch true
    errorStrategy 'ignore'

    input:
        tuple val(id), path(aln), path(tree)
    
    output:
        path "csubst_*"
        
    """
    csubst analyze \
    --alignment_file $aln \
    --rooted_tree_file $tree \
    --foreground ${params.foreground} \
    --fg_exclude_wg no \
    --max_arity 10 \
    --exhaustive_until 2 &> csubst_${id}.log
    """
}


process combineData {

    publishDir params.outdir, mode : 'copy' 

    conda "$conda_env_path"
    scratch true
    errorStrategy 'ignore'

    input:
        val id
    
    // output:
    //     path("combined/combined_$id")

    """
    mkdir -p combined
    Rscript ${params.combineRscript} ${params.cusbstdir} ${id} 8
    """
}

process combineData_for2 {

    // publishDir params.outdir, mode : 'copy' 

    conda "$conda_env_path"
    scratch true
    errorStrategy 'ignore'

    input:
        val id
    
    output:
        path 'csubst_cb_2a_*.tsv'

    """
    for INDEX in ${id.join(' ')}; do
        echo "task \${INDEX}"
        awk -v var="\$INDEX" 'BEGIN{FS=OFS="\t"} {print (NR>1?var:"transcript"), \$0}' ${params.cusbstdir}\${INDEX}/csubst_cb_2.tsv > csubst_cb_2a_\${INDEX}.tsv
    done
    """
}