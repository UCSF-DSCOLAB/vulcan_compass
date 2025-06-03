configfile: "config.yaml"

rule all:
    input:
        thumbnail="output/red_blue.png",
        compass_tgz="output/compass.tar.gz"

rule get_dataset_and_summarize:
    params:
        dataset_name=config["dataset_name"]
    output:
        scdata="output/scdata.Rds",
        plotting_options="output/scdata_plotting_options.json",
        discrete_metadata_summary="output/scdata_discrete_metadata_summary.json",
        all_opts="output/scdata_all_opts.json",
        continuous_opts="output/scdata_continuous_opts.json",
        discrete_opts="output/scdata_discrete_opts.json",
        reduction_opts="output/scdata_reduction_opts.json"
    resources:
        mem_mb: 100000 # 100gbs
        runtime: "2h"
    singularity:
        "/dscolab/vulcan/containers/archimedes-r.sif"
    script:
        "scripts/get_dataset_and_summarize.R"
rule pseudobulk_dataset:
    params:
        min_cells=config["min_cells"],
        sample_id_column=config["sample_id_column"],
        cell_type_column=config["cell_type_column"],
        norm_method=config["norm_method"]
    input:
        target_genes="module_compass_targets/genes_targeted.txt",
        scd_rds="output/scdata.Rds"
    output:
        pseudo_matrix="output/delog_pseudobulk_matrix.tsv",
        pseudo_metadata="output/pseudobulk_metadata.tsv",
        cell_types="output/pseudobulk_celltypes.json",
        samples="output/pseudobulk_samples.json",
        discrete_metadata_summary="output/pseudobulk_discrete_metadata_summary.json"
    resources:
        mem_mb: 100000 # 100gbs
        runtime: "2h"
    singularity:
        "/dscolab/vulcan/containers/archimedes-r.sif"
    script:
        "scripts/pseudobulk_data_elements.R"

rule run_compass:
    params:
        species=config["species"]
    input:
        pseudo_matrix="output/delog_pseudobulk_matrix.tsv",
        meta_subsystems="module_compass_targets/meta_subsystems.txt"
    threads: 30
    resources:
        mem_mb: 150000 # 150gbs
        runtime: "7d"
    output:
        carbon_reactions="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/reactions.tsv",
        lipid_reactions="output/compass_output/LIPID_META_SUBSYSTEM/reactions.tsv",
        AA_reactions="output/compass_output/AA_META_SUBSYSTEM/reactions.tsv"
    singularity:
        "/dscolab/vulcan/containers/compass.sif"
    shell:
        """
        if [ "{params.species}" = "human" ]; then
            species="homo_sapiens"
        elif [ "{params.species}" = "mouse" ]; then
            species="mus_musculus"
        fi

        compass \
            --num-processes {threads} \
            --temp-dir output/compass_tmp \
            --output-dir output/compass_output \
            --data {input.pseudo_matrix} \
            --species ${{species}} \
            --model Human1 \
            --select-meta-subsystems {input.meta_subsystems}
        """
rule after_compass:
    input:
        carbon_reactions="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/reactions.tsv",
        lipid_reactions="output/compass_output/LIPID_META_SUBSYSTEM/reactions.tsv",
        AA_reactions="output/compass_output/AA_META_SUBSYSTEM/reactions.tsv"
    output:
        carbon_reaction_scores="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/reaction_scores.tsv",
        carbon_reaction_norm_sum="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/reactions_norm_sum.tsv",
        carbon_reaction_norm_rank="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/reactions_norm_rank.tsv",
        carbon_readme="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/readme.tsv",
        lipid_reaction_scores="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/reaction_scores.tsv",
        lipid_reaction_norm_sum="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/reactions_norm_sum.tsv",
        lipid_reaction_norm_rank="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/reactions_norm_rank.tsv",
        lipid_readme="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/readme.tsv",
        AA_reaction_scores="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/reaction_scores.tsv",
        AA_reaction_norm_sum="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/reactions_norm_sum.tsv",
        AA_reaction_norm_rank="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/reactions_norm_rank.tsv",
        AA_readme="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/readme.tsv",
        compass_tgz="output/compass.tar.gz"
    singularity:
        "/dscolab/vulcan/containers/archimedes-py.sif"
    script:
        "scripts/normalize_compass_reactions_and_compress.py"

rule ui_diff_targets_cells: #UI
    input:
        "output/pseudobulk_discrete_metadata_summary.json",
        "output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/reactions.tsv"
    params:
        ui=True
    output:
        ["output/diff_group_1__formula.json", "output/diff_group_2__formula.json"]
rule ui_diff_targets_subsystem: #UI
    input:
        "output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/reactions.tsv"
    params:
        ui=True
    output:
        "output/diff_subsystem_target.txt"

rule parse_groupings:
    params:
        sample_id_column=config["sample_id_column"],
        cell_type_column=config["cell_type_column"]
    input:
        group_def_1="output/diff_group_1__formula.json",
        group_def_2="output/diff_group_2__formula.json",
        pseudo_metadata="output/pseudobulk_metadata.tsv"
    output:
        group_1_inds="output/diff_group_1__indexes.json",
        group_2_inds="output/diff_group_2__indexes.json",
        group_1_metadata="output/diff_group_1__metadata.csv",
        group_2_metadata="output/diff_group_2__metadata.csv",
        groups_summary="output/diff_groups_summary.md"
    singularity:
        "/dscolab/vulcan/containers/archimedes-py.sif"
    script:
        "scripts/interpret_diff_groups.py"

rule plot_red_blue:
    params:
        norm_method=config["norm_method"]
    input:
        subsystem="output/diff_subsystem_target.txt",
        group_1_inds="output/diff_group_1__indexes.json",
        group_2_inds="output/diff_group_2__indexes.json",
        carbon_reaction_scores="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/reaction_scores.tsv",
        carbon_reaction_norm_sum="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/reactions_norm_sum.tsv",
        carbon_reaction_norm_rank="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/reactions_norm_rank.tsv",
        lipid_reaction_scores="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/reaction_scores.tsv",
        lipid_reaction_norm_sum="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/reactions_norm_sum.tsv",
        lipid_reaction_norm_rank="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/reactions_norm_rank.tsv",
        AA_reaction_scores="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/reaction_scores.tsv",
        AA_reaction_norm_sum="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/reactions_norm_sum.tsv",
        AA_reaction_norm_rank="output/compass_output/CENTRAL_CARBON_META_SUBSYSTEM/reactions_norm_rank.tsv"
    output:
        plot="output/red_blue.png"
    singularity:
        "/dscolab/vulcan/containers/compass.sif"
    script:
        "scripts/plot_red_blue.py"
