configfile: "config.yaml"

rule all:
    input:
        thumbnail="output/red_blue.png",
        red_blue_plot="output/red_blue.png"

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
        scd_rds="output/scdata.Rds"
    output:
        pseudo_matrix="output/delog_pseudobulk_matrix.tsv",
        pseudo_metadata="output/pseudobulk_metadata.tsv",
        cell_types="output/pseudobulk_celltypes.json",
        samples="output/pseudobulk_samples.json",
        discrete_metadata_summary="output/pseudobulk_discrete_metadata_summary.json"
    singularity:
        "/dscolab/vulcan/containers/archimedes-r.sif"
    script:
        "scripts/pseudobulk_data_elements.R"

rule run_compass:
    params:
        species=config["species"],
        norm_method=config["norm_method"]
    input:
        pseudo_matrix="output/delog_pseudobulk_matrix.tsv",
        meta_subsystems="module_compass_targets/meta_subsystems.txt"
    output:
        carbon_reactions="output/compass_outputs/CENTRAL_CARBON_META_SUBSYSTEM/reactions.tsv",
        lipid_reactions="output/compass_outputs/LIPID_META_SUBSYSTEM/reactions.tsv",
        AA_reactions="output/compass_outputs/AA_META_SUBSYSTEM/reactions.tsv"
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
            --num-processes 30 \
            --temp-dir output/compass_outputs/tmp \
            --output-dir output/compass_outputs \
            --data {input.pseudo_matrix} \
            --species ${species} \
            --model Human1 \
            --select-meta-subsystems {input.meta_subsystems}
        
        if [ "{params.norm_method}" = "late__by_reaction_sum" ]; then
            echo "NORMALIZATION NOT YET IMPLEMENTED"
        elif [ "{params.norm_method}" = "late__by_reaction_rank" ]; then
            echo "NORMALIZATION NOT YET IMPLEMENTED"
        fi
        """

rule ui_diff_targets_cells: #UI
    input: [
        "output/pseudobulk_discrete_metadata_summary.json",
        "output/compass_outputs/CENTRAL_CARBON_META_SUBSYSTEM/reactions.tsv"
    ]
    output:
        "output/diff_groups_selections.json"
rule ui_diff_targets_subsystem: #UI
    input:
        "output/compass_outputs/CENTRAL_CARBON_META_SUBSYSTEM/reactions.tsv"
    output:
        "output/diff_subsystem_target.txt"
rule plot_red_blue:
    input:
        subsystem="output/diff_subsystem_target.txt",
        group_defs="output/diff_groups_selections.json",
        pseudo_metadata="output/pseudobulk_metadata.tsv"
    output:
        plot="outputs/red_blue.png"
    singularity:
        "/dscolab/vulcan/containers/compass.sif"
    script:
        "scripts/plot_red_blue.py"

# ### Optional standard scViz plotting
# rule plot_setup_ui:
#     input:
#         plotting_options="output/plotting_options.json"
#     output:
#         plot_settings="output/plot_setup.json"

# rule make_plot:
#     input:
#         scdata="output/scdata.Rds",
#         plot_setup="output/plot_setup.json",
#         plotting_options="output/plotting_options.json"
#     output:
#         plot_out="output/plot.out",
#         thumbnail="output/thumbnail.png",
#         plot_Rds="output/plot.Rds"
#     singularity:
#         "/app/archimedes-r.sif"
#     script:
#         "scripts/make_dittoSeq_plot.R"

