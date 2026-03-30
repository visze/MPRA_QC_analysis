rule activity_main:
    input:
        activity=lookup(
            within=activity_files,
            query="file == 'activity_df'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_main_plots(activity_files),
            file_type=["pdf", "eps", "png", "svg"],
        ),
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
         python {input.script} activity main --activity {input.activity} --output-path {params.outdir}
        """


rule activity_control_boxplots:
    input:
        activity=lookup(
            within=activity_files,
            query="file == 'activity_df'",
            cols=["path"],
        ),
        controls=lookup(
            within=activity_files,
            query="file == 'control_df'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_control_boxplots_plots(activity_files),
            file_type=["pdf", "eps", "png", "svg"],
        ),
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
         python {input.script} activity control-boxplots --activity {input.activity} --controls {input.controls} --output-path {params.outdir}
        """


rule activity_replicability_by_activity:
    input:
        activity=lookup(
            within=activity_files,
            query="file == 'activity_df'",
            cols=["path"],
        ),
        activity_per_rep=lookup(
            within=activity_files,
            query="file == 'activity_per_rep'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_replicability_by_activity_plots(activity_files),
            file_type=["pdf", "eps", "png", "svg"],
        ),
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
         python {input.script} activity replicability-by-activity --activity {input.activity} --activity-per-rep {input.activity_per_rep} --output-path {params.outdir}
        """


rule activity_gc_content_bias:
    input:
        design=lookup(
            within=activity_files,
            query="file == 'cCRE_fasta'",
            cols=["path"],
        ),
        activity=lookup(
            within=activity_files,
            query="file == 'activity_df'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_gc_content_bias_plots(activity_files),
            file_type=["pdf", "eps", "png", "svg"],
        ),
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
         python {input.script} activity gc-content-bias --activity {input.activity} --design {input.design} --output-path {params.outdir}
        """


rule activity_ratio_correlation_between_replicates:
    input:
        activity_per_rep=lookup(
            within=activity_files,
            query="file == 'activity_per_rep'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_ratio_correlation_between_replicates_plots(
                activity_files
            ),
            file_type=["pdf", "eps", "png", "svg"],
        ),
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
         python {input.script} activity ratio-correlation-between-replicates --activity-per-rep {input.activity_per_rep} --output-path {params.outdir}
        """


rule activity_ratio_correlation_with_controls:
    input:
        controls=lookup(
            within=activity_files,
            query="file == 'control_df'",
            cols=["path"],
        ),
        activity_per_rep=lookup(
            within=activity_files,
            query="file == 'activity_per_rep'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_ratio_correlation_with_controls_plots(activity_files),
            file_type=["pdf", "eps", "png", "svg"],
        ),
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
         python {input.script} activity ratio-correlation-with-controls --activity-per-rep {input.activity_per_rep} --controls {input.controls} --output-path {params.outdir}
        """


rule activity_downsampling:
    input:
        downsampling_activity_path=lookup(
            within=activity_files,
            query="file == 'downsampling_activity_path'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_downsampling_plots(activity_files),
            file_type=["pdf", "eps", "png", "svg"],
        ),
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
         python {input.script} activity downsampling --downsampling-activity-path {input.downsampling_activity_path} --output-path {params.outdir}
        """


rule activity_reproducibility_by_sequencing_depth:
    input:
        downsampling_ratio_path=lookup(
            within=activity_files,
            query="file == 'downsampling_ratio_path'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_reproducibility_by_sequencing_depth_plots(
                activity_files
            ),
            file_type=["pdf", "eps", "png", "svg"],
        ),
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
         python {input.script} activity reproducibility-by-sequencing-depth --downsampling-ratio-path {input.downsampling_ratio_path} --output-path {params.outdir}
        """


rule activity_mimimise_noise:
    input:
        sdt_thresholds=lookup(
            within=activity_files,
            query="file == 'different_std_threshold_analysis'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_mimimise_noise_plots(activity_files),
            file_type=["pdf", "eps", "png", "svg"],
        ),
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
         python {input.script} activity mimimise-noise --sdt-thresholds {input.sdt_thresholds} --output-path {params.outdir}
        """


rule activity_screen_annotations:
    input:
        screen=lookup(
            within=activity_files,
            query="file == 'screen_df'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        association_script=getScript("association_analysis.py"),
        activity_script=getScript("activity_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const_lib=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/activity/{plot}.{file_type}",
            plot=get_activity_screen_annotations_plots(activity_files),
            file_type=["pdf", "eps", "png", "svg"],
        ),
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/activity/"),
    shell:
        """
         python {input.script} activity screen-annotations --screen {input.screen} --output-path {params.outdir}
        """
