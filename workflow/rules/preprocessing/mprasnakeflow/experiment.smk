rule preprocessing_mprasnakeflow_experiment_bcalm_elements:
    input:
        reporter_experiment_barcode=config.get("mprasnakeflow", {})
        .get("experiment", {})
        .get("reporter_experiment_barcode", []),
        labels=config.get("mprasnakeflow", {}).get("experiment", {}).get("labels", []),
        script=getScript("preprocessing/mprasnakeflow/bcalm_elements.R"),
    output:
        bcalm_result="results/{project}/preprocessing/mprasnakeflow/activity/bcalm/elements.tsv.gz",
        output_volcano_plot="results/{project}/preprocessing/mprasnakeflow/activity/bcalm/volcano_plot.png",
        output_density_plot="results/{project}/preprocessing/mprasnakeflow/activity/bcalm/density_plot.png",
    log:
        "logs/preprocessing/mprasnakeflow/experiment/bcalm_elements.{project}.log",
    conda:
        getCondaEnv("bcalm.yml")
    params:
        test_label=config.get("mprasnakeflow", {})
        .get("experiment", {})
        .get("test_label", ""),
        control_label=config.get("mprasnakeflow", {})
        .get("experiment", {})
        .get("control_label", ""),
        percentile=config.get("mprasnakeflow", {})
        .get("experiment", {})
        .get("percentile", 0.95),
        normalize=(
            "TRUE"
            if config.get("mprasnakeflow", {})
            .get("experiment", {})
            .get("normalize", True)
            else "FALSE"
        ),
    shell:
        """
        Rscript {input.script} --counts {input.reporter_experiment_barcode} \
        --labels {input.labels} --test-label {params.test_label} --control-label {params.control_label} \
        --percentile {params.percentile} --normalize {params.normalize} \
        --output {output.bcalm_result} --output-volcano-plot {output.output_volcano_plot} --output-density-plot {output.output_density_plot} > {log} 2>&1
        """


rule preprocessing_mprasnakeflow_experiment_activity_df:
    input:
        reporter_experiment_barcode=config.get("mprasnakeflow", {})
        .get("experiment", {})
        .get("reporter_experiment_barcode", []),
        bcalm_result="results/{project}/preprocessing/mprasnakeflow/activity/bcalm/elements.tsv.gz",
        script=getScript("preprocessing/mprasnakeflow/create_activity_df.py"),
    output:
        activity_df="results/{project}/preprocessing/mprasnakeflow/activity/activity_df.csv.gz",
    log:
        "logs/preprocessing/mprasnakeflow/experiment/activity_df.{project}.log",
    conda:
        getCondaEnv("mpralib.yml")
    shell:
        """
        python {input.script} --counts {input.reporter_experiment_barcode} \
        --activity {input.bcalm_result} \
        --output {output.activity_df} > {log} 2>&1
        """
