rule activity:
    """
    Activity workflow to perform the activity step of the MPRA QC analysis.
    """
    conda:
        getCondaEnv("default.yml")
    input:
        activity_file=config["activity"] if "activity" in config else [],
        script=getScript("activity_analysis.py"),
    output:
        touch("results/{project}/activity.done"),
    shell:
        """
        python {input.script} {input.activity_file}
        """
