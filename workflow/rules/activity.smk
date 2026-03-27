rule activity:
    input:
        activity_file=config["activity"] if "activity" in config else [],
        script=getScript("activity_analysis.py"),
        const=getScript("const.py"),
    output:
        touch("results/{project}/activity.done"),
    conda:
        getCondaEnv("default.yml")
    shell:
        """
        python {input.script} run run-all --config {input.activity_file}
        """
