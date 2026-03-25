rule association:
    """
    Association workflow to perform the association step of the MPRA QC analysis.
    """
    conda:
        getCondaEnv("default.yml")
    input:
        association_file=config["association"] if "association" in config else [],
        script=getScript("association_analysis.py"),
    output:
        touch("results/{project}/association.done"),
    shell:
        """
        python {input.script} {input.association_file}
        """
