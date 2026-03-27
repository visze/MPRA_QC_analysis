rule association_final:
    input:
        design=lookup(
            within=association_files, query="file == 'cCRE_fasta'", cols=["path"]
        ),
        associations=lookup(
            within=association_files,
            query="file == 'final_associations'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const=getScript("const.py"),
    output:
        touch("results/{project}/association.done"),
        expand(
            "results/{{project}}/association/{plot}.{file_type}",
            plot=get_association_final_plots(association_files),
            file_type=["pdf", "eps", "png", "svg"],
        ),
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/association/"),
    shell:
        """
         python {input.script}  association final --design {input.design} --associations {input.associations} --output-path {params.outdir}
        """


rule association_before_minimum_observations:
    input:
        associations=lookup(
            within=association_files,
            query="file == 'associations_before_minimum_observations'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/association/{plot}.{file_type}",
            plot=get_association_before_minimum_observations_plots(association_files),
            file_type=["pdf", "eps", "png", "svg"],
        ),
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/association/"),
    shell:
        """
         python {input.script}  association before_minimum_observations --associations {input.associations} --output-path {params.outdir}
        """


rule association_before_promiscuity:
    input:
        associations=lookup(
            within=association_files,
            query="file == 'associations_before_promiscuity'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/association/{plot}.{file_type}",
            plot=get_association_before_promiscuity_plots(association_files),
            file_type=["pdf", "eps", "png", "svg"],
        ),
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/association/"),
    shell:
        """
         python {input.script}  association before_promiscuity --associations {input.associations} --output-path {params.outdir}
        """


rule association_downsampling:
    input:
        design=lookup(
            within=association_files, query="file == 'cCRE_fasta'", cols=["path"]
        ),
        downsampling_path=lookup(
            within=association_files,
            query="file == 'associations_downsampling_path'",
            cols=["path"],
        ),
        script=getScript("mpra_qc_analysis.py"),
        plot_lib=getScript("plot_lib.py"),
        const=getScript("const.py"),
    output:
        expand(
            "results/{{project}}/association/{plot}.{file_type}",
            plot=get_association_downsampling_plots(association_files),
            file_type=["pdf", "eps", "png", "svg"],
        ),
    conda:
        getCondaEnv("default.yml")
    params:
        outdir=directory("results/{project}/association/"),
    shell:
        """
         python {input.script}  association downsampling --design {input.design} --downsampling-path {input.downsampling_path} --output-path {params.outdir}
        """
