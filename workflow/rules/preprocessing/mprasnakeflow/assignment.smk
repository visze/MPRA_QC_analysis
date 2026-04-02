rule preprocessing_mprasnakeflow_assignment_associations_before_promiscuity:
    input:
        barcodes_incl_other=config.get("mprasnakeflow", {})
        .get("assignment", {})
        .get("barcodes_incl_other", []),
    output:
        "results/{project}/preprocessing/mprasnakeflow/association/associations_before_promiscuity.csv.gz",
    log:
        "logs/preprocessing/mprasnakeflow/assignment/associations_before_promiscuity.{project}.log",
    conda:
        getCondaEnv("default.yml")
    shell:
        """
        (
            echo "barcode,cCRE,match_count"; 
            zcat {input.barcodes_incl_other} | egrep -v "\\sother\\sNA" | \
            cut -f 1,2 | uniq -c | awk -v "OFS=," '{{print $2,$3,$1}}';
        ) | gzip -c > {output} 2> {log}
        """


rule preprocessing_mprasnakeflow_assignment_associations_before_minimum_observations:
    input:
        association_before_promiscuity="results/{project}/preprocessing/mprasnakeflow/association/associations_before_promiscuity.csv.gz",
        script=getScript(
            "preprocessing/mprasnakeflow/create_before_minimum_associations.py"
        ),
    output:
        "results/{project}/preprocessing/mprasnakeflow/association/associations_before_minimum_observations.csv.gz",
    log:
        "logs/preprocessing/mprasnakeflow/assignment/associations_before_minimum_observations.{project}.log",
    conda:
        getCondaEnv("default.yml")
    params:
        threshold=config.get("mprasnakeflow", {})
        .get("assignment", {})
        .get("fraction", 0.75),
    shell:
        """
        python {input.script} --input {input.association_before_promiscuity} --output {output} --threshold {params.threshold} > {log} 2>&1
        """


rule preprocessing_mprasnakeflow_assignment_final_associations:
    input:
        assignment=config.get("mprasnakeflow", {})
        .get("assignment", {})
        .get("assignment_barcodes_with_ambiguous", []),
    output:
        "results/{project}/preprocessing/mprasnakeflow/association/final_associations.csv.gz",
    log:
        "logs/preprocessing/mprasnakeflow/assignment/final_associations.{project}.log",
    conda:
        getCondaEnv("default.yml")
    shell:
        """
        (
            echo "barcode,cCRE,match_count";
            zcat {input.assignment} | egrep -v "\\sother\\s" |  egrep -v "\\sambiguous\\s" | cut -f 1,2,4 | tr '\\t' ',' | sed 's|/[^/]*$||';
        ) | gzip -c > {output} 2> {log}
        """


# TODO downsampling associations for QC plots
