################################
#### Global functions       ####
################################

SCRIPTS_DIR = "../../scripts"
ENVS_DIR = "../envs"


def getWorkflowFile(dir_name, name):
    return workflow.source_path("%s/%s" % (dir_name, name))


def getScript(name):
    return getWorkflowFile(SCRIPTS_DIR, name)


def getCondaEnv(name):
    return getWorkflowFile(ENVS_DIR, name)


##### load config and sample sheets #####
from snakemake.utils import validate
import pandas as pd

# Workaround: validate() is broken from Snakemake 9.5.1 to snakemake 9.14.7 in remote jobs
if version.parse(snakemake.__version__) >= version.parse("9.5.1") and version.parse(
    snakemake.__version__
) <= version.parse("9.14.7"):
    from snakemake_interface_executor_plugins.settings import ExecMode

    # Use the global 'workflow' variable directly as recommended by Snakemake
    if workflow.remote_exec:
        old_exec_mode = workflow.exec_mode
        workflow.workflow_settings.exec_mode = ExecMode.DEFAULT
        validate(config, schema="../schemas/config.schema.yml")
        workflow.workflow_settings.exec_mode = old_exec_mode
    else:
        validate(config, schema="../schemas/config.schema.yml")
else:
    validate(config, schema="../schemas/config.schema.yml")

# load sample sheets
association_files = pd.DataFrame()
activity_files = pd.DataFrame()
if "association" in config:
    association_files = pd.read_csv(config["association"], delimiter="\t")
    validate(association_files, schema="../schemas/association_file.schema.yml")
if "activity" in config:
    activity_files = pd.read_csv(config["activity"], delimiter="\t")
    validate(activity_files, schema="../schemas/activity_file.schema.yml")


################################
#### HELPERS AND EXCEPTIONS ####
################################


##### Exceptions #####


##### Helpers #####


def get_association_final_plots(association_df: pd.DataFrame) -> list:

    if association_df.empty:
        return []
    else:
        if (
            "cCRE_fasta" in association_df["file"].values
            and "final_associations" in association_df["file"].values
        ):
            return [
                "BCs_per_cCRE",
                "Retained_cCREs",
                "PCR_bias_GC",
                "PCR_bias_G_stretches",
            ]
    return []


def get_association_before_minimum_observations_plots(
    association_df: pd.DataFrame,
) -> list:

    if association_df.empty:
        return []
    else:
        if "associations_before_minimum_observations" in association_df["file"].values:
            return ["Reads_per_association"]
    return []


def get_association_before_promiscuity_plots(association_df: pd.DataFrame) -> list:

    if association_df.empty:
        return []
    else:
        if "associations_before_promiscuity" in association_df["file"].values:
            return ["cCREs_per_BC"]
    return []


def get_association_downsampling_plots(association_df: pd.DataFrame) -> list:

    if association_df.empty:
        return []
    else:
        if "associations_downsampling_path" in association_df["file"].values:
            return ["Downsampling_Retained_cCREs", "Downsampling_Barcodes_per_cCRE"]
    return []
