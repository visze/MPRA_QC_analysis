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
if "association" in config:
    experiment = pd.read_csv(config["association"], delimiter="\t")
    validate(experiment, schema="../schemas/association_file.schema.yml")
if "activity" in config:
    experiment = pd.read_csv(config["activity"], delimiter="\t")
    validate(experiment, schema="../schemas/activity_file.schema.yml")


################################
#### HELPERS AND EXCEPTIONS ####
################################


##### Exceptions #####
class MissingAssignmentInConfigException(Exception):
    """
    Exception raised for if no assignment file is set in the config.

    Args:
        Exception ([type]): Exception class cast

    Attributes:
        config_name (string): name of the configuration which assignment is missing.
    """

    def __init__(self, config_name):
        self.config_name = config_name

    def __str__(self):
        return "Config %s has no assignment file defined!" % (self.config_name)


class MissingVariantInConfigException(Exception):
    """
    Exception raised for if no variants config.

    Args:
        Exception ([type]): Exception class cast

    Attributes:
        config_name (string): name of the configuration which assignment is missing.
    """

    def __init__(self, config_name):
        self.config_name = config_name

    def __str__(self):
        return "Config %s has no variants defined!" % (self.config_name)
