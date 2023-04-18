import pathlib

PACKAGE = "FAM_Pipeline"

MODULE_DIR = pathlib.Path(__file__).parent

from . import famp_simulation
from . import famp_modeling
from . import famp_data_analysis
