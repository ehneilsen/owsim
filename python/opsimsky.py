"""A wrapper around the opsim4 calculation of sky brightness
"""

import sys
import os
import importlib.util
from types import ModuleType

import numpy as np

# The LSST code depends on a couple procudecs defined in the LSST
# stack. To avoid needing to load and install the entire stack, create
# adequate implementations and fiddle python's sys.modules to make it
# look like they came from the stack.

#
# create getPackageDir, and make it look like it came from lsst.utils
#

def getPackageDir(*args, **kwargs):
    raise Exception("I am fake")

lsst_module = ModuleType("lsst")
utils_module = ModuleType("utils")
lsst_module.utils = utils_module
utils_module.getPackageDir = getPackageDir
sys.modules[lsst_module.__name__] = lsst_module
sys.modules['lsst.utils'] = utils_module

#
# create _angularSeparation, and make it look like it came from lsst.sims.utils
# 

def _angularSeparation(long1, lat1, long2, lat2):
# Copied from sims_utils/python/lsst/sims/utils/CoordinateTransformations
    t1 = np.sin(lat2/2.0 - lat1/2.0)**2
    t2 = np.cos(lat1)*np.cos(lat2)*np.sin(long2/2.0 - long1/2.0)**2
    _sum = t1 + t2

    if isinstance(_sum, numbers.Number):
        if _sum<0.0:
            _sum = 0.0
    else:
        _sum = np.where(_sum<0.0, 0.0, _sum)

    return 2.0*np.arcsin(np.sqrt(_sum))

sims_module = ModuleType("sims")
lsst_module.sims = sims_module
simsutils_module = ModuleType("utils")
lsst_module.sims.utils = simsutils_module
lsst_module.sims.utils._angularSeparation = _angularSeparation
sys.modules['lsst.sims'] = sims_module
sys.modules['lsst.sims.utils'] = simsutils_module

opsim_skymodle_python = os.environ['OPSIM_SKYMODEL_PYTHON']
spec = importlib.util.spec_from_file_location("SkyModelPre", opsim_skymodle_python)
SkyModelPre = importlib.util.module_from_spec(spec)
spec.loader.exec_module(SkyModelPre)
