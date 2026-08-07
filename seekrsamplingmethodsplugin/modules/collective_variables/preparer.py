"""
modules/collective_variables/preparer.py

Collective variable objects for the prepare stage of seekr.
"""

import typing

import seekr.modules.collective_variables.base as seekr_cv_base

# The file imported below must be created in this plugin
import seekrsamplingmethodsplugin.modules.collective_variables\
    .rmsd_no_align_cv as rmsd_no_align_cv
import seekrsamplingmethodsplugin.modules.collective_variables\
    .min_pair_dist_cv as min_pair_dist_cv

# This object connects the name of CV object to the CV_input object:
# The former is defined for the prepare stage, and is simpler in general.
# The latter is defined for the model, and tends to include more attributes
#  and methods.

cv_inputs: typing.Dict[str, typing.Type[seekr_cv_base.Base_CV_input]] \
     = {"rmsd_no_align": rmsd_no_align_cv.RMSD_no_align_CV_input,
        "min_pair_dist": min_pair_dist_cv.Min_pair_dist_CV_input} 
