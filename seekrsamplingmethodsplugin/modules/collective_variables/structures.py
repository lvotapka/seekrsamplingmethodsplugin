"""
modules/collective_variables/structures.py

Structures for collective variables within this plugin in the model.
"""

import typing

import seekr.modules.collective_variables.base as seekr_cv_base

# This file should be created and filled out in this directory
import seekrsamplingmethodsplugin.modules.collective_variables\
    .rmsd_no_align_cv as rmsd_no_align_cv
import seekrsamplingmethodsplugin.modules.collective_variables\
    .min_pair_dist_cv as min_pair_dist_cv

#Collective_variable = typing.Union[seekr_cv_structures.Collective_variable,
#                                   YOUR_cv.YOUR_CV]
# This dictionary is needed by seekr to structure/unstructure data for the model file,
#  and its key is the .type attribute of the object, while the value is the class of
#  the object.

cvs: typing.Dict[str, typing.Type[seekr_cv_base.Base_CV]] \
    = {"rmsd_no_align_cv": rmsd_no_align_cv.RMSD_no_align_CV,
       "min_pair_dist_cv": min_pair_dist_cv.Min_pair_dist_CV}
