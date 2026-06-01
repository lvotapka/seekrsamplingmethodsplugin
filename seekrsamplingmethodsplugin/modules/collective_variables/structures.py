"""
modules/collective_variables/structures.py

Structures for collective variables within this plugin in the model.
"""

import typing

import seekr.modules.collective_variables.base as seekr_cv_base

# This file should be created and filled out in this directory
#import seekrsamplingmethodsplugin.modules.collective_variables.YOUR_cv as YOUR_cv

#Collective_variable = typing.Union[seekr_cv_structures.Collective_variable,
#                                   YOUR_cv.YOUR_CV]
# This dictionary is needed by seekr to structure/unstructure data for the model file,
#  and its key is the .type attribute of the object, while the value is the class of
#  the object.
#cvs: typing.Dict[str, typing.Type[seekr_cv_base.Base_CV]] \
#    = {"YOUR_cv": YOUR_cv.YOUR_CV, }
# Leave empty for now
cvs: typing.Dict[str, typing.Type[seekr_cv_base.Base_CV]] \
    = {}
