"""
modules/collective_variables/openmm_interface.py

The OpenMM interface for the collective variables defined in this
plugin.
"""

import typing

from attrs import define

import seekr.modules.collective_variables.base as cv_base 
import seekr.modules.collective_variables.openmm_interface as seekr_cv_openmm_interface
import seekr.modules.engines.structures as engines_structures

# One or more such files such modules as imported below must exist and 
#  define a collective variable object
# import seekrsamplingmethodsplugin.modules.collective_variables.YOUR_cv as YOUR_cv
# import seekrsamplingmethodsplugin.modules.structures as plugin_cv_structures

# @define
# class Your_cv_openmm_interface:
#     """
#     Put your CV's OpenMM interface here.
#     """
#     ...

# Create a typing union between the OpenMM interface objects in seekr, and those
# defined here. For now, just include those in seekr.
#CV_openmm_interface = typing.Union[seekr_cv_openmm_interface.CV_openmm_interface,
#                                   Your_cv_openmm_interface]
CV_openmm_interface = typing.Union[seekr_cv_openmm_interface.CV_openmm_interface]

def assign_cv_interface(
        collective_variable: cv_base.Base_CV,
        md_engine: engines_structures.Base_engine,
        root_directory: str) -> CV_openmm_interface | None:
    """
    Given this plugin's collective variable object, assign the appropriate 
    OpenMM interface.
    """
    is_found = False
    openmm_interface = None
    # Detect if it's a collective variable defined in this plugin
    # if collective_variable.type == "SOMETHING":
    #     ...
    #     is_found = True
    #     openmm_interface = Your_cv_openmm_interface(...)

    return is_found, openmm_interface

