"""
modules/reporters/openmm_interface.py

Reporters that interface with OpenMM - shown as an example reporter -
so-called YOUR trajectory reporter.
"""

import typing

from attrs import define, field

import openmm.app as openmm_app
import seekr.modules.reporters.structures as reporters_structures
import seekr.modules.reporters.openmm_interface as reporters_openmm_interface

import seekrsamplingmethodsplugin.modules.reporters.structures as plugin_reporters_structures

#@define
#class YOUR_reporter_openmm_interface(reporters_openmm_interface.Base_reporter_openmm_interface):
#    """
#    A custom reporter for SEEKR.
#    """
#    _reporter: plugin_reporters_structures.YOUR_reporter = field()
#    def add_reporter_to_simulation(
#            self, 
#            **kwargs: typing.Any
#            ) -> None:
#        """
#        Add the reporter to the OpenMM simulation object.
#        """
#       ...
#        return
    
#    def get_reporter_extension(self) -> str:
#        """
#        Get the file extension for the reporter.
#        """
#        return self._reporter.get_reporter_extension()
        
def assign_reporter_interface(
        reporter: reporters_structures.Base_reporter, 
        ) -> tuple[bool, typing.Optional[reporters_openmm_interface.Base_reporter_openmm_interface]]:
    """
    Given a sampling_method and an md_engine, confirm that the OpenMM
    MD engine is being used, and return the corresponding interface
    """
    is_found = False
    openmm_interface = None
    #if reporter.type == "your_reporter":
    #    is_found = True
    #    openmm_interface = YOUR_reporter_openmm_interface(reporter)
    
    return is_found, openmm_interface
