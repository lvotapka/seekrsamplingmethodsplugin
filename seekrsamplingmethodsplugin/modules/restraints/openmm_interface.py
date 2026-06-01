"""
modules/restraints/openmm_interface.py

Restraints that interface with OpenMM.
"""

import typing

from attrs import define, field

import openmm
import openmm.unit as unit
import seekr.modules.restraints.structures as restraints_structures
import seekr.modules.restraints.openmm_interface as restraints_openmm_interface

import seekrsamplingmethodsplugin.modules.restraints.structures as plugin_restraints_structures

#@define
#class YOUR_restraint_openmm_interface(restraints_openmm_interface.Base_restraint_openmm_interface):
#    """
#    OpenMM interface for YOUR_restraint.
#    """
#    _restraint: plugin_restraints_structures.YOUR_restraint = field()
#    def add_forces_to_system(
#            self, 
#            **kwargs: typing.Any
#            ) -> None:
#        """
#        Add force(s) to an OpenMM System object.
#        """
#        ...
#        return None

def assign_restraints_interface(
        restraint: restraints_structures.Base_restraint
        ) -> tuple[bool, typing.Optional[restraints_openmm_interface.Base_restraint_openmm_interface]]:
    """
    Given a sampling_method and an md_engine, confirm that the OpenMM
    MD engine is being used, and return the corresponding interface
    """
    is_found: bool = False
    openmm_interface = None
    #if restraint.type == "your_restraint":
    #    is_found = True
    #    openmm_interface = YOUR_restraint_openmm_interface(restraint)
        
    return is_found, openmm_interface
