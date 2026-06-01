"""
modules/restraints/preparer.py

Prepare the restraints inputs for simulations for this example restraint
of the z-coordinate of a group of atoms.
"""

import typing
from attrs import define, field, validators
import seekr.modules.restraints.preparer as restraints_preparer
import seekr.modules.engines.preparer as engines_preparer

import seekrsamplingmethodsplugin.modules.restraints.structures \
    as plugin_restraints_structures

#@define
#class YOUR_restraint_input(restraints_preparer.Base_restraint_input):
#    """
#    Input for custom restraint in plugin.
#    """
#    type: typing.Literal["your_restraint"] = "your_restraint"
#   ...
#    def make_restraint(
#            self, 
#            unpartitioned_settings: engines_preparer.MD_unpartitioned_settings_input
#            ) -> plugin_restraints_structures.YOUR_restraint:
#        """
#        Make the restraint object for the model.
#        """
#        restraint = plugin_restraints_structures.YOUR_restraint(
#            # Pass any parameters needed for your restraint here.
#        )
#       ...
#        return restraint
        
# This dictionary is needed by seekr to structure/unstructure data for the input file,
#  and its key is the .type attribute of the object, while the value is the class of
#  the object.
#restraint_inputs: typing.Dict[str, typing.Type[restraints_preparer.Base_restraint_input]] \
#    = {"your_restraint": YOUR_restraint_input, }
restraint_inputs: typing.Dict[str, typing.Type[restraints_preparer.Base_restraint_input]] \
    = {}
