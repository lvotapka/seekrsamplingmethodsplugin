"""
modules/restraints/structures.py

Restraint structure(s) that could be defined for a plugin.
"""

import typing
from attrs import define, field, validators, Factory

import seekr.modules.restraints.structures as restraints_structures

#@define
#class YOUR_restraint(restraints_structures.Base_restraint):
#    """
#    Restraint that keeps a subset of atoms close to some desired positions
#    """
#    type: typing.Literal["your_restraint"] = "your_restraint"
#    # Define any parameters needed for your restraint here.
#    ...

                                            
# This dictionary is needed by seekr to structure/unstructure data for the model file,
#  and its key is the .type attribute of the object, while the value is the class of
#  the object.
#restraints: typing.Dict[str, typing.Type[restraints_structures.Base_restraint]] \
#    = {"your_restraint": YOUR_restraint, }
restraints: typing.Dict[str, typing.Type[restraints_structures.Base_restraint]] \
    = {}