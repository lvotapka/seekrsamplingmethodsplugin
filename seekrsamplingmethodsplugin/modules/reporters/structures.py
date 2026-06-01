"""
modules/reporters/structures.py

Reporter structure that could be defined for a plugin.
"""

import typing

from attrs import define, field, validators

import seekr.modules.reporters.structures as reporters_structures

#@define
#class YOUR_reporter:
#    """
#    A custom reporter for SEEKR.
#    """
#    type: typing.Literal["your_reporter"] = "your_reporter"
#    interval: int = field(default=1000,
#                          validator=[validators.instance_of(int),
#                                     validators.gt(0)])
#    ...
#    def get_reporter_extension(self) -> str:
#        """
#        Get the file extension for the reporter.
#        """
#        return "WHATEVER"

# This dictionary is needed by seekr to structure/unstructure data for the model file,
#  and its key is the .type attribute of the object, while the value is the class of
#  the object.
#reporters: typing.Dict[str, typing.Type[reporters_structures.Base_reporter]] \
#    = {"your_reporter": YOUR_reporter, }
# Empty for now
reporters: typing.Dict[str, typing.Type[reporters_structures.Base_reporter]] \
    = {}
