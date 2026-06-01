"""
modules/reporters/preparer.py

Prepare the reporter inputs for simulations for the plugin
"""

import typing
from attrs import define, field, validators, Factory
import seekr.modules.reporters.preparer as reporters_preparer

import seekrsamplingmethodsplugin.modules.reporters.structures \
    as plugin_reporters_structures


#@define
#class YOUR_reporter_input(reporters_preparer.Base_reporter_input):
#    """
#    Input for YOUR reporter.
#    """
#    type: typing.Literal["your_reporter"] = "your_reporter"
#    interval: int = field(default=1000,
#                          validator=[validators.instance_of(int),
#                                     validators.gt(0)])
#   ...
#    def make_reporter(
#            self,
#            unpartitioned_settings: engines_preparer.MD_unpartitioned_settings_input
#        ) -> plugin_reporters_structures.YOUR_reporter:
#        """
#        Make the reporter object for the model file.
#        """
#        reporter = plugin_reporters_structures.YOUR_reporter()
#        reporter.interval = self.interval
#        return reporter
        
# This dictionary is needed by seekr to structure/unstructure data for the input file,
#  and its key is the .type attribute of the object, while the value is the class of
#  the object.
#reporter_inputs: typing.Dict[str, typing.Type[reporters_preparer.Base_reporter_input]] \
#    = {"xtc": XTC_reporter_input, }
# Empty for now
reporter_inputs: \
    typing.Dict[str, typing.Type[reporters_preparer.Base_reporter_input]] \
        = {}
