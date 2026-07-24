"""
modules/reporters/structures.py

Reporter structures defined for this plugin.
"""

import typing

from attrs import define, field, validators

import seekr.modules.reporters.structures as reporters_structures

import seekrsamplingmethodsplugin.modules.reporters.features as plugin_features


@define
class Feature_reporter(reporters_structures.Base_reporter):
    """
    Reporter that writes a CSV table of geometric features vs time.
    """
    type: typing.Literal["feature"] = "feature"
    name: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)))
    interval: int = field(
        default=1000,
        validator=[validators.instance_of(int), validators.gt(0)])
    # Annotated as Any so cattrs leaves dicts for structure_feature_list.
    features: list = field(
        factory=list,
        converter=plugin_features.structure_feature_list,
        validator=validators.min_len(1))

    def get_reporter_extension(self) -> str:
        """
        Get the file extension for the reporter.
        """
        return "csv"


# This dictionary is needed by seekr to structure/unstructure data for the model file,
#  and its key is the .type attribute of the object, while the value is the class of
#  the object.
reporters: typing.Dict[str, typing.Type[reporters_structures.Base_reporter]] \
    = {"feature": Feature_reporter}
