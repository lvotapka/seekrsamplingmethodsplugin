"""
modules/reporters/structures.py

Reporter structures defined for this plugin.
"""

import typing

from attrs import define, field, validators

import seekr.modules.reporters.structures as reporters_structures

import seekrsamplingmethodsplugin.modules.reporters.features as plugin_features
import seekrsamplingmethodsplugin.modules.reporters.conditions as plugin_conditions


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


@define
class Conditional_position_reporter(reporters_structures.Base_reporter):
    """
    Position reporter that writes frames only while CV conditions hold.

    Frames for a sojourn are stored in a temporary file and discarded if the
    condition becomes false. The temp file is promoted to the final trajectory
    only if the condition is still true when the stage ends.
    """
    type: typing.Literal["conditional_position"] = "conditional_position"
    name: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)))
    interval: int = field(
        default=1000,
        validator=[validators.instance_of(int), validators.gt(0)])
    extension: str = field(
        default="dcd",
        validator=validators.in_(["dcd", "pdb"]))
    # False by default so egress trajectories are continuous (no PBC wraps)
    # for pathway clustering / RMSD analysis.
    enforce_periodic_box: bool = field(
        default=False,
        validator=validators.instance_of(bool))
    atom_subset: typing.Optional[typing.List[int]] = field(
        default=None,
        validator=validators.optional(
            validators.deep_iterable(
                member_validator=validators.instance_of(int),
                iterable_validator=validators.instance_of(list))))
    operator: typing.Literal["and", "or"] = field(
        default="and",
        validator=validators.in_(["and", "or"]))
    # Assign-after-init friendly; length checked when the OpenMM reporter starts.
    conditions: list = field(
        factory=list,
        converter=plugin_conditions.structure_condition_list)

    def get_reporter_extension(self) -> str:
        """
        Get the file extension for the reporter.
        """
        return self.extension


# This dictionary is needed by seekr to structure/unstructure data for the model file,
#  and its key is the .type attribute of the object, while the value is the class of
#  the object.
reporters: typing.Dict[str, typing.Type[reporters_structures.Base_reporter]] \
    = {
        "feature": Feature_reporter,
        "conditional_position": Conditional_position_reporter,
    }
