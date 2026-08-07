"""
modules/reporters/preparer.py

Prepare the reporter inputs for simulations for the plugin.
"""

import typing

from attrs import define, field, validators

import seekr.modules.reporters.preparer as reporters_preparer
import seekr.modules.engines.preparer as engines_preparer

import seekrsamplingmethodsplugin.modules.reporters.structures \
    as plugin_reporters_structures
import seekrsamplingmethodsplugin.modules.reporters.features as plugin_features
import seekrsamplingmethodsplugin.modules.reporters.conditions as plugin_conditions


@define
class Feature_reporter_input(reporters_preparer.Base_reporter_input):
    """
    Input for the feature reporter.
    """
    type: typing.Literal["feature"] = "feature"
    name: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)))
    interval: int = field(
        default=1000,
        validator=[validators.instance_of(int), validators.gt(0)])
    # Annotated as Any so cattrs leaves dicts for structure_feature_list.
    # min_len is enforced in make_reporter so callers can assign features after init.
    features: list = field(
        factory=list,
        converter=plugin_features.structure_feature_list)

    def make_reporter(
            self,
            unpartitioned_settings: engines_preparer.MD_unpartitioned_settings_input
            ) -> plugin_reporters_structures.Feature_reporter:
        """
        Make the reporter object for the model file.
        """
        assert len(self.features) > 0, \
            "Feature_reporter_input must define at least one feature."
        return plugin_reporters_structures.Feature_reporter(
            name=self.name,
            interval=self.interval,
            features=list(self.features),
        )


@define
class Conditional_position_reporter_input(reporters_preparer.Base_reporter_input):
    """
    Input for the conditional position reporter.
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
    # Selection name string (resolved in make_reporter), or None.
    atom_subset: str | None = field(
        default=None,
        validator=validators.optional(validators.instance_of(str)))
    operator: typing.Literal["and", "or"] = field(
        default="and",
        validator=validators.in_(["and", "or"]))
    conditions: list = field(
        factory=list,
        converter=plugin_conditions.structure_condition_list)

    def make_reporter(
            self,
            unpartitioned_settings: engines_preparer.MD_unpartitioned_settings_input
            ) -> plugin_reporters_structures.Conditional_position_reporter:
        """
        Make the reporter object for the model file.
        """
        assert len(self.conditions) > 0, \
            "Conditional_position_reporter_input must define at least one condition."
        if self.atom_subset is not None:
            atom_indices = unpartitioned_settings.get_atom_indices_by_selection_name(
                self.atom_subset)
        else:
            atom_indices = None
        return plugin_reporters_structures.Conditional_position_reporter(
            name=self.name,
            interval=self.interval,
            extension=self.extension,
            enforce_periodic_box=self.enforce_periodic_box,
            atom_subset=atom_indices,
            operator=self.operator,
            conditions=list(self.conditions),
        )


# This dictionary is needed by seekr to structure/unstructure data for the input file,
#  and its key is the .type attribute of the object, while the value is the class of
#  the object.
reporter_inputs: \
    typing.Dict[str, typing.Type[reporters_preparer.Base_reporter_input]] \
        = {
            "feature": Feature_reporter_input,
            "conditional_position": Conditional_position_reporter_input,
        }
