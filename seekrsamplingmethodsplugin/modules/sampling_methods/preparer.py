"""
modules/sampling_methods/preparer.py

Prepare the plugin sampling method(s)
"""

import typing
from attrs import define, field, validators
import seekr.modules.sampling_methods.preparer as sampling_methods_preparer

import seekrsamplingmethodsplugin.modules.sampling_methods.structures \
    as plugin_sampling_methods_structures

@define
class Random_accelerated_md_sampling_method_input(
    sampling_methods_preparer.Base_md_sampling_method_input):
    """
    A Random accelerated molecular dynamics (RAMD) sampling method input.
    """
    type: typing.Literal["random_accelerated"] = "random_accelerated"
    scale_type: str = "molecular_dynamics"
    ligand_selection_name: str = field(default="ligand", validator=validators.instance_of(str))
    receptor_selection_name: str | None = field(default="receptor", validator=validators.instance_of(str))
    steps_per_update: int = field(default=50,
                                  validator=[validators.instance_of(int),
                                             validators.gt(0)])
    force_magnitude: float = field(default=14.0,
                                   validator=[validators.instance_of(float),
                                              validators.gt(0)])
    cutoff_distance: float = field(default=0.0025,
                                   validator=[validators.instance_of(float),
                                              validators.gt(0)])
    # TODO: going to use a completion criteria for max distance
    ramp_start_step: int = field(default=0,
                                validator=[validators.instance_of(int),
                                           validators.ge(0)])
    ramp_force_increase_per_step: float = field(default=4.0e-7,
                                                validator=validators.instance_of(float))

    def make_model_sampling_method(self):
        """
        Make the model sampling method object for this input.
        """
        sampling_method = plugin_sampling_methods_structures.Random_accelerated_md_sampling_method()
        sampling_method.ligand_selection_name = self.ligand_selection_name
        sampling_method.receptor_selection_name = self.receptor_selection_name
        sampling_method.steps_per_update = self.steps_per_update
        sampling_method.force_magnitude = self.force_magnitude
        sampling_method.cutoff_distance = self.cutoff_distance
        #sampling_method.max_distance = self.max_distance
        sampling_method.ramp_start_step = self.ramp_start_step
        sampling_method.ramp_force_increase_per_step = self.ramp_force_increase_per_step
        return sampling_method

sampling_methods: typing.Dict[str, typing.Dict[str, typing.Any]] \
    = {"molecular_dynamics": {"random_accelerated": Random_accelerated_md_sampling_method_input, }}
