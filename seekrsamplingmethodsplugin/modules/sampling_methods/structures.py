"""
modules/sampling_methods/structures.py

Sampling method that could be defined for a plugin.
"""

import typing
from attrs import define, field, validators

import seekr.modules.sampling_methods.structures \
    as sampling_methods_structures

@define
class Random_accelerated_md_sampling_method(sampling_methods_structures.Base_md_sampling_method):
    """
    Implements the random accelerated molecular dynamics (RAMD) sampling method. 
    In this method, a random force of a constant magnitude is applied to the ligand 
    in a random direction, and the direction of the force is changed when the ligand 
    moves a certain distance. This method is useful for exploring the conformational 
    space of the ligand and finding possible binding/unbinding pathways.
    """
    type: typing.Literal["random_accelerated"] = "random_accelerated"
    scale_type = "molecular_dynamics"
    ligand_selection_name: str = field(default="ligand", validator=validators.instance_of(str))
    receptor_selection_name: str = field(default="receptor", validator=validators.instance_of(str))
    steps_per_update: int = field(default=50,
                                  validator=[validators.instance_of(int),
                                             validators.gt(0)])
    force_magnitude: float = field(default=14.0,
                                   validator=[validators.instance_of(float),
                                              validators.gt(0)])
    # Minimum COM travel between force retargets, in nanometers (openmm_ramd
    # RAMD_cutoff_distance convention). Default 0.0025 nm = 0.025 Å.
    cutoff_distance: float = field(default=0.0025,
                                   validator=[validators.instance_of(float),
                                              validators.gt(0)])
    # TODO: going to use a completion criteria
    #max_distance: float = field(default=3.0,
    #                            validator=[validators.instance_of(float),
    #                                       validators.gt(0)])
    ramp_start_step: int = field(default=0,
                                validator=[validators.instance_of(int),
                                           validators.ge(0)])
    ramp_force_increase_per_step: float = field(default=4.0e-7,
                                                validator=validators.instance_of(float))
    def input_from_unpartitioned_settings_output_to_anchors(self) -> bool:
        return True
    
    def input_from_anchors_output_to_unpartitioned_settings(self) -> bool:
        return False
    
    def needs_to_be_in_anchor(self) -> bool:
        return False
    
    def get_output_positions_basename_noext(self) -> str:
        return "random_accelerated"
    
    def get_output_positions_format(self) -> str:
        return "random_accelerated.{}.swarm_{}"
    
    def get_output_positions_glob_list(
            self,
            swarm_id: int = 0
            ) -> list[str]:
        swarm_string = f"swarm_{swarm_id}"
        return [f"random_accelerated*{swarm_string}*.pdb", f"random_accelerated*{swarm_string}*.txt"]
    
    def get_final_pdb_output_basename(
            self, 
            toy_mode: bool,
            frame_number: int | None = None,
            ) -> str:
        if toy_mode:
            extension = "txt"
        else:
            extension = "pdb"
        if frame_number is None:
            return f"random_accelerated.{extension}"
        else:
            return f"random_accelerated.swarm_{frame_number}.{extension}"
    
    def get_output_md_parameters_topology(self) -> None:
        return None
    
                                            
# This dictionary is needed by seekr to structure/unstructure data for the model file,
#  and its key is the .type attribute of the object, while the value is the class of
#  the object.
sampling_methods: typing.Dict[str, typing.Dict[str, typing.Any]] \
    = {"molecular_dynamics": {"random_accelerated": Random_accelerated_md_sampling_method, }}
