"""
modules/collective_variables/min_pair_dist_cv.py

Model and input classes for a minimum pair distance collective 
variable where two groups are provided, and the minimum distance 
between any two atoms across the two groups is calculated.
"""

import typing

from attrs import define, field, validators
import mdtraj

import seekr.modules.base as seekr_base
import seekr.modules.collective_variables.base as seekr_cv_base

if typing.TYPE_CHECKING:
    import seekr.modules.structures as structures

@define
class Min_pair_dist_CV(seekr_cv_base.Base_CV):
    """
    Minimum pair distance collective variable where two groups are 
    provided, and the minimum distance between any two atoms across 
    the two groups is calculated.
    """
    index: int = field(validator=[validators.instance_of(int),
                                  validators.ge(0)])
    group1_name: str = field(default="",
                            validator=validators.instance_of(str))
    group2_name: str = field(default="",
                            validator=validators.instance_of(str))
    exponent: int = field(default=50)
    cutoff_distance: float = field(default=1.0)
    name: str = field(default="",
                        validator=validators.instance_of(str))
    unit: str = field(default="nm",
                      validator=validators.instance_of(str))
    num_groups: typing.Literal[2] = 2
    variable_name = "r"
    periodic: bool = False
    min_value: float = field(default=0.0,
                             validator=validators.instance_of(float))
    max_value: float = field(default=1.0,
                             validator=validators.instance_of(float))

    type: typing.Literal["min_pair_dist"] = "min_pair_dist"
    
    def is_periodic(self) -> bool:
        """
        Return whether this CV is periodic.
        """
        return self.periodic

    def get_min_value(self) -> float:
        """
        Return the minimum value of this CV.
        """
        return self.min_value
    
    def get_max_value(self) -> float:
        """
        Return the maximum value of this CV.
        """
        return self.max_value
    
    def get_mdtraj_cv_value(
            self, 
            model: "structures.Seekr_model", 
            traj: mdtraj.Trajectory, 
            frame_index: int = 0
            ) -> float:
        """
        Determine the current CV value for an mdtraj object.
        """
        md_engine = model.get_engine_by_scale_type("molecular_dynamics")
        group1_atom_indices = md_engine.get_atom_indices_by_selection_name(
            self.group1_name)
        group2_atom_indices = md_engine.get_atom_indices_by_selection_name(
            self.group2_name)
        positions = traj.xyz[frame_index]
        pos_subset_group1 = []
        pos_subset_group2 = []
        for atom_index in group1_atom_indices:
            pos_subset_group1.append(positions[atom_index].value_in_unit(unit.nanometers))
            
        for atom_index in group2_atom_indices:
            pos_subset_group2.append(positions[atom_index].value_in_unit(unit.nanometers))

        min_dist = 9e9
        for i in range(len(group1_atom_indices)):
            for j in range(len(group2_atom_indices)):
                dist = np.linalg.norm(pos_subset_group2[j] - pos_subset_group1[i])
                if dist < min_dist:
                    min_dist = dist
        return min_dist
    
    def get_atom_group_names(self) -> typing.List[str]:
        """
        Get the atom group names for this CV.
        """
        return [self.group1_name, self.group2_name]
        
@define
class Min_pair_dist_CV_input(seekr_cv_base.Base_CV_input):
    """
    Min_pair_dist_CV input which is used as an object for the 
    input JSON or API within the prepare stage of a seekr calculation.
    """
    name: str = field(default="",
                        validator=validators.instance_of(str))
    group1_name: str = field(default="",
                             validator=validators.instance_of(str))
    group2_name: str = field(default="",
                             validator=validators.instance_of(str))
    type: typing.Literal["min_pair_dist"] = "min_pair_dist"
    _index: int = field(default=0,
                          validator=validators.instance_of(int))
    min_value: float = field(default=0.0,
                             validator=validators.instance_of(float))
    max_value: float = field(default=3.0,
                                validator=validators.instance_of(float))
    
    def make_cv_object(
            self, 
            index: int,
            root_directory: str) -> Min_pair_dist_CV:
        """
        Create a Min_pair_dist_CV object from the Input to be placed into the Model.
        """
        self._index = index
        if self.name == "":
            name = f"min_pair_dist_{self.group1_name}_{self.group2_name}"
        else:
            name = self.name
        cv = Min_pair_dist_CV(index, self.group1_name, self.group2_name)
        cv.name = name
        cv.min_value = self.min_value
        cv.max_value = self.max_value
        return cv