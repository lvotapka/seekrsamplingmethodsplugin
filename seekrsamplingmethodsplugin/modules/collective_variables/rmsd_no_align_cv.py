"""
modules/collective_variables/rmsd_no_align_cv.py

Model and input classes for a RMSD collective variable where no 
alignment is performed. This precludes the need for the RMSD_plus 
seekr plugin, which is hard to install.
"""

import typing

from attrs import define, field, validators
import mdtraj

import seekr.modules.base as seekr_base
import seekr.modules.collective_variables.base as seekr_cv_base

if typing.TYPE_CHECKING:
    import seekr.modules.structures as structures
    
@define
class RMSD_no_align_CV(seekr_cv_base.Base_CV):
    """
    Root mean squared deviation (RMSD) collective variable where 
    no alignment is performed.
    """
    # Required fields must come before any field with a default (attrs rule).
    index: int = field(validator=[validators.instance_of(int),
                                  validators.ge(0)])
    group_name: str = field(validator=validators.instance_of(str))
    directory: str = field(default="",
                           validator=validators.instance_of(str))
    ref_structure: str = field(default="",
                               validator=validators.instance_of(str))
    name: str = field(default="",
                        validator=validators.instance_of(str))
    unit: str = field(default="nm",
                      validator=validators.instance_of(str))
    num_groups: typing.Literal[1] = 1
    variable_name = "v"
    periodic: bool = False
    min_value: float = 0.0
    max_value: float = 4.0
    type: typing.Literal["rmsd_no_align"] = "rmsd_no_align"
    
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
        group_atom_indices = md_engine.get_atom_indices_by_selection_name(
            self.group_name)
        ref_traj_filename = os.path.join(
            model.directory, "collective_variable_sets", "main", 
            self.directory, self.ref_structure)
        ref_traj = mdtraj.load(ref_traj_filename)
        ref_xyz = ref_traj.xyz[0]
        value = 0.0
        for index in group_atom_indices:
            increment = (traj_xyz[index,0] - ref_xyz[index,0])**2 \
                + (traj_xyz[index,1] - ref_xyz[index,1])**2 \
                + (traj_xyz[index,2] - ref_xyz[index,2])**2
            value += increment
        
        rmsd = np.sqrt(value/len(group_atom_indices))
        return rmsd
    
    def get_atom_group_names(self) -> typing.List[str]:
        """
        Get the atom group names for this CV.
        """
        return [self.group_name]
        
@define
class RMSD_no_align_CV_input(seekr_cv_base.Base_CV_input):
    """
    RMSD_no_align_CV input which is used as an object for the 
    input JSON or API within the prepare stage of a seekr calculation.
    """
    name: str = field(default="",
                        validator=validators.instance_of(str))
    group_name: str = field(default="",
                             validator=validators.instance_of(str))
    ref_structure: str = field(default="",
                             validator=validators.instance_of(str))
    type: typing.Literal["rmsd_no_align"] = "rmsd_no_align"
    _index: int = field(default=0,
                          validator=validators.instance_of(int))
    min_value: float = field(default=0.0,
                             validator=validators.instance_of(float))
    max_value: float = field(default=4.0,
                                validator=validators.instance_of(float))
    
    def make_cv_object(
            self, 
            index: int,
            root_directory: str) -> RMSD_no_align_CV:
        """
        Create a Z_CV object from the Input to be placed into the Model.
        """
        self._index = index
        directory = f"collective_variable_{index}"
        if self.name == "":
            name = f"rmsd_no_align_{self.group_name}"
        else:
            name = self.name
        cv_sets_main_directory = os.path.join(
            root_directory, "collective_variable_sets", "main")
        cv_tree = {directory: {}}
        seekr_base.make_tree(cv_sets_main_directory, cv_tree)
        rmsd_ref_pdb = f"rmsd_reference_cv_{index}.pdb"
        absolute_RMSD_ref_pdb = os.path.join(
            cv_sets_main_directory, directory, rmsd_ref_pdb)
        cv = RMSD_no_align_CV(
            index, self.group_name,
            directory=directory, ref_structure=rmsd_ref_pdb)
        cv.name = name
        cv.min_value = self.min_value
        cv.max_value = self.max_value
        assert os.path.exists(self.ref_structure), \
            f"RMSD_no_align_CV Reference structure file {self.ref_structure} "\
            f"does not exist."
        shutil.copyfile(self.ref_structure, absolute_RMSD_ref_pdb)
        return cv