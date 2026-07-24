"""
modules/collective_variables/openmm_interface.py

The OpenMM interface for the collective variables defined in this
plugin.
"""

import typing

from attrs import define

import seekr.modules.collective_variables.base as cv_base 
import seekr.modules.collective_variables.openmm_interface \
    as seekr_cv_openmm_interface
import seekr.modules.engines.structures as engines_structures

import seekrsamplingmethodsplugin.modules.collective_variables\
    .rmsd_no_align_cv as rmsd_no_align_cv
import seekrsamplingmethodsplugin.modules.collective_variables\
    .min_pair_dist_cv as min_pair_dist_cv
import seekrsamplingmethodsplugin.modules.structures as plugin_cv_structures

@define
class RMSD_no_align_CV_openmm_interface:
    collective_variable: rmsd_no_align_cv.RMSD_no_align_CV = field()
    group_atom_indices: typing.List[int] = field()
    ref_structure: str = field()
    CV_EXPRESSION: str = field(default="periodicdistance(x, y, z, x0, y0, z0)^2")
    _cv_force = field(default=None)
    restraining_force = field(default=None)
    _cv_names = field(default=[])
    
    def get_num_groups(self) -> int:
        """
        Get the number of groups in the CV.
        """
        num_groups = len(self.groups_atom_indices)
        return num_groups

    def make_cv_force(self) -> typing.Any:
        """
        Create the OpenMM force object for this CV.
        """
        external_force = openmm.CustomExternalForce(CV_EXPRESSION)
        pdb_file = openmm_app.PDBFile(self.ref_structure)
        external_force.addPerParticleParameter("x0")
        external_force.addPerParticleParameter("y0")
        external_force.addPerParticleParameter("z0")
        for i, index in enumerate(self.group_atom_indices):
            position = pdb_file.positions[index]
            external_force.addParticle(index, [position[0], position[1], position[2]])
        cv_force = openmm.CustomCVForce("sqrt(cv/num_atoms)")
        cv_force.addCollectiveVariable("cv", external_force)
        cv_force.addGlobalParameter("num_atoms", len(self.group_atom_indices))
        self._cv_force = cv_force
        return cv_force
        
    def add_cv_groups_and_variables(
            self,
            force: typing.Any
            ) -> None:
        """
        Add the groups and variables to the OpenMM force object
        for when this is being used just as a CV.
        """
        return

    def make_voronoi_cv_boundary_forces(
            self,
            me_val: float,
            neighbor_val: float
            ) -> tuple[openmm.CustomCVForce, openmm.CustomCVForce]:
        """
        Make the force objects to track proximity to current and adjacent anchors.
        """
        pdb_file = openmm_app.PDBFile(self.ref_structure)
        me_expr = f"(me_val - sqrt(cv/num_atoms))^2"
        me_external_force = openmm.CustomExternalForce(CV_EXPRESSION)
        me_external_force.addPerParticleParameter("x0")
        me_external_force.addPerParticleParameter("y0")
        me_external_force.addPerParticleParameter("z0")
        for i, index in enumerate(self.group_atom_indices):
            position = pdb_file.positions[index]
            me_external_force.addParticle(index, [position[0], position[1], position[2]])
        me_cv_force = openmm.CustomCVForce(me_expr)
        me_cv_force.addCollectiveVariable("cv", me_external_force)
        me_cv_force.addGlobalParameter("num_atoms", len(self.group_atom_indices))

        neighbor_expr = f"(neighbor_val - sqrt(cv/num_atoms))^2"
        neighbor_external_force = openmm.CustomExternalForce(CV_EXPRESSION)
        neighbor_external_force.addPerParticleParameter("x0")
        neighbor_external_force.addPerParticleParameter("y0")
        neighbor_external_force.addPerParticleParameter("z0")
        for i, index in enumerate(self.group_atom_indices):
            position = pdb_file.positions[index]
            neighbor_external_force.addParticle(index, [position[0], position[1], position[2]])
        neighbor_cv_force = openmm.CustomCVForce(neighbor_expr)
        neighbor_cv_force.addCollectiveVariable("cv", neighbor_external_force)
        neighbor_cv_force.addGlobalParameter("num_atoms", len(self.group_atom_indices))
        return me_cv_force, neighbor_cv_force

    def make_boundary_force_expr(
            self,
            index: int
            ) -> str:
        """
        Create the OpenMM force object for this CV.
        """
        boundary_expression = f"step(k_{index}*(sqrt(cv_{index}/num_atoms)-value_{index}))"
        self._cv_names.append(f"cv_{index}")
        return boundary_expression

    def make_boundary_force(
            self,
            full_expr: str
            ) -> openmm.CustomCVForce:
        """
        Create the OpenMM boundary force object for this CV.
        """
        force = openmm.CustomCVForce(full_expr)
        return force

    def add_boundary_groups_and_variables(
            self, 
            force: openmm.CustomCVForce,
            boundary_parameters: dict[str, unit.Quantity]
            ) -> None:
        """
        Add the groups and variables to the OpenMM force object
        for when this is being used just as a CV.
        """
        pdb_file = openmm_app.PDBFile(self.ref_structure)
        force.addGlobalParameter("num_atoms", len(self.group_atom_indices))
        for i, cv_name in enumerate(self._cv_names):
            external_force = openmm.CustomExternalForce(CV_EXPRESSION)
            external_force.addPerParticleParameter("x0")
            external_force.addPerParticleParameter("y0")
            external_force.addPerParticleParameter("z0")
            for index in self.group_atom_indices:
                position = pdb_file.positions[index]
                external_force.addParticle(index, [position[0], position[1], position[2]])

            force.addCollectiveVariable(cv_name, external_force)
            
        for key, value in boundary_parameters.items():
            force.addGlobalParameter(key, value)
        return

    def make_restraining_force(self) -> openmm.CustomCVForce:
        """
        Create the OpenMM restraining force object for this CV.
        """
        pdb_file = openmm_app.PDBFile(self.ref_structure)
        restraint_expression = f"0.5*k*sqrt(cv/num_atoms)-value)^2"
        self.restraining_force = openmm.CustomCVForce(restraint_expression)
        return self.restraining_force
    
    def make_voronoi_cv_restraining_forces(
            self,
            #toggle: float,
            #value: float,
            ) -> openmm.CustomCVForce:
        """
        Make the force objects to track proximity to current and adjacent anchors.
        """
        pdb_file = openmm_app.PDBFile(self.ref_structure)
        restraint_expression = f"sqrt(cv/num_atoms)"
        self.restraining_force = openmm.CustomCVForce(restraint_expression)
        external_force = openmm.CustomExternalForce(CV_EXPRESSION)
        external_force.addPerParticleParameter("x0")
        external_force.addPerParticleParameter("y0")
        external_force.addPerParticleParameter("z0")
        for index in self.group_atom_indices:
            position = pdb_file.positions[index]
            external_force.addParticle(index, [position[0], position[1], position[2]])
        restraint_force.addCollectiveVariable("cv", external_force)
        restraint_force.addGlobalParameter("num_atoms", len(self.group_atom_indices))
        return self.restraining_force
    
    def add_restraining_groups_and_variables(
            self, 
            k: unit.Quantity, 
            value: unit.Quantity
            ) -> None:
        """
        Add the groups and variables to the OpenMM force object
        for when this is being used as a restraint.
        """
        assert self.restraining_force is not None
        force = self.restraining_force
        external_force = openmm.CustomExternalForce(CV_EXPRESSION)
        external_force.addPerParticleParameter("x0")
        external_force.addPerParticleParameter("y0")
        external_force.addPerParticleParameter("z0")
        for index in self.group_atom_indices:
            position = pdb_file.positions[index]
            external_force.addParticle(index, [position[0], position[1], position[2]])
        force.addCollectiveVariable("cv", external_force)
        force.addGlobalParameter("num_atoms", len(self.group_atom_indices))
        force.addGlobalParameter("k", k)
        force.addGlobalParameter("value", value)
        return

@define
class Min_pair_dist_CV_openmm_interface:
    collective_variable: min_pair_dist_cv.Min_pair_dist_CV = field()
    group1_atom_indices: typing.List[int] = field()
    group2_atom_indices: typing.List[int] = field()
    CV_EXPRESSION: str = field(default="periodicdistance(x, y, z, x0, y0, z0)^2")
    _cv_force = field(default=None)
    restraining_force = field(default=None)
    _cv_names = field(default=[])
    
    def get_num_groups(self) -> int:
        """
        Get the number of groups in the CV.
        """
        num_groups = 2
        return num_groups

    def make_cv_force(self) -> typing.Any:
        """
        Create the OpenMM force object for this CV.
        """
        external_force = openmm.CustomExternalForce(CV_EXPRESSION)
        pdb_file = openmm_app.PDBFile(self.ref_structure)
        external_force.addPerParticleParameter("x0")
        external_force.addPerParticleParameter("y0")
        external_force.addPerParticleParameter("z0")
        for i, index in enumerate(self.group_atom_indices):
            position = pdb_file.positions[index]
            external_force.addParticle(index, [position[0], position[1], position[2]])
        cv_force = openmm.CustomCVForce("sqrt(cv/num_atoms)")
        cv_force.addCollectiveVariable("cv", external_force)
        cv_force.addGlobalParameter("num_atoms", len(self.group_atom_indices))
        self._cv_force = cv_force
        return cv_force
        
    def add_cv_groups_and_variables(
            self,
            force: typing.Any
            ) -> None:
        """
        Add the groups and variables to the OpenMM force object
        for when this is being used just as a CV.
        """
        return

    def make_voronoi_cv_boundary_forces(
            self,
            me_val: float,
            neighbor_val: float
            ) -> tuple[openmm.CustomCVForce, openmm.CustomCVForce]:
        """
        Make the force objects to track proximity to current and adjacent anchors.
        """
        pdb_file = openmm_app.PDBFile(self.ref_structure)
        me_expr = f"(me_val - sqrt(cv/num_atoms))^2"
        me_external_force = openmm.CustomExternalForce(CV_EXPRESSION)
        me_external_force.addPerParticleParameter("x0")
        me_external_force.addPerParticleParameter("y0")
        me_external_force.addPerParticleParameter("z0")
        for i, index in enumerate(self.group_atom_indices):
            position = pdb_file.positions[index]
            me_external_force.addParticle(index, [position[0], position[1], position[2]])
        me_cv_force = openmm.CustomCVForce(me_expr)
        me_cv_force.addCollectiveVariable("cv", me_external_force)
        me_cv_force.addGlobalParameter("num_atoms", len(self.group_atom_indices))

        neighbor_expr = f"(neighbor_val - sqrt(cv/num_atoms))^2"
        neighbor_external_force = openmm.CustomExternalForce(CV_EXPRESSION)
        neighbor_external_force.addPerParticleParameter("x0")
        neighbor_external_force.addPerParticleParameter("y0")
        neighbor_external_force.addPerParticleParameter("z0")
        for i, index in enumerate(self.group_atom_indices):
            position = pdb_file.positions[index]
            neighbor_external_force.addParticle(index, [position[0], position[1], position[2]])
        neighbor_cv_force = openmm.CustomCVForce(neighbor_expr)
        neighbor_cv_force.addCollectiveVariable("cv", neighbor_external_force)
        neighbor_cv_force.addGlobalParameter("num_atoms", len(self.group_atom_indices))
        return me_cv_force, neighbor_cv_force

    def make_boundary_force_expr(
            self,
            index: int
            ) -> str:
        """
        Create the OpenMM force object for this CV.
        """
        boundary_expression = f"step(k_{index}*(sqrt(cv_{index}/num_atoms)-value_{index}))"
        self._cv_names.append(f"cv_{index}")
        return boundary_expression

    def make_boundary_force(
            self,
            full_expr: str
            ) -> openmm.CustomCVForce:
        """
        Create the OpenMM boundary force object for this CV.
        """
        force = openmm.CustomCVForce(full_expr)
        return force

    def add_boundary_groups_and_variables(
            self, 
            force: openmm.CustomCVForce,
            boundary_parameters: dict[str, unit.Quantity]
            ) -> None:
        """
        Add the groups and variables to the OpenMM force object
        for when this is being used just as a CV.
        """
        pdb_file = openmm_app.PDBFile(self.ref_structure)
        force.addGlobalParameter("num_atoms", len(self.group_atom_indices))
        for i, cv_name in enumerate(self._cv_names):
            external_force = openmm.CustomExternalForce(CV_EXPRESSION)
            external_force.addPerParticleParameter("x0")
            external_force.addPerParticleParameter("y0")
            external_force.addPerParticleParameter("z0")
            for index in self.group_atom_indices:
                position = pdb_file.positions[index]
                external_force.addParticle(index, [position[0], position[1], position[2]])

            force.addCollectiveVariable(cv_name, external_force)
            
        for key, value in boundary_parameters.items():
            force.addGlobalParameter(key, value)
        return

    def make_restraining_force(self) -> openmm.CustomCVForce:
        """
        Create the OpenMM restraining force object for this CV.
        """
        pdb_file = openmm_app.PDBFile(self.ref_structure)
        restraint_expression = f"0.5*k*sqrt(cv/num_atoms)-value)^2"
        self.restraining_force = openmm.CustomCVForce(restraint_expression)
        return self.restraining_force
    
    def make_voronoi_cv_restraining_forces(
            self,
            #toggle: float,
            #value: float,
            ) -> openmm.CustomCVForce:
        """
        Make the force objects to track proximity to current and adjacent anchors.
        """
        pdb_file = openmm_app.PDBFile(self.ref_structure)
        restraint_expression = f"sqrt(cv/num_atoms)"
        self.restraining_force = openmm.CustomCVForce(restraint_expression)
        external_force = openmm.CustomExternalForce(CV_EXPRESSION)
        external_force.addPerParticleParameter("x0")
        external_force.addPerParticleParameter("y0")
        external_force.addPerParticleParameter("z0")
        for index in self.group_atom_indices:
            position = pdb_file.positions[index]
            external_force.addParticle(index, [position[0], position[1], position[2]])
        restraint_force.addCollectiveVariable("cv", external_force)
        restraint_force.addGlobalParameter("num_atoms", len(self.group_atom_indices))
        return self.restraining_force
    
    def add_restraining_groups_and_variables(
            self, 
            k: unit.Quantity, 
            value: unit.Quantity
            ) -> None:
        """
        Add the groups and variables to the OpenMM force object
        for when this is being used as a restraint.
        """
        assert self.restraining_force is not None
        force = self.restraining_force
        external_force = openmm.CustomExternalForce(CV_EXPRESSION)
        external_force.addPerParticleParameter("x0")
        external_force.addPerParticleParameter("y0")
        external_force.addPerParticleParameter("z0")
        for index in self.group_atom_indices:
            position = pdb_file.positions[index]
            external_force.addParticle(index, [position[0], position[1], position[2]])
        force.addCollectiveVariable("cv", external_force)
        force.addGlobalParameter("num_atoms", len(self.group_atom_indices))
        force.addGlobalParameter("k", k)
        force.addGlobalParameter("value", value)
        return

CV_openmm_interface = typing.Union[seekr_cv_openmm_interface.CV_openmm_interface,
                                   RMSD_no_align_CV_openmm_interface,
                                   Min_pair_dist_CV_openmm_interface]

def assign_cv_interface(
        collective_variable: cv_base.Base_CV,
        md_engine: engines_structures.Base_engine,
        root_directory: str) -> CV_openmm_interface | None:
    """
    Given this plugin's collective variable object, assign the appropriate 
    OpenMM interface.
    """
    is_found = False
    openmm_interface = None
    if collective_variable.type == "rmsd_no_align":
        group_atom_indices = md_engine.get_atom_indices_by_selection_name(
            collective_variable.group_name)
        ref_structure_full_path = os.path.join(
            root_directory, collective_variable.directory,
            collective_variable.ref_structure)
        return True, RMSD_no_align_CV_openmm_interface(
            collective_variable=collective_variable, 
            group_atom_indices=group_atom_indices,
            ref_structure=ref_structure_full_path)
    elif collective_variable.type == "min_pair_dist":
        group1_atom_indices = md_engine.get_atom_indices_by_selection_name(
            collective_variable.group1_name)
        group2_atom_indices = md_engine.get_atom_indices_by_selection_name(
            collective_variable.group2_name)
        return True, Min_pair_dist_CV_openmm_interface(
            collective_variable=collective_variable,
            group1_atom_indices=group1_atom_indices,
            group2_atom_indices=group2_atom_indices)
    return is_found, openmm_interface

