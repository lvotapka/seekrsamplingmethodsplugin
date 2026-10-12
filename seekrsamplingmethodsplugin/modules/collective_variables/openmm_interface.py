"""
modules/collective_variables/openmm_interface.py

The OpenMM interface for the collective variables defined in this
plugin.
"""

import os
import typing

from attrs import define, field
import openmm
import openmm.app as openmm_app
import openmm.unit as unit

import seekr.modules.collective_variables.base as cv_base 
import seekr.modules.collective_variables.openmm_interface \
    as seekr_cv_openmm_interface
import seekr.modules.engines.structures as engines_structures

import seekrsamplingmethodsplugin.modules.collective_variables\
    .rmsd_no_align_cv as rmsd_no_align_cv
import seekrsamplingmethodsplugin.modules.collective_variables\
    .min_pair_dist_cv as min_pair_dist_cv
import seekrsamplingmethodsplugin.modules.collective_variables\
    .structures as plugin_cv_structures
import seekrsamplingmethodsplugin.modules.collective_variables\
    .xy_cv as xy_cv
import seekrsamplingmethodsplugin.modules.collective_variables\
    .z_cv as z_cv

@define
class RMSD_no_align_CV_openmm_interface:
    collective_variable: rmsd_no_align_cv.RMSD_no_align_CV = field()
    group_atom_indices: typing.List[int] = field()
    ref_structure: str = field()
    CV_EXPRESSION: str = field(default="periodicdistance(x, y, z, x0, y0, z0)^2")
    _cv_force = field(default=None)
    restraining_force = field(default=None)
    _cv_names = field(factory=list)
    
    def get_num_groups(self) -> int:
        """
        Get the number of groups in the CV.
        """
        num_groups = len(self.group_atom_indices)
        return num_groups

    def make_cv_force(self) -> typing.Any:
        """
        Create the OpenMM force object for this CV.
        """
        external_force = openmm.CustomExternalForce(self.CV_EXPRESSION)
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
        me_external_force = openmm.CustomExternalForce(self.CV_EXPRESSION)
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
        neighbor_external_force = openmm.CustomExternalForce(self.CV_EXPRESSION)
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
            external_force = openmm.CustomExternalForce(self.CV_EXPRESSION)
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
        restraint_expression = f"0.5*k*(sqrt(cv/num_atoms)-value)^2"
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
        external_force = openmm.CustomExternalForce(self.CV_EXPRESSION)
        external_force.addPerParticleParameter("x0")
        external_force.addPerParticleParameter("y0")
        external_force.addPerParticleParameter("z0")
        for index in self.group_atom_indices:
            position = pdb_file.positions[index]
            external_force.addParticle(index, [position[0], position[1], position[2]])
        self.restraining_force.addCollectiveVariable("cv", external_force)
        self.restraining_force.addGlobalParameter("num_atoms", len(self.group_atom_indices))
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
        pdb_file = openmm_app.PDBFile(self.ref_structure)
        force = self.restraining_force
        external_force = openmm.CustomExternalForce(self.CV_EXPRESSION)
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
    CV_EXPRESSION: str = field(default="(1/r)^exponent")
    num_system_particles: int = field(default=0)
    # Set by assign_system_parameters(); distinct from CV.periodic (metadynamics grid).
    is_periodic_system: bool = field(default=False)
    _cv_force = field(default=None)
    restraining_force = field(default=None)
    _cv_names = field(factory=list)
    
    def get_num_groups(self) -> int:
        """
        Get the number of groups in the CV.
        """
        num_groups = 2
        return num_groups

    def _get_nonbonded_method(self) -> int:
        """
        Choose CustomNonbondedForce cutoff method from system periodicity.
        """
        if self.is_periodic_system:
            return openmm.CustomNonbondedForce.CutoffPeriodic
        return openmm.CustomNonbondedForce.CutoffNonPeriodic

    def _make_softmin_nonbonded_force(self) -> openmm.CustomNonbondedForce:
        """
        Build the CustomNonbondedForce that accumulates softmin pair terms.
        """
        assert self.num_system_particles > 0, \
            "assign_system_parameters() must be called before building "\
            "softmin nonbonded forces."
        cutoff_distance = self.collective_variable.cutoff_distance
        exponent = self.collective_variable.exponent
        nonbonded_force = openmm.CustomNonbondedForce(self.CV_EXPRESSION)
        nonbonded_force.addGlobalParameter("exponent", exponent)
        nonbonded_force.addInteractionGroup(self.group1_atom_indices,
                                            self.group2_atom_indices)
        nonbonded_force.setUseSwitchingFunction(False)
        nonbonded_force.setNonbondedMethod(self._get_nonbonded_method())
        nonbonded_force.setCutoffDistance(cutoff_distance*unit.nanometer)
        nonbonded_force.setUseLongRangeCorrection(False)
        for _ in range(self.num_system_particles):
            nonbonded_force.addParticle([])
        return nonbonded_force

    def _softmin_cv_expression(self, cv_name: str = "cv") -> str:
        """
        Softmin transform of the summed pair terms into a scalar CV value.
        Currently min(r) (or cutoff).
        """
        cutoff_distance = self.collective_variable.cutoff_distance
        return f"({cv_name} + (1/{cutoff_distance})^exponent)^(-1/exponent)"

    def make_cv_force(self) -> typing.Any:
        """
        Create the OpenMM force object for this CV.
        """
        exponent = self.collective_variable.exponent
        cv_force = openmm.CustomCVForce(self._softmin_cv_expression())
        cv_force.addCollectiveVariable("cv", self._make_softmin_nonbonded_force())
        cv_force.addGlobalParameter("exponent", exponent)
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
        softmin = self._softmin_cv_expression()
        me_cv_force = openmm.CustomCVForce(f"(me_val - {softmin})^2")
        me_cv_force.addCollectiveVariable("cv", self._make_softmin_nonbonded_force())

        neighbor_cv_force = openmm.CustomCVForce(f"(neighbor_val - {softmin})^2")
        neighbor_cv_force.addCollectiveVariable(
            "cv", self._make_softmin_nonbonded_force())

        return me_cv_force, neighbor_cv_force

    def make_boundary_force_expr(
            self,
            index: int
            ) -> str:
        """
        Create the OpenMM force object for this CV.
        """
        softmin = self._softmin_cv_expression(cv_name=f"cv_{index}")
        boundary_expression = f"step(k_{index}*({softmin}-value_{index}))"
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
        force.addGlobalParameter("exponent", self.collective_variable.exponent)
        for cv_name in self._cv_names:
            force.addCollectiveVariable(cv_name, self._make_softmin_nonbonded_force())
        for key, value in boundary_parameters.items():
            force.addGlobalParameter(key, value)
        return

    def make_restraining_force(self) -> openmm.CustomCVForce:
        """
        Create the OpenMM restraining force object for this CV.
        """
        restraint_expression = f"0.5*k*({self._softmin_cv_expression()}-value)^2"
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
        self.restraining_force = openmm.CustomCVForce(
            self._softmin_cv_expression())
        self.restraining_force.addCollectiveVariable(
            "cv", self._make_softmin_nonbonded_force())
        self.restraining_force.addGlobalParameter(
            "exponent", self.collective_variable.exponent)
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
        force.addCollectiveVariable("cv", self._make_softmin_nonbonded_force())
        force.addGlobalParameter("exponent", self.collective_variable.exponent)
        force.addGlobalParameter("k", k)
        force.addGlobalParameter("value", value)
        return

    def assign_system_parameters(
            self,
            openmm_system: openmm.System
            ) -> None:
        """
        Cache system size and periodicity for CustomNonbondedForce setup.

        is_periodic_system is about the MD box (PBC), not CV.periodic used for
        metadynamics grid wrapping of CV values like dihedrals.
        """
        self.num_system_particles = openmm_system.getNumParticles()
        self.is_periodic_system = openmm_system.usesPeriodicBoundaryConditions()
        return

@define
class Z_cv_openmm_interface:
    collective_variable: z_cv.Z_CV = field()
    group1_atom_indices: typing.List[int] = field()
    group2_atom_indices: typing.List[int] = field()
    NUM_GROUPS: int = field(default=2)
    CV_EXPRESSION: str = field(default="z2 - z1")
    BOUNDARY_EXPRESSION: str = field(default="step(k*(z2 - z1 - value))")
    RESTRAINING_EXPRESSION: str = field(default="0.5*k*(z2 - z1 - value)^2")
    _cv_force = field(default=None)
    restraining_force = field(default=None)
    megroup1 = field(default=0)
    megroup2 = field(default=1)

    def get_num_groups(self) -> int:
        """
        Get the number of groups in the CV.
        """
        return self.NUM_GROUPS

    def make_cv_force(self) -> openmm.CustomCentroidBondForce:
        """
        Create the OpenMM force object for this CV.
        """
        self._cv_force = openmm.CustomCentroidBondForce(
            self.NUM_GROUPS, self.CV_EXPRESSION)
        return self._cv_force

    def add_cv_groups_and_variables(
            self,
            force: openmm.CustomCentroidBondForce
            ) -> None:
        """
        Add the groups and variables to the OpenMM force object
        for when this is being used just as a CV.
        """
        group1 = force.addGroup(self.group1_atom_indices)
        group2 = force.addGroup(self.group2_atom_indices)
        group_list = [group1, group2]
        force.addBond(group_list, [])
        # NOTE: without this, the CVs will sometimes suffer from
        # jumps caused by periodic boundary affects.
        force.setUsesPeriodicBoundaryConditions(True)
        return

    def make_voronoi_cv_boundary_forces(
            self,
            me_val: float,
            neighbor_val: float
            ) -> tuple[openmm.CustomCentroidBondForce, openmm.CustomCentroidBondForce]:
        """
        Make the force objects to track proximity to current and adjacent anchors.
        """
        me_expr = f"(me_val - {self.CV_EXPRESSION})^2"
        me_force = openmm.CustomCentroidBondForce(self.NUM_GROUPS, me_expr)
        megroup1 = me_force.addGroup(self.group1_atom_indices)
        megroup2 = me_force.addGroup(self.group2_atom_indices)
        me_force.addPerBondParameter("me_val")
        me_force.addBond([megroup1, megroup2], [me_val])

        neighbor_expr = f"(neighbor_val - {self.CV_EXPRESSION})^2"
        neighbor_force = openmm.CustomCentroidBondForce(self.NUM_GROUPS, neighbor_expr)
        neighborgroup1 = neighbor_force.addGroup(self.group1_atom_indices)
        neighborgroup2 = neighbor_force.addGroup(self.group2_atom_indices)
        neighbor_force.addPerBondParameter("neighbor_val")
        neighbor_force.addBond([neighborgroup1, neighborgroup2], [neighbor_val])
        return me_force, neighbor_force

    def make_boundary_force_expr(
            self,
            index: int
            ) -> openmm.CustomCentroidBondForce:
        """
        Create the OpenMM force object for this CV.
        """
        boundary_expression = f"step(k_{index}*("+self.CV_EXPRESSION+f"-value_{index}))"
        return boundary_expression

    def make_boundary_force(
            self,
            full_expr: str
            ) -> openmm.CustomCentroidBondForce:
        """
        Create the OpenMM boundary force object for this CV.
        """
        num_groups = self.get_num_groups()
        force = openmm.CustomCentroidBondForce(num_groups, full_expr)
        return force

    def add_boundary_groups_and_variables(
            self, 
            force: openmm.CustomCentroidBondForce,
            boundary_parameters: dict[str, unit.Quantity]
            ) -> None:
        """
        Add the groups and variables to the OpenMM force object
        for when this is being used just as a CV.
        """
        group_list = []
        group1 = force.addGroup(self.group1_atom_indices)
        group_list.append(group1)
        group2 = force.addGroup(self.group2_atom_indices)
        group_list.append(group2)
        variable_values_list = []
        for key, value in boundary_parameters.items():
            force.addPerBondParameter(key)
            variable_values_list.append(value)
        force.addBond(group_list, variable_values_list)
        return

    def make_restraining_force(self):
        """
        Create the OpenMM restraining force object for this CV.
        """
        self.restraining_force = openmm.CustomCentroidBondForce(
            self.NUM_GROUPS, self.RESTRAINING_EXPRESSION)
        return self.restraining_force
    
    def make_voronoi_cv_restraining_forces(
            self,
            #toggle: float,
            #value: float,
            ) -> openmm.CustomCentroidBondForce:
        """
        Make the force objects to track proximity to current and adjacent anchors.
        """
        expr = self.CV_EXPRESSION
        force = openmm.CustomCentroidBondForce(self.NUM_GROUPS, expr)
        self.megroup1 = force.addGroup(self.group1_atom_indices)
        self.megroup2 = force.addGroup(self.group2_atom_indices)
        force.addBond([self.megroup1, self.megroup2], [])
        self.restraining_force = force
        return force
    
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
        group_list = []
        group1 = force.addGroup(self.group1_atom_indices)
        group_list.append(group1)
        group2 = force.addGroup(self.group2_atom_indices)
        group_list.append(group2)
        force.addPerBondParameter("k")
        force.addPerBondParameter("value")
        variable_values_list = [k, value]
        force.addBond(group_list, variable_values_list)
        return

@define
class XY_cv_openmm_interface:
    collective_variable: xy_cv.XY_CV = field()
    group1_atom_indices: typing.List[int] = field()
    group2_atom_indices: typing.List[int] = field()
    NUM_GROUPS: int = field(default=2)
    CV_EXPRESSION: str = field(default="sqrt((x2 - x1)^2 + (y2 - y1)^2)")
    BOUNDARY_EXPRESSION: str = field(default="step(k*(sqrt((x2 - x1)^2 + (y2 - y1)^2) - value))")
    RESTRAINING_EXPRESSION: str = field(default="0.5*k*(sqrt((x2 - x1)^2 + (y2 - y1)^2) - value)^2")
    _cv_force = field(default=None)
    restraining_force = field(default=None)
    megroup1 = field(default=0)
    megroup2 = field(default=1)

    def get_num_groups(self) -> int:
        """
        Get the number of groups in the CV.
        """
        return self.NUM_GROUPS

    def make_cv_force(self) -> openmm.CustomCentroidBondForce:
        """
        Create the OpenMM force object for this CV.
        """
        self._cv_force = openmm.CustomCentroidBondForce(
            self.NUM_GROUPS, self.CV_EXPRESSION)
        return self._cv_force

    def add_cv_groups_and_variables(
            self,
            force: openmm.CustomCentroidBondForce
            ) -> None:
        """
        Add the groups and variables to the OpenMM force object
        for when this is being used just as a CV.
        """
        group1 = force.addGroup(self.group1_atom_indices)
        group2 = force.addGroup(self.group2_atom_indices)
        group_list = [group1, group2]
        force.addBond(group_list, [])
        force.setUsesPeriodicBoundaryConditions(True)
        return

    def make_voronoi_cv_boundary_forces(
            self,
            me_val: float,
            neighbor_val: float
            ) -> tuple[openmm.CustomCentroidBondForce, openmm.CustomCentroidBondForce]:
        """
        Make the force objects to track proximity to current and adjacent anchors.
        """
        me_expr = f"(me_val - {self.CV_EXPRESSION})^2"
        me_force = openmm.CustomCentroidBondForce(self.NUM_GROUPS, me_expr)
        megroup1 = me_force.addGroup(self.group1_atom_indices)
        megroup2 = me_force.addGroup(self.group2_atom_indices)
        me_force.addPerBondParameter("me_val")
        me_force.addBond([megroup1, megroup2], [me_val])

        neighbor_expr = f"(neighbor_val - {self.CV_EXPRESSION})^2"
        neighbor_force = openmm.CustomCentroidBondForce(self.NUM_GROUPS, neighbor_expr)
        neighborgroup1 = neighbor_force.addGroup(self.group1_atom_indices)
        neighborgroup2 = neighbor_force.addGroup(self.group2_atom_indices)
        neighbor_force.addPerBondParameter("neighbor_val")
        neighbor_force.addBond([neighborgroup1, neighborgroup2], [neighbor_val])
        return me_force, neighbor_force

    def make_boundary_force_expr(
            self,
            index: int
            ) -> openmm.CustomCentroidBondForce:
        """
        Create the OpenMM force object for this CV.
        """
        boundary_expression = f"step(k_{index}*("+self.CV_EXPRESSION+f"-value_{index}))"
        return boundary_expression

    def make_boundary_force(
            self,
            full_expr: str
            ) -> openmm.CustomCentroidBondForce:
        """
        Create the OpenMM boundary force object for this CV.
        """
        num_groups = self.get_num_groups()
        force = openmm.CustomCentroidBondForce(num_groups, full_expr)
        return force

    def add_boundary_groups_and_variables(
            self, 
            force: openmm.CustomCentroidBondForce,
            boundary_parameters: dict[str, unit.Quantity]
            ) -> None:
        """
        Add the groups and variables to the OpenMM force object
        for when this is being used just as a CV.
        """
        group_list = []
        group1 = force.addGroup(self.group1_atom_indices)
        group_list.append(group1)
        group2 = force.addGroup(self.group2_atom_indices)
        group_list.append(group2)
        variable_values_list = []
        for key, value in boundary_parameters.items():
            force.addPerBondParameter(key)
            variable_values_list.append(value)
        force.addBond(group_list, variable_values_list)
        return

    def make_restraining_force(self):
        """
        Create the OpenMM restraining force object for this CV.
        """
        self.restraining_force = openmm.CustomCentroidBondForce(
            self.NUM_GROUPS, self.RESTRAINING_EXPRESSION)
        return self.restraining_force
    
    def make_voronoi_cv_restraining_forces(
            self,
            #toggle: float,
            #value: float,
            ) -> openmm.CustomCentroidBondForce:
        """
        Make the force objects to track proximity to current and adjacent anchors.
        """
        expr = self.CV_EXPRESSION
        force = openmm.CustomCentroidBondForce(self.NUM_GROUPS, expr)
        self.megroup1 = force.addGroup(self.group1_atom_indices)
        self.megroup2 = force.addGroup(self.group2_atom_indices)
        force.addBond([self.megroup1, self.megroup2], [])
        self.restraining_force = force
        return force
    
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
        group_list = []
        group1 = force.addGroup(self.group1_atom_indices)
        group_list.append(group1)
        group2 = force.addGroup(self.group2_atom_indices)
        group_list.append(group2)
        force.addPerBondParameter("k")
        force.addPerBondParameter("value")
        variable_values_list = [k, value]
        force.addBond(group_list, variable_values_list)
        return

CV_openmm_interface = typing.Union[seekr_cv_openmm_interface.CV_openmm_interface,
                                   RMSD_no_align_CV_openmm_interface,
                                   Min_pair_dist_CV_openmm_interface,
                                   XY_cv_openmm_interface,
                                   Z_cv_openmm_interface]

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
    elif collective_variable.type == "xy":
        group1_atom_indices = md_engine.get_atom_indices_by_selection_name(
            collective_variable.group1_name)
        group2_atom_indices = md_engine.get_atom_indices_by_selection_name(
            collective_variable.group2_name)
        return True, XY_cv_openmm_interface(
            collective_variable=collective_variable,
            group1_atom_indices=group1_atom_indices,
            group2_atom_indices=group2_atom_indices)
    elif collective_variable.type == "z":
        group1_atom_indices = md_engine.get_atom_indices_by_selection_name(
            collective_variable.group1_name)
        group2_atom_indices = md_engine.get_atom_indices_by_selection_name(
            collective_variable.group2_name)
        return True, Z_cv_openmm_interface(
            collective_variable=collective_variable,
            group1_atom_indices=group1_atom_indices,
            group2_atom_indices=group2_atom_indices)
    return is_found, openmm_interface

