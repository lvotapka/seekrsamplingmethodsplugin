"""
The openmm interface for a Random Accelerated Molecular Dynamics (RAMD) 
sampling method for this plugin.
"""

import os
import typing

from attrs import define
import numpy as np
import mdtraj
import openmm
import openmm.app as openmm_app
import openmm.unit as unit

import seekr.modules.structures as structures
import seekr.modules.sampling_methods.structures \
    as sampling_methods_structures
import seekr.modules.sampling_methods.openmm_interface \
    as sampling_methods_openmm_interface
import seekr.modules.collective_variables.openmm_interface \
    as collective_variables_openmm_interface

import seekrsamplingmethodsplugin.modules.sampling_methods.structures \
    as plugin_sampl_meth_struct

kcal_per_mole_per_angstrom = unit.kilocalories_per_mole / unit.angstrom
RAMD_FORCE_EXPRESSION \
    = "-ramd_force_magnitude*(ramd_f_x*x1 + ramd_f_y*y1 + ramd_f_z*z1)"
RAMD_FORCE_EXPRESSION_RELATIVE \
    = "-ramd_force_magnitude*(ramd_f_x*(x1-x2) + ramd_f_y*(y1-y2) "\
        + " + ramd_f_z*(z1-z2))"

@define
class Random_accelerated_md_openmm_interface(
        sampling_methods_openmm_interface.Base_md_openmm_interface):
    """
    Interface for the random accelerated molecular dynamics (RAMD) MD sampling 
    method using the OpenMM engine
    """
    _stage_index: int = 0
    sampling_method: plugin_sampl_meth_struct\
        .Random_accelerated_md_sampling_method | None = None
    _system_obj: openmm.System | None = None
    _simulation_obj: openmm_app.Simulation | None = None
    _steps_between_evals: int = 0
    _random_seed: int = 0
    _model: structures.Seekr_model | None = None
    _directory: str = ""
    _openmm_cv_list: list[typing.Any] | None = None
    _cv_interface_list: list[collective_variables_openmm_interface\
                             .CV_openmm_interface] | None = None
    _openmm_cv_object: typing.Any | None = None
    _force_vector: np.ndarray | None = None
    _random_vector: np.ndarray | None = None
    _start_force_magnitude: float = 0.0
    _current_force_magnitude: float = 0.0
    # Minimum ligand (or lig-rec) COM travel between updates, in nm.
    _cutoff_distance: float = 0.0
    _ramp_start_step: int = 0.0
    _ramp_force_increase_per_step: float = 0.0
    _old_lig_com: np.ndarray | None = None
    _lig_masses: np.ndarray | None = None
    _rec_masses: np.ndarray | None = None
    _force_object: typing.Any = None
    _receptor_atom_indices: list[int] | None = None
    _ligand_atom_indices: list[int] | None = None
    _group_index_lig: int = 0
    _group_index_rec: int = 0
    _bond_index: int = 0

    def get_evaluation_interval(self) -> int:
        """
        RAMD force retargeting interval (must match completion eval_interval
        when both are nonzero; SEEKR asserts that in the OpenMM runner).
        """
        assert self.sampling_method is not None
        return self.sampling_method.steps_per_update

    def _set_new_force_direction(self) -> None:
        """
        Sample a new unit-vector force direction and apply it to the context.
        """
        assert self._simulation_obj is not None
        f_gauss_x = np.random.normal()
        f_gauss_y = np.random.normal()
        f_gauss_z = np.random.normal()
        f_gauss_len = np.linalg.norm([f_gauss_x, f_gauss_y, f_gauss_z])
        force_vector_unscaled = np.array(
            [f_gauss_x / f_gauss_len, f_gauss_y / f_gauss_len,
             f_gauss_z / f_gauss_len])
        self._random_vector = force_vector_unscaled
        self._force_vector = force_vector_unscaled * self._current_force_magnitude
        self._simulation_obj.context.setParameter(
            "ramd_f_x", force_vector_unscaled[0])
        self._simulation_obj.context.setParameter(
            "ramd_f_y", force_vector_unscaled[1])
        self._simulation_obj.context.setParameter(
            "ramd_f_z", force_vector_unscaled[2])
        return

    def _compute_lig_com_baseline(self, positions: np.ndarray) -> np.ndarray:
        """
        Ligand COM, or ligand-minus-receptor COM vector, in nm.
        """
        ligand_positions = positions[self._ligand_atom_indices]
        if self._lig_masses is not None:
            ligand_com = np.average(
                ligand_positions, axis=0, weights=self._lig_masses)
        else:
            ligand_com = np.mean(ligand_positions, axis=0)
        if self._receptor_atom_indices is None \
                or len(self._receptor_atom_indices) == 0:
            return ligand_com
        receptor_positions = positions[self._receptor_atom_indices]
        if self._rec_masses is not None:
            receptor_com = np.average(
                receptor_positions, axis=0, weights=self._rec_masses)
        else:
            receptor_com = np.mean(receptor_positions, axis=0)
        return ligand_com - receptor_com
    
    def pre_simulation_preparations(
            self, 
            **kwargs: typing.Any
            ) -> None:
        """
        Perform any necessary preparations before the simulation starts.
        """
        system = self._system_obj
        assert system is not None, \
            "The system object must be assigned before calling this method."
        model: structures.Seekr_model = kwargs["model"]
        self._model = model
        self._directory: str = kwargs["directory"]
        md_engine = self._model.get_engine_by_scale_type("molecular_dynamics")
        assert md_engine is not None
        assert md_engine.type == "openmm", \
            "Only OpenMM engines can be used in the openmm_interface module."
        
        lig_selection = None
        rec_selection = None
        for selection in md_engine.unpartitioned_settings.selections:
            if selection.name == self.sampling_method.ligand_selection_name:
                lig_selection = selection
            elif selection.name == self.sampling_method.receptor_selection_name:
                    rec_selection = selection

        assert lig_selection is not None, \
            f"The ligand selection {self.sampling_method.ligand_selection_name} "\
            "was not found"
        if self.sampling_method.receptor_selection_name is not None:
            assert rec_selection is not None, \
            f"The receptor selection {self.sampling_method.receptor_selection_name} "\
            "was not found"
        lig_indices = lig_selection.atom_indices
        if self.sampling_method.receptor_selection_name is None:
            rec_indices = None
            
        else:
            rec_indices = rec_selection.atom_indices

        self._steps_between_evals = self.sampling_method.steps_per_update
        self._current_force_magnitude = self.sampling_method.force_magnitude
        self._start_force_magnitude = self.sampling_method.force_magnitude
        self._cutoff_distance = self.sampling_method.cutoff_distance
        self._ramp_start_step = self.sampling_method.ramp_start_step
        self._ramp_force_increase_per_step \
            = self.sampling_method.ramp_force_increase_per_step

        freeGroups = set(range(32)) \
            - set(myforce.getForceGroup() \
                  for myforce in system.getForces())
        if len(freeGroups) == 0:
            raise RuntimeError(
                "Cannot assign a force group to the steered MD stage. "\
                "The maximum number (32) of the force groups is "\
                "already used.")
        force_group = max(freeGroups)
        if rec_indices is None \
                or len(rec_indices) == 0:
            rec_indices = None
            force = openmm.CustomCentroidBondForce(
                1, RAMD_FORCE_EXPRESSION)
            force.setForceGroup(force_group)
            self._group_index_lig = force.addGroup(lig_indices)
            force.addGlobalParameter(
                "ramd_force_magnitude", 0.0*kcal_per_mole_per_angstrom)
            force.addGlobalParameter("ramd_f_x", 0.0)
            force.addGlobalParameter("ramd_f_y", 0.0)
            force.addGlobalParameter("ramd_f_z", 0.0)
            self._bond_index = force.addBond([self._group_index_lig])
            

        else:
            force = openmm.CustomCentroidBondForce(
                2, RAMD_FORCE_EXPRESSION_RELATIVE)
            force.setForceGroup(force_group)
            self._group_index_lig = force.addGroup(lig_indices)
            self._group_index_rec = force.addGroup(rec_indices)
            force.addGlobalParameter(
                "ramd_force_magnitude", 0.0*kcal_per_mole_per_angstrom)
            force.addGlobalParameter("ramd_f_x", 0.0)
            force.addGlobalParameter("ramd_f_y", 0.0)
            force.addGlobalParameter("ramd_f_z", 0.0)
            self._bond_index = force.addBond([self._group_index_lig, 
                               self._group_index_rec])
        
        self._force_vector = np.array([0.0, 0.0, 0.0])
        self._force_object = force
        system.addForce(force)
        assert self._simulation_obj is not None
        self._simulation_obj.context.reinitialize(preserveState=True)
        self._simulation_obj.context.setParameter(
            "ramd_force_magnitude", self._current_force_magnitude\
                *kcal_per_mole_per_angstrom)
        
        # compute the masses for ligand and receptor
        if rec_indices is None or len(rec_indices) == 0:
            self._rec_masses = None
            self._receptor_atom_indices = None
        else:
            rec_masses = []
            for index in rec_indices:
                rec_masses.append(system.getParticleMass(index)\
                                  .value_in_unit(unit.dalton))
            self._rec_masses = np.array(rec_masses)
            self._receptor_atom_indices = rec_indices
        lig_masses = []
        for index in lig_indices:
            lig_masses.append(system.getParticleMass(index)\
                              .value_in_unit(unit.dalton))
        self._lig_masses = np.array(lig_masses)
        self._ligand_atom_indices = lig_indices

        # Seed direction RNG (engine already called set_random_seed).
        np.random.seed(self._random_seed)
        # Apply a nonzero force direction before any MD steps, and record the
        # starting COM baseline (same convention as openmm_ramd RAMD_start).
        self._set_new_force_direction()
        state = self._simulation_obj.context.getState(getPositions=True)
        positions = state.getPositions(asNumpy=True)\
            .value_in_unit(unit.nanometer)
        self._old_lig_com = self._compute_lig_com_baseline(positions)
        return
    
    def post_simulation_procedures(
            self, 
            **kwargs: typing.Any
            ) -> None:
        """
        Perform any necessary procedures after the simulation ends.
        """
        swarm_id = kwargs["swarm_id"]
        step_count = int(self._simulation_obj.context.getStepCount())
        # Write the RAMD simulation length
        output_basename = f"ramd_swarm_{swarm_id}.out"
        output_filename = os.path.join(self._directory, output_basename)
        block_str = \
              f"RAMD Simulation Ending\n" \
            + f"step_count: {step_count}\n" 
        with open(output_filename, "w") as f:
            f.write(block_str)
        return
    
    def perform_eval(
            self, 
            _restart_index: int, 
            _swarm_id: int
            ) -> bool:
        """
        Perform an evaluation at a certain point in the simulation.
        """
        assert self._model is not None
        assert self._simulation_obj is not None
        if self._force_vector is None \
                or np.linalg.norm(self._force_vector) == 0.0:
            change_force = True
        else:
            change_force = False

        step_count = self._simulation_obj.context.getStepCount()
        if (step_count > self._ramp_start_step) \
                and self._ramp_force_increase_per_step > 0.0:
            # Modify force magnitude if ramping is being used
            new_force_magnitude = self._start_force_magnitude \
                + self._ramp_force_increase_per_step \
                * (step_count - self._ramp_start_step)
            old_force_magnitude = self._current_force_magnitude
            if old_force_magnitude == 0.0:
                assert self._random_vector is not None
                self._force_vector = new_force_magnitude \
                    * self._random_vector
            else:
                self._force_vector \
                    *= (new_force_magnitude/old_force_magnitude)
            self._simulation_obj.context.setParameter(
                "ramd_force_magnitude", new_force_magnitude*kcal_per_mole_per_angstrom)
            self._current_force_magnitude = new_force_magnitude

        state = self._simulation_obj.context.getState(getPositions=True)
        positions = state.getPositions(asNumpy=True)\
            .value_in_unit(unit.nanometer)
        lig_com_baseline = self._compute_lig_com_baseline(positions)
        if self._old_lig_com is not None:
            distance_moved = np.linalg.norm(
                lig_com_baseline - self._old_lig_com)
            if distance_moved <= self._cutoff_distance:
                change_force = True
        self._old_lig_com = lig_com_baseline

        if change_force:
            self._set_new_force_direction()
        
        return False
    
def assign_sampling_method_interface(
        stage_index: int,
        sampling_method: sampling_methods_structures\
            .Base_md_sampling_method, 
        ) -> tuple[bool, typing.Optional[
            sampling_methods_openmm_interface\
                .Base_md_openmm_interface]]:
    """
    Given a sampling_method and an md_engine, return the 
    corresponding interface
    """
    if sampling_method.type == "random_accelerated":
        return True, Random_accelerated_md_openmm_interface(
            stage_index, sampling_method=sampling_method)
        
    return False, None
