"""
modules/reporters/openmm_interface.py

OpenMM interfaces for seekrsamplingmethodsplugin reporters.
"""

from __future__ import annotations

import csv
import os
import typing

import numpy as np
from attrs import define, field
import openmm.app as openmm_app
from openmm import unit

import seekr.modules.reporters.structures as reporters_structures
import seekr.modules.reporters.openmm_interface as reporters_openmm_interface

import seekrsamplingmethodsplugin.modules.reporters.structures \
    as plugin_reporters_structures
import seekrsamplingmethodsplugin.modules.reporters.features as plugin_features
import seekrsamplingmethodsplugin.modules.reporters.conditions as plugin_conditions


class FeatureOpenMMReporter:
    """
    Custom OpenMM reporter that writes geometric features to a CSV file.

    The ``append`` argument mirrors OpenMM's DCDReporter: when True, rows are
    appended to an existing file (no new header); when False, the file is
    overwritten and a header row is written.
    """

    def __init__(
            self,
            filename: str,
            reportInterval: int,
            features: list[plugin_features.Base_feature],
            append: bool = False,
            ) -> None:
        assert len(features) > 0, \
            "FeatureOpenMMReporter requires at least one feature."
        self._filename = filename
        self._reportInterval = reportInterval
        self._features = list(features)
        self._append = append
        self._masses: np.ndarray | None = None
        self._uses_pbc: bool | None = None
        self._column_names = [feat.column_name() for feat in self._features]

        mode = "a" if append else "w"
        self._out = open(filename, mode, newline="")
        self._writer = csv.writer(self._out)
        if not append:
            self._writer.writerow(["step", "time_ps", *self._column_names])
            self._out.flush()

    def describeNextReport(self, simulation):
        steps = self._reportInterval \
            - simulation.currentStep % self._reportInterval
        if self._uses_pbc is None:
            self._uses_pbc = simulation.system.usesPeriodicBoundaryConditions()
        return {
            "steps": steps,
            "periodic": bool(self._uses_pbc),
            "include": ["positions"],
        }

    def report(self, simulation, state) -> None:
        if self._masses is None:
            n_particles = simulation.system.getNumParticles()
            self._masses = np.array([
                simulation.system.getParticleMass(i).value_in_unit(unit.dalton)
                for i in range(n_particles)
            ], dtype=float)
        if self._uses_pbc is None:
            self._uses_pbc = simulation.system.usesPeriodicBoundaryConditions()

        positions = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        box_vectors = None
        if self._uses_pbc:
            box_vectors = np.array([
                v.value_in_unit(unit.nanometer)
                for v in state.getPeriodicBoxVectors()
            ], dtype=float)

        step = simulation.currentStep
        dt_ps = simulation.context.getIntegrator().getStepSize()\
            .value_in_unit(unit.picoseconds)
        time_ps = step * dt_ps

        values = [
            feat.compute(
                positions, box_vectors=box_vectors, masses=self._masses)
            for feat in self._features
        ]
        self._writer.writerow([step, time_ps, *values])
        self._out.flush()

    def __del__(self):
        try:
            if hasattr(self, "_out") and self._out is not None \
                    and not self._out.closed:
                self._out.close()
        except Exception:
            pass


def _topology_subset(
        topology: openmm_app.Topology,
        atom_indices: list[int],
        ) -> openmm_app.Topology:
    """
    Build a topology containing only the given atom indices (in that order).
    """
    keep_order = list(atom_indices)
    keep_set = set(keep_order)
    atoms_by_index = {atom.index: atom for atom in topology.atoms()}
    new_top = openmm_app.Topology()
    # One chain/residue bucket is enough for trajectory writing.
    new_chain = new_top.addChain()
    new_res = new_top.addResidue("SUB", new_chain)
    for idx in keep_order:
        assert idx in atoms_by_index, f"atom_subset index {idx} not in topology."
        atom = atoms_by_index[idx]
        new_top.addAtom(atom.name, atom.element, new_res)
    assert len(keep_set) == len(keep_order), \
        "atom_subset must not contain duplicate indices."
    return new_top


class ConditionalPositionOpenMMReporter:
    """
    Write positions to a temp trajectory while CV conditions hold.

    On condition exit the temp file is deleted (sojourn forgotten). On finalize,
    the temp file is promoted to the final path only if still active.
    """

    def __init__(
            self,
            filename: str,
            reportInterval: int,
            conditions: list[plugin_conditions.CV_threshold_condition],
            operator: str,
            extension: str,
            enforcePeriodicBox: bool,
            atom_subset: list[int] | None = None,
            sampling_method_interface=None,
            ) -> None:
        assert len(conditions) > 0, \
            "ConditionalPositionOpenMMReporter requires at least one condition."
        assert extension in ("dcd", "pdb"), \
            f"Unsupported extension {extension}."
        assert operator in ("and", "or"), f"Unsupported operator {operator}."
        assert sampling_method_interface is not None, \
            "sampling_method_interface is required."

        self._filepath = filename
        self._reportInterval = reportInterval
        self._conditions = list(conditions)
        self._operator = operator
        self._extension = extension
        self._enforcePeriodicBox = enforcePeriodicBox
        self._atom_subset = atom_subset
        self._sampling_method_interface = sampling_method_interface

        self._temp_path = filename + ".tmp"
        self._active = False
        self._n_frames = 0
        self._file = None
        self._dcd = None
        self._topology = None
        self._write_topology = None
        self._dt = None
        self._cv_indices: list[int] | None = None
        self._finalized = False

    def describeNextReport(self, simulation):
        steps = self._reportInterval \
            - simulation.currentStep % self._reportInterval
        return {
            "steps": steps,
            "periodic": self._enforcePeriodicBox,
            "include": ["positions"],
        }

    def _resolve_cv_indices(self) -> None:
        if self._cv_indices is not None:
            return
        model = self._sampling_method_interface._model
        indices = []
        for condition in self._conditions:
            cv = model.get_cv_by_name(condition.cv_name)
            assert cv is not None, \
                f"CV '{condition.cv_name}' not found in model."
            indices.append(cv.index)
        self._cv_indices = indices

    def _ensure_cv_forces(self) -> None:
        smi = self._sampling_method_interface
        smi.create_openmm_cv_list()
        smi.create_openmm_cv_object()

    def _conditions_satisfied(self, cv_values: list[float]) -> bool:
        results = []
        for condition, cv_index in zip(self._conditions, self._cv_indices):
            results.append(condition.is_satisfied(float(cv_values[cv_index])))
        if self._operator == "and":
            return all(results)
        return any(results)

    def _close_temp_writer(self) -> None:
        if self._extension == "pdb" and self._file is not None \
                and self._write_topology is not None and self._n_frames > 0:
            openmm_app.PDBFile.writeFooter(self._write_topology, self._file)
        if self._file is not None and not self._file.closed:
            self._file.close()
        self._file = None
        self._dcd = None

    def _discard_temp(self) -> None:
        self._close_temp_writer()
        if os.path.exists(self._temp_path):
            os.remove(self._temp_path)
        self._active = False
        self._n_frames = 0

    def _open_temp(self) -> None:
        self._discard_temp()
        assert self._write_topology is not None
        assert self._dt is not None
        if self._extension == "dcd":
            self._file = open(self._temp_path, "wb")
            self._dcd = openmm_app.DCDFile(
                self._file, self._write_topology, dt=self._dt,
                firstStep=0, interval=1, append=False)
        else:
            self._file = open(self._temp_path, "w")
            openmm_app.PDBFile.writeHeader(self._write_topology, self._file)
        self._active = True
        self._n_frames = 0

    def _positions_for_write(self, positions):
        if self._atom_subset is None:
            return positions
        # Quantity array indexing preserves units.
        return positions[self._atom_subset]

    def _write_frame(self, positions, box_vectors) -> None:
        assert self._active
        write_pos = self._positions_for_write(positions)
        if self._extension == "dcd":
            assert self._dcd is not None
            self._dcd.writeModel(
                write_pos, periodicBoxVectors=box_vectors)
        else:
            assert self._file is not None
            assert self._write_topology is not None
            openmm_app.PDBFile.writeModel(
                self._write_topology, write_pos, self._file,
                self._n_frames + 1)
        self._n_frames += 1

    def report(self, simulation, state) -> None:
        if self._topology is None:
            self._topology = simulation.topology
            self._dt = simulation.integrator.getStepSize()
            if self._atom_subset is None:
                self._write_topology = self._topology
            else:
                self._write_topology = _topology_subset(
                    self._topology, self._atom_subset)
            self._ensure_cv_forces()
            self._resolve_cv_indices()

        cv_values = self._sampling_method_interface.get_cv_values()
        if cv_values is None:
            return

        satisfied = self._conditions_satisfied(cv_values)
        positions = state.getPositions(asNumpy=True)
        box_vectors = state.getPeriodicBoxVectors()

        if satisfied:
            if not self._active:
                self._open_temp()
            self._write_frame(positions, box_vectors)
        else:
            if self._active:
                self._discard_temp()
        return

    def finalize(self) -> None:
        """
        Promote the temp sojourn to the final path if still active; otherwise
        discard. Safe to call multiple times.
        """
        if self._finalized:
            return
        self._finalized = True
        if self._active and self._n_frames > 0:
            self._close_temp_writer()
            os.replace(self._temp_path, self._filepath)
            self._active = False
            self._n_frames = 0
        else:
            self._discard_temp()
        return

    def __del__(self):
        try:
            self.finalize()
        except Exception:
            pass


@define
class Feature_reporter_openmm_interface(
        reporters_openmm_interface.Base_reporter_openmm_interface):
    """
    Interface for the Feature_reporter using the OpenMM engine.
    """
    _reporter: plugin_reporters_structures.Feature_reporter = field()
    _sampling_method_interface: typing.Any = field(default=None)

    def add_reporter_to_simulation(
            self,
            **kwargs: typing.Any
            ) -> None:
        """
        Add the feature reporter to the OpenMM simulation object.
        """
        simulation = kwargs["simulation"]
        filename = kwargs["filename"]
        openmm_reporter = FeatureOpenMMReporter(
            filename=filename,
            reportInterval=self._reporter.interval,
            features=self._reporter.features,
            append=False,
        )
        simulation.reporters.append(openmm_reporter)
        return

    def get_reporter_extension(self) -> str:
        """
        Get the file extension for the reporter.
        """
        return self._reporter.get_reporter_extension()

    def get_progress(self) -> float:
        """
        Feature reporters have no notion of partial progress.
        """
        return 1.0

    def add_sampling_method_interface(
            self,
            sampling_method_interface
            ) -> None:
        self._sampling_method_interface = sampling_method_interface
        return


@define
class Conditional_position_reporter_openmm_interface(
        reporters_openmm_interface.Base_reporter_openmm_interface):
    """
    Interface for the Conditional_position_reporter using OpenMM.
    """
    # Declare all fields explicitly (stable attrs order for construction).
    _reporter: plugin_reporters_structures.Conditional_position_reporter = field()
    _sampling_method_interface: typing.Any = field(default=None)

    def add_reporter_to_simulation(
            self,
            **kwargs: typing.Any
            ) -> None:
        """
        Add the conditional position reporter to the OpenMM simulation.
        """
        assert self._sampling_method_interface is not None, \
            "sampling_method_interface must be set before adding the reporter."
        simulation = kwargs["simulation"]
        filename = kwargs["filename"]
        toy_mode = kwargs.get("toy_mode", False)
        if toy_mode:
            enforce_periodic_box = False
        else:
            enforce_periodic_box = self._reporter.enforce_periodic_box
        openmm_reporter = ConditionalPositionOpenMMReporter(
            filename=filename,
            reportInterval=self._reporter.interval,
            conditions=self._reporter.conditions,
            operator=self._reporter.operator,
            extension=self._reporter.extension,
            enforcePeriodicBox=enforce_periodic_box,
            atom_subset=self._reporter.atom_subset,
            sampling_method_interface=self._sampling_method_interface,
        )
        simulation.reporters.append(openmm_reporter)
        return

    def get_reporter_extension(self) -> str:
        """
        Get the file extension for the reporter.
        """
        return self._reporter.get_reporter_extension()

    def get_progress(self) -> float:
        """
        Conditional position reporters have no notion of partial progress.
        """
        return 1.0

    def add_sampling_method_interface(
            self,
            sampling_method_interface
            ) -> None:
        self._sampling_method_interface = sampling_method_interface
        return


def assign_reporter_interface(
        reporter: reporters_structures.Base_reporter,
        ) -> tuple[bool, typing.Optional[
            reporters_openmm_interface.Base_reporter_openmm_interface]]:
    """
    Assign the OpenMM interface for plugin reporters.
    """
    is_found = False
    openmm_interface = None
    if reporter.type == "feature":
        is_found = True
        openmm_interface = Feature_reporter_openmm_interface(reporter)
    elif reporter.type == "conditional_position":
        is_found = True
        openmm_interface = Conditional_position_reporter_openmm_interface(
            reporter)

    return is_found, openmm_interface
