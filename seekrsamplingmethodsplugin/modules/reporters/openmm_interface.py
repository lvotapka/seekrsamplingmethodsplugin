"""
modules/reporters/openmm_interface.py

OpenMM interfaces for seekrsamplingmethodsplugin reporters.
"""

from __future__ import annotations

import csv
import typing

import numpy as np
from attrs import define, field
from openmm import unit

import seekr.modules.reporters.structures as reporters_structures
import seekr.modules.reporters.openmm_interface as reporters_openmm_interface

import seekrsamplingmethodsplugin.modules.reporters.structures \
    as plugin_reporters_structures
import seekrsamplingmethodsplugin.modules.reporters.features as plugin_features


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

    return is_found, openmm_interface
