"""
Unit tests for seekrsamplingmethodsplugin XY collective variable.
"""

from __future__ import annotations

from unittest.mock import MagicMock

import mdtraj
import numpy as np
import pytest

import seekrsamplingmethodsplugin.modules.collective_variables.xy_cv as xy_cv
import seekrsamplingmethodsplugin.modules.collective_variables.openmm_interface \
    as plugin_cv_openmm


def _two_atom_traj(xyz: np.ndarray) -> mdtraj.Trajectory:
    top = mdtraj.Topology()
    chain = top.add_chain()
    res = top.add_residue("DUM", chain)
    top.add_atom("C", mdtraj.element.carbon, res)
    top.add_atom("C", mdtraj.element.carbon, res)
    return mdtraj.Trajectory(np.asarray(xyz, dtype=np.float32), top)


def _mock_model(group1_indices, group2_indices, g1_name="host", g2_name="guest"):
    mock_engine = MagicMock()
    mock_engine.get_atom_indices_by_selection_name.side_effect = (
        lambda name: group1_indices if name == g1_name else group2_indices
    )
    mock_model = MagicMock()
    mock_model.get_engine_by_scale_type.return_value = mock_engine
    return mock_model, mock_engine


def test_get_atom_group_names():
    cv = xy_cv.XY_CV(index=0, group1_name="host", group2_name="guest")
    assert cv.get_atom_group_names() == ["host", "guest"]


def test_get_mdtraj_cv_value():
    """rho = 5 for atoms at (0,0,0) and (3,4,0)."""
    cv = xy_cv.XY_CV(index=0, group1_name="host", group2_name="guest")
    traj = _two_atom_traj([[[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]]])
    mock_model, _ = _mock_model([0], [1])
    result = cv.get_mdtraj_cv_value(mock_model, traj, frame_index=0)
    assert result == pytest.approx(5.0, abs=1e-5)


def test_get_mdtraj_cv_value_zero():
    cv = xy_cv.XY_CV(index=0, group1_name="host", group2_name="guest")
    traj = _two_atom_traj([[[0.5, 0.5, 0.5], [0.5, 0.5, 0.5]]])
    mock_model, _ = _mock_model([0], [1])
    result = cv.get_mdtraj_cv_value(mock_model, traj, frame_index=0)
    assert result == pytest.approx(0.0, abs=1e-5)


def test_make_cv_object_default_name_and_bounds():
    cv_input = xy_cv.XY_CV_input()
    cv_input.group1_name = "host"
    cv_input.group2_name = "guest"
    cv_input.min_value = 0.0
    cv_input.max_value = 2.0
    cv = cv_input.make_cv_object(index=2, root_directory="/tmp")
    assert cv.index == 2
    assert cv.name == "xy_host_guest"
    assert cv.min_value == 0.0
    assert cv.max_value == 2.0
    assert cv.type == "xy"
    assert cv.unit == "nm"


def test_make_cv_object_custom_name():
    cv_input = xy_cv.XY_CV_input()
    cv_input.name = "xy_guest"
    cv_input.group1_name = "host"
    cv_input.group2_name = "guest"
    cv = cv_input.make_cv_object(index=0, root_directory="/tmp")
    assert cv.name == "xy_guest"


def test_openmm_cv_force_uses_pbc():
    cv = xy_cv.XY_CV(index=0, group1_name="host", group2_name="guest")
    iface = plugin_cv_openmm.XY_cv_openmm_interface(
        collective_variable=cv,
        group1_atom_indices=[0],
        group2_atom_indices=[1],
    )
    force = iface.make_cv_force()
    iface.add_cv_groups_and_variables(force)
    assert force.getEnergyFunction() == "sqrt((x2 - x1)^2 + (y2 - y1)^2)"
    assert force.usesPeriodicBoundaryConditions() is True


def test_assign_cv_interface_xy():
    cv = xy_cv.XY_CV(index=0, group1_name="host", group2_name="guest")
    mock_engine = MagicMock()
    mock_engine.get_atom_indices_by_selection_name.side_effect = (
        lambda name: [0] if name == "host" else [1]
    )
    found, iface = plugin_cv_openmm.assign_cv_interface(
        cv, mock_engine, root_directory="/tmp")
    assert found is True
    assert isinstance(iface, plugin_cv_openmm.XY_cv_openmm_interface)
    assert iface.group1_atom_indices == [0]
    assert iface.group2_atom_indices == [1]
