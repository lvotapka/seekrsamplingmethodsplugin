"""
create_seekr_input.py

Create SEEKR input objects for use in tests. The two public creator 
functions mirror the signatures of make_input.py so they can be 
the authoritative reference for all test scenarios.
"""

import os

import numpy as np

import seekr
import seekr.modules.base as base
import seekr.modules.engines.structures as engines_structures
import seekr.modules.completion_criteria.structures as compl_crit_struct
import seekr.modules.engines.openmm.preparer as openmm_preparer
import seekr.modules.engines.browndye.preparer as browndye_preparer
import seekr.modules.engines.sda.preparer as sda_preparer
import seekr.modules.sampling_methods.preparer as sampling_preparer
import seekr.modules.engines.preparer as engines_preparer
import seekr.modules.reporters.preparer as reporters_preparer
import seekr.modules.collective_variables.distance_cv as distance_cv
import seekr.modules.collective_variables.external_cv as external_cv
import seekr.modules.scales.molecular_dynamics.preparer as md_preparer
import seekr.modules.scales.brownian_dynamics.preparer as bd_preparer
import seekr.modules.anchor_schemes.preparer as anchor_scheme_preparer
import seekr.modules.scales.logistic.reassign_anchors_logistic_stage \
    as reassign_anchors_logistic_stage
import seekr.modules.scales.logistic.anchor_starting_structures_from_reporter_logistic_stage \
    as anchor_starting_structures_from_reporter_logistic_stage
import seekr.modules.preparer as preparer

#import seekrsamplingmethodsplugin.modules.something.something...


TEST_DIRECTORY = os.path.dirname(__file__)
SEEKR_TOOLS_DIRECTORY = os.path.dirname(seekr.__file__)
DATA_DIRECTORY = os.path.abspath(os.path.join(SEEKR_TOOLS_DIRECTORY, "data"))
HOSTGUEST_DIRECTORY = os.path.join(DATA_DIRECTORY, "hostguest_files")


# ---------------------------------------------------------------------------
# Topology / parameter helpers (kept for future ff-type fixtures)
# ---------------------------------------------------------------------------

def assign_amber_params(
        unpartitioned: engines_preparer.MD_unpartitioned_settings_input,
        prmtop_filename: str,
        pdb_filename: str,
        ) -> None:
    unpartitioned.md_parameters_topology = engines_structures.Amber_parameters_topology()
    unpartitioned.md_parameters_topology.prmtop_filename = prmtop_filename
    unpartitioned.positions = pdb_filename

def assign_forcefield_params(
        unpartitioned: engines_preparer.MD_unpartitioned_settings_input,
        built_in_ff_list: list[str],
        custom_ff_list: list[str],
        pdb_filename: str,
        ) -> None:
    unpartitioned.md_parameters_topology = engines_structures.Forcefield_parameters()
    unpartitioned.md_parameters_topology.built_in_forcefield_filenames = built_in_ff_list
    unpartitioned.md_parameters_topology.custom_forcefield_filenames = custom_ff_list
    unpartitioned.positions = pdb_filename

def assign_system_params(
        unpartitioned: engines_preparer.MD_unpartitioned_settings_input,
        system_filename: str,
        pdb_filename: str,
        ) -> None:
    unpartitioned.md_parameters_topology = engines_structures.Openmm_system()
    unpartitioned.md_parameters_topology.system_filename = system_filename
    unpartitioned.positions = pdb_filename

def assign_charmm_params(
        unpartitioned: engines_preparer.MD_unpartitioned_settings_input,
        psf_filename: str,
        charmm_ff_filenames: list[str],
        pdb_filename: str,
        ) -> None:
    unpartitioned.md_parameters_topology = engines_structures.Charmm_parameters_topology()
    unpartitioned.md_parameters_topology.psf_filename = psf_filename
    unpartitioned.md_parameters_topology.param_filename_list = charmm_ff_filenames
    unpartitioned.positions = pdb_filename


# ---------------------------------------------------------------------------
# Shared stage-building helpers (reduce duplication between creators)
# ---------------------------------------------------------------------------

def _make_state_reporter(interval: int, to_stdout: bool = True
                         ) -> reporters_preparer.State_reporter_input:
    r = reporters_preparer.State_reporter_input()
    r.to_stdout = to_stdout
    r.interval = interval
    return r

def _make_position_reporter(interval: int, atom_subset: str | None = None,
                             extension: str = "dcd",
                             enforce_periodic_box: bool = True,
                             ) -> reporters_preparer.Position_reporter_input:
    r = reporters_preparer.Position_reporter_input()
    r.interval = interval
    r.atom_subset = atom_subset
    r.extension = extension
    r.enforce_periodic_box = enforce_periodic_box
    return r

def _make_anchor_reporter(name: str, interval: int, extension: str = "dcd",
                          enforce_periodic_box: bool = True,
                          ) -> reporters_preparer.Anchor_reporter_input:
    r = reporters_preparer.Anchor_reporter_input()
    r.name = name
    r.interval = interval
    r.extension = extension
    r.enforce_periodic_box = enforce_periodic_box
    return r

# ---------------------------------------------------------------------------
# Public creator 1: toy system
# ---------------------------------------------------------------------------

def create_toy_input(
        root_dir: str,
        use_2cvs: bool = False,
        ) -> preparer.Seekr_input:
    """
    Create a toy (entropy-barrier) SEEKR input. Mirrors make_input.py's
    create_toy_model_input() signature.
    """
    if use_2cvs:
        values_list = [
            [0.1, -0.7], [-0.1, -0.5], [0.1, -0.3], [-0.1, -0.1],
            [0.1,  0.1], [-0.1,  0.3], [0.1,  0.5], [-0.1,  0.7],
        ]
        positions_list = [[[0.0, v[1], 0.0]] for v in values_list]
    else:
        values_list = [[-0.7], [-0.5], [-0.3], [-0.1], [0.1], [0.3], [0.5], [0.7]]
        positions_list = [[[0.0, v[0], 0.0]] for v in values_list]

    seekr_input = preparer.Seekr_input()
    seekr_input.directory = root_dir
    seekr_input.plugins = ["seekrsamplingmethodsplugin"]
    seekr_input.temperature = 298.15
    seekr_input.random_seed = 111

    md_settings = openmm_preparer.Openmm_settings_input()
    md_settings.platform_type = "reference"
    md_settings.unpartitioned_settings = engines_preparer.MD_unpartitioned_settings_input()
    sel = engines_structures.Selection()
    sel.name = "particle0"
    sel.atom_indices = [0]
    md_settings.unpartitioned_settings.selections = [sel]
    md_settings.unpartitioned_settings.positions = [[[0.0, -0.7, 0.0]]]
    seekr_input.engine_settings_list = [md_settings]

    seekr_input.initial_scopes = {
        "molecular_dynamics": "unpartitioned",
        "brownian_dynamics":  "unpartitioned",
        "logistic":           "unpartitioned",
    }

    cv_name = "distance_host_guest"

    state_point1 = base.State_point()
    state_point1.name = "stateA"
    state_point1.value = [0.0, -0.7] if use_2cvs else [-0.7]
    state_point1.end_state = True
    state_point2 = base.State_point()
    state_point2.name = "stateB"
    state_point2.value = [0.0, 0.7] if use_2cvs else [0.7]
    state_point2.end_state = True
    state_points = [state_point1, state_point2]

    anchor_scheme = anchor_scheme_preparer.Anchors_as_voronoi_cells_anchor_scheme()
    anchor_scheme.voronoi_cells = []
    for value, position in zip(values_list, positions_list):
        vc = anchor_scheme_preparer.Voronoi_cell()
        vc.value = value
        vc.positions = [position]
        vc.scale_types = ??? # scale type here. Example: ["molecular_dynamics"]
        anchor_scheme.voronoi_cells.append(vc)

    # --- stage 1: equilibration ---
    stage1 = md_preparer.Input_stage_md()
    stage1.name = "equilibration"
    stage1.run_minimization = False
    stage1.checkpoint_interval = 10
    stage1.ensemble = base.NVT_ensemble()
    stage1.sampling_method = sampling_preparer.Conventional_md_input()
    stage1.completion_criteria = compl_crit_struct.Number_of_steps_completion_criteria()
    stage1.completion_criteria.number_of_steps = 10000
    stage1.reporters = [
        _make_state_reporter(100),
        #_make_custom_reporter(100, enforce_periodic_box=False),
    ]
    # custom_restraint = ...
    #stage1.restraints = [custom_restraint]
    stage1.restraints = []

    # --- stage 2: metadynamics or steered MD ---
    ANCHOR_REPORTER = "my_anchor_reporter"
    stage2 = md_preparer.Input_stage_md()
    stage2.run_minimization = False
    stage2.checkpoint_interval = 10
    stage2.ensemble = base.NVT_ensemble()

    stage2.name = "metadynamics"
    stage2.sampling_method = sampling_preparer.Metadynamics_md_input()
    stage2.sampling_method.gaussian_height = 1.0
    stage2.sampling_method.cv_names = ["progress_coordinate"]

    stage2_anchor_reporter = _make_anchor_reporter(ANCHOR_REPORTER, 10, enforce_periodic_box=False)
    stage2.completion_criteria = compl_crit_struct.Reporter_progress_completion_criteria()
    stage2.completion_criteria.reporter_name = ANCHOR_REPORTER
    stage2.completion_criteria.interval = 10
    stage2.reporters = [stage2_anchor_reporter, _make_position_reporter(500, enforce_periodic_box=False)]

    stage3 = anchor_starting_structures_from_reporter_logistic_stage\
        .Input_stage_anchor_starting_structures_from_reporter()
    stage3.name = "anchor_structures_from_metaD"
    stage3.set_name = "main"
    stage3.reporter_name = ANCHOR_REPORTER

    # --- stage 4: MMVT ---
    stage4 = md_preparer.Input_stage_md()
    stage4.name = "MMVT"
    stage4.run_minimization = False
    stage4.checkpoint_interval = 5000
    stage4.ensemble = base.NVT_ensemble()
    stage4.sampling_method = sampling_preparer.MMVT_md_input()
    stage4.completion_criteria = compl_crit_struct.Number_of_steps_completion_criteria()
    stage4.completion_criteria.number_of_steps = 200000
    stage4.reporters = [
        _make_state_reporter(5000),
        _make_position_reporter(5000, enforce_periodic_box=False),
    ]

    seekr_input.stages = [stage1, stage2, stage3, stage4]

    # --- CVs ---
    external_cv_input = external_cv.External_cv_input()
    external_cv_input.name = "progress_coordinate"
    external_cv_input.expression = "y1"
    external_cv_input.group_names = ["particle0"]
    external_cv_input.min_value = -1.0
    external_cv_input.max_value = 1.0

    seekr_input.partitions = [preparer.Partition_input()]
    if use_2cvs:
        external_cv_input2 = external_cv.External_cv_input()
        external_cv_input2.name = "additional_coordinate"
        external_cv_input2.expression = "x1"
        external_cv_input2.group_names = ["particle0"]
        external_cv_input2.min_value = -1.0
        external_cv_input2.max_value = 1.0
        seekr_input.partitions[0].cvs = [external_cv_input2, external_cv_input]
    else:
        seekr_input.partitions[0].cvs = [external_cv_input]

    seekr_input.partitions[0].state_points = state_points
    seekr_input.partitions[0].anchor_scheme = anchor_scheme

    seekr_input.toy_settings = preparer.Toy_settings_input()
    seekr_input.toy_settings.potential_energy_expression = (
        "5*(x1^6+y1^6+exp(-(10*y1)^2)*(1-exp(-(10*x1)^2)))+0.5*z1^2"
    )
    seekr_input.toy_settings.num_particles = 1
    seekr_input.toy_settings.masses = [10.0]

    return seekr_input


# ---------------------------------------------------------------------------
# Public creator 2: host-guest molecular system
# ---------------------------------------------------------------------------

def create_host_guest_input(
        root_dir: str,
        ff: str = "amber",
        ) -> preparer.Seekr_input:
    """
    Create a host-guest (beta-cyclodextrin / 1-butanol) SEEKR input.
    All data file paths are absolute so the caller's cwd does not matter.
    """
    amber_prmtop_filename  = os.path.join(HOSTGUEST_DIRECTORY, "hostguest.parm7")
    #forcefield_built_in_ff = ["amber14/tip3pfb.xml"]
    #forcefield_custom_ff   = [os.path.join(HOSTGUEST_DIRECTORY, "hostguest.xml")]
    #system_filename        = os.path.join(HOSTGUEST_DIRECTORY, "hostguest_system.xml")

    pdb_basenames = [
        "hostguest_at0.5.pdb",  "hostguest_at1.5.pdb",  "hostguest_at2.5.pdb",
        "hostguest_at3.5.pdb",  "hostguest_at4.5.pdb",  "hostguest_at5.5.pdb",
        "hostguest_at6.5.pdb",  "hostguest_at7.5.pdb",  "hostguest_at8.5.pdb",
        "hostguest_at9.5.pdb",  "hostguest_at10.5.pdb", "hostguest_at11.5.pdb",
        "hostguest_at12.5.pdb", "",
    ]
    values_list = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75,
                   0.85, 0.95, 1.05, 1.15, 1.25, 1.35]

    receptor_indices_md  = list(range(147))
    ligand_indices_md    = list(range(147, 162))
    ligand_noh_indices   = [147, 150, 152, 155, 158]
    receptor_indices_bd  = list(range(147))
    ligand_indices_bd    = list(range(15))
    o5_group_indices     = [15, 36, 57, 78, 99, 120, 141]
    single_group_index   = [15]

    seekr_input = preparer.Seekr_input()
    seekr_input.directory = root_dir
    seekr_input.plugins = ["seekrsamplingmethodsplugin"]
    seekr_input.temperature = 298.15
    seekr_input.random_seed = 111

    md_settings = openmm_preparer.Openmm_settings_input()
    md_settings.unpartitioned_settings = engines_preparer.MD_unpartitioned_settings_input()

    if ff == "amber":
        assign_amber_params(
            md_settings.unpartitioned_settings,
            amber_prmtop_filename,
            os.path.join(HOSTGUEST_DIRECTORY, "hostguest_at0.5.pdb"),
        )
    #elif ff == "forcefield":
    #    assign_forcefield_params(
    #        md_settings.unpartitioned_settings,
    #        forcefield_built_in_ff, forcefield_custom_ff,
    #        os.path.join(HOSTGUEST_DIRECTORY, "hostguest_at0.5.pdb"),
    #    )
    #elif ff == "system":
    #    assign_system_params(
    #        md_settings.unpartitioned_settings,
    #        system_filename,
    #        os.path.join(HOSTGUEST_DIRECTORY, "hostguest_at0.5.pdb"),
    #    )
    else:
        raise ValueError(f"ff option not implemented: {ff}")

    def _sel(name, indices):
        s = engines_structures.Selection()
        s.name = name
        s.atom_indices = indices
        return s

    md_settings.unpartitioned_settings.selections = [
        _sel("host",        receptor_indices_md),
        _sel("guest",       ligand_indices_md),
        _sel("guest_noH",   ligand_noh_indices),
        _sel("solute",      receptor_indices_md + ligand_indices_md),
        _sel("o5_group",    o5_group_indices),
        _sel("single_index", single_group_index),
    ]
    seekr_input.engine_settings_list = [md_settings]

    seekr_input.initial_scopes = {
        "molecular_dynamics": "unpartitioned",
        "brownian_dynamics":  "unpartitioned",
        "logistic":           "unpartitioned",
    }

    # --- CV name depends on cv type ---
    cv_name = "distance_host_guest"

    # --- stage 1: equilibration ---
    from seekr.modules.restraints.preparer import Pairwise_restraint_input
    pairwise_restraint = Pairwise_restraint_input()
    pairwise_restraint.selection_name = "o5_group"
    stage1 = md_preparer.Input_stage_md()
    stage1.name = "equilibration"
    stage1.run_minimization = True
    stage1.checkpoint_interval = 1000
    stage1.ensemble = base.NVT_ensemble()
    stage1.sampling_method = sampling_preparer.Conventional_md_input()
    stage1.completion_criteria = compl_crit_struct.Number_of_steps_completion_criteria()
    stage1.completion_criteria.number_of_steps = 10000
    stage1.restraints = [pairwise_restraint]
    stage1.reporters = [
        _make_state_reporter(1000),
        _make_position_reporter(1000),
    ]

    # --- stage 2: metadynamics or steered MD ---
    ANCHOR_REPORTER = "my_anchor_reporter"
    stage2 = md_preparer.Input_stage_md()
    stage2.run_minimization = False
    stage2.checkpoint_interval = 1000
    stage2.ensemble = base.NVT_ensemble()

    stage2.name = "metadynamics"
    stage2.sampling_method = sampling_preparer.Metadynamics_md_input()
    stage2.sampling_method.gaussian_height = 1.5
    stage2.sampling_method.cv_names = [cv_name]

    from seekr.modules.restraints.preparer import (
        Proximity_restraint_input, Cylinder_restraint_input)
    proximity_restraint = Proximity_restraint_input()
    proximity_restraint.selection1_name = "host"
    proximity_restraint.selection2_name = "guest"
    proximity_restraint.radius = 1.05
    cylinder_restraint = Cylinder_restraint_input()
    cylinder_restraint.selection1_name = "guest"
    cylinder_restraint.selection2_name = "host"
    cylinder_restraint.selection3_name = "single_index"
    cylinder_restraint.radius = 0.4
    stage2.restraints = [proximity_restraint, cylinder_restraint]

    stage2_anchor_reporter = _make_anchor_reporter(ANCHOR_REPORTER, 100)
    stage2.completion_criteria = compl_crit_struct.Reporter_progress_completion_criteria()
    stage2.completion_criteria.reporter_name = ANCHOR_REPORTER
    stage2.completion_criteria.interval = 100
    stage2.reporters = [
        stage2_anchor_reporter,
        _make_state_reporter(500),
        _make_position_reporter(500, atom_subset="solute"),
    ]

    # --- stage 3: logistic ---
    stage3 = anchor_starting_structures_from_reporter_logistic_stage\
        .Input_stage_anchor_starting_structures_from_reporter()
    stage3.name = "anchor_structures_from_metaD"
    stage3.reporter_name = ANCHOR_REPORTER

    # --- stage 4: MMVT ---
    from seekr.modules.restraints.preparer import Positional_restraint_input
    positional_restraint = Positional_restraint_input()
    positional_restraint.selection_name = "host"
    stage4 = md_preparer.Input_stage_md()
    stage4.name = "MMVT"
    stage4.run_minimization = False
    stage4.checkpoint_interval = 1000
    stage4.ensemble = base.NVT_ensemble()
    stage4.sampling_method = sampling_preparer.MMVT_md_input()
    stage4.completion_criteria = compl_crit_struct.Number_of_steps_completion_criteria()
    stage4.completion_criteria.number_of_steps = 1000
    stage4.restraints = [positional_restraint]
    stage4.reporters = [
        _make_state_reporter(250),
        _make_position_reporter(250, atom_subset="solute"),
    ]

    seekr_input.stages = [stage1, stage2, stage3, stage4]

    # --- CVs ---
    seekr_input.partitions = [preparer.Partition_input()]
    dist_cv = distance_cv.Distance_cv_input()
    dist_cv.name = "distance_host_guest"
    dist_cv.group1_name = "host"
    dist_cv.group2_name = "guest"
    seekr_input.partitions[0].cvs = [dist_cv]

    bound = base.State_point()
    bound.name = "bound"; bound.value = [0.05]; bound.end_state = True
    bulk  = base.State_point()
    bulk.name  = "bulk";  bulk.value  = [1.35]; bulk.end_state  = True
    seekr_input.partitions[0].state_points = [bound, bulk]

    anchor_scheme = anchor_scheme_preparer.Anchors_as_voronoi_cells_anchor_scheme()
    anchor_scheme.voronoi_cells = []
    for value, pdb_basename in zip(values_list, pdb_basenames):
        vc = anchor_scheme_preparer.Voronoi_cell()
        vc.value = [value]
        vc.positions = ""
        vc.scale_types = ??? # Example: ["molecular_dynamics"]
        anchor_scheme.voronoi_cells.append(vc)
    seekr_input.partitions[0].anchor_scheme = anchor_scheme

    return seekr_input
