"""
Take the position reporter outputs from a previous stage 
representing a set of egress events (like from RAMD) and 
cluster the pathways in order to identify likely escape pathways.
"""

import os
import glob
import logging
import typing
import shutil

from attrs import define, field, validators, Factory
import numpy as np
from scipy.spatial import transform
from scipy.interpolate import splprep, splev
from sklearn.cluster import AffinityPropagation
import mdtraj
import networkx as nx

import seekr.modules.base as base
import seekr.modules.scales.logistic.base as logistic_scale_base

MIN_TRAJ_FRAMES = 20 #50
RMSD_SPACING = 0.15
# The mindist that both path waypoints must be less than in order to have their
# RMSDs clustered. The manuscript says, 0.3, but this seems very low, and
# also, they mention that the purpose of this is to make sure that at least one
# waypoint is not solvent exposed. They also say that they want to only 
# consider the the difference in pairs when contact with the protein is made.
# These three statements are definitely not consistent, but with the 
# intent in mind, we are choosing to only count pairs when both are within 
# 0.8 nm of the protein surface, as we are mainly interested to find the
# channels close to the bound state. However, if one want to find the routes
# to the channel, one could consider counting pairs when only ONE or the OTHER
# is within 0.8 of the protein surface.
RMSD_PATH_CLUSTER_MINDIST = 0.8

@define
class Stage_egress_pathway(logistic_scale_base.Stage_logistic):
    """
    Take position reporter outputs from a previous stage 
    representing a set of egress events (like from RAMD) and 
    cluster the pathways in order to identify likely escape pathways.
    """
    type: typing.Literal["egress_pathway_stage"] = "egress_pathway_stage"
    scale_type: typing.Literal["logistic"] = "logistic"
    index: int = field(default=0,
                        validator=[validators.instance_of(int),
                                      validators.ge(0)])
    name: str = field(default="egress_pathway", 
                      validator=validators.instance_of(str))
    comment: str = field(default="This is a logistic stage.",
                         validator=validators.instance_of(str))
    output_positions_basename_noext: str = field(
        default="centroid", 
        validator=validators.instance_of(str))
    output_positions_glob_list: list[str] = field(
        default=Factory(lambda: ["centroid*.dcd"]), 
        validator=validators.instance_of(list))
    input_stage_index: int = field(
        default=0,
        validator=[validators.instance_of(int),
                    validators.ge(0)])
    set_name: typing.Optional[str] = field(
        default="main", 
        validator=validators.optional(validators.instance_of(str)))
    reporter_name: str = field(
        default="",
        validator=validators.instance_of(str))
    ligand_selection_name: str = field(default="ligand", validator=validators.instance_of(str))
    receptor_selection_name: str = field(default="receptor", validator=validators.instance_of(str))
    scope: str = field(default="unpartitioned",
                       validator=validators.in_(["unpartitioned"]))
                           #["partitioned", "unpartitioned"]))

    def get_directory_name(self) -> str:
        """
        Get the directory of this stage.
        """
        directory_name: str = f"stage_{self.index}_{self.name}"
        return directory_name

    def add_stage_parameters_to_anchors(
            self, 
            model: "Seekr_model",
            ) -> None:
        """
        Add the stage directory to all anchors in the model.
        """
        for anchor in model.anchors:
            directory_name: str = self.get_directory_name()
            stage_parameter: dict[str, int | str] = {}
            stage_parameter["index"] = self.index
            stage_parameter["directory"] = directory_name
            stage_parameter["positions"] \
                = self.output_positions_basename_noext
            anchor.stage_parameters[self.name] = stage_parameter
        return

    def add_stage_parameters_to_unpartitioned_settings(
            self, 
            model: "Seekr_model"
            ) -> None:
        """
        Add the stage parameters to the unpartitioned settings.
        """
        directory_name: str = self.get_directory_name()
        stage_parameter: dict[str, int | str] = {}
        stage_parameter["index"] = self.index
        stage_parameter["directory"] = directory_name
        stage_parameter["positions"] \
            = self.output_positions_basename_noext
        
        md_engine = model.get_engine_by_scale_type("molecular_dynamics")
        assert md_engine.unpartitioned_settings is not None
        md_engine.unpartitioned_settings.stage_parameters[self.name] \
            = stage_parameter
        return

    def has_files(
            self, 
            model: "Seekr_model",
            anchor_index: int | None = None,
            ) -> bool:
        """
        Returns true if any of the stage directories contain
        outputs from any of the stages.
        """
        
        stage_output_globs: list[str] = self.output_positions_glob_list

        assert anchor_index is None, \
            "This stage must be run from an unpartitioned scope"
        stage_directory = self.get_directory_name()
        for stage_output_glob in stage_output_globs:
            output_files_glob = os.path.join(
                model.directory, base.UNPARTITIONED_DIRECTORY, 
                stage_directory, stage_output_glob)
            if len(glob.glob(output_files_glob)) > 0:
                return True
        return False

    def remove_directory(
            self, 
            model: "Seekr_model"
            ) -> None:
        """
        Remove entire directory in all anchors and unpartitioned 
        settings for this stage.
        """
        my_directory = self.get_directory_name()
        full_directory = os.path.join(
            model.directory, base.UNPARTITIONED_DIRECTORY, my_directory)
        if os.path.exists(full_directory):
            shutil.rmtree(full_directory)
        return
    
    def is_logistic_stage(self) -> bool:
        """
        Returns true since this is a logistic stage.
        """
        return True
    
    def load_sets_into_model(
            self, 
            model: "Seekr_model"
            ) -> None:
        """
        Load the CVs, selections, and anchor set for this stage 
        into the model.
        """
        model.load_cvs(self.set_name)
        model.load_selections(self.set_name)
        model.load_anchor_set(self.set_name)
        return
    
    def get_output_set_name(self) -> str | None:
        """
        Get the output set name for this stage.
        """
        return self.set_name
    
    def write_to_log(
            self,
            model: "Seekr_model",
            instructions: list[logistic_scale_base\
                               .Run_instruction_logistic]
        ) -> None:
        """
        Write the given instructions to every anchor's log file 
        within this stage's directory.
        """
        for instruction in instructions:
            log_filename = model.get_synopsis_filename()
            base.configure_logging(
                log_filename, overwrite=False)
            logging.info(instruction.message)

@define
class Input_stage_egress_pathway(logistic_scale_base.Input_stage_logistic):
    """
    This stage takes the output of a previous stage 
    representing a set of egress events (like from RAMD) and 
    clusters the pathways in order to identify likely escape pathways.
    """
    type: typing.Literal["egress_pathway_stage"] = "egress_pathway_stage"
    name: str = field(default="egress_pathway", 
                      validator=validators.instance_of(str))
    comment: str = field(default="Assign egress pathways from a previous stage.",
                         validator=validators.instance_of(str))
    input_stage_name: str = field(default="",
                                  validator=validators.instance_of(str))
    set_name: typing.Optional[str] = field(
        default="main", 
        validator=validators.optional(validators.instance_of(str)))
    reporter_name: str = field(
        default="",
        validator=validators.instance_of(str))
    ligand_selection_name: str = field(default="ligand", validator=validators.instance_of(str))
    receptor_selection_name: str = field(default="receptor", validator=validators.instance_of(str))
    
    def get_scale_type(self) -> str:
        return "logistic"

    def make_model_stage(
            self, 
            index: int,
            scope: str,
            dummy_settings: typing.Any = None,
            dummy_directory: str = ""
            ) -> typing.Tuple[Stage_egress_pathway, str]:
        """
        Given this input stage object, create the model's Stage object.
        """
        stage = Stage_egress_pathway()
        stage.index = index
        stage.name = self.name
        stage.comment = self.comment
        stage.set_name = self.set_name
        stage.reporter_name = self.reporter_name
        stage.scope = scope
        stage.ligand_selection_name = self.ligand_selection_name
        stage.receptor_selection_name = self.receptor_selection_name
        return stage, scope

def _get_topology_file(
        model: "Seekr_model",
        md_engine: typing.Any,
        input_stage: typing.Any,
        running_dir: str
        ) -> str:
    """
    Get the path to a topology/structure file for loading DCD 
    trajectories, mirroring the logic in 
    modules/engines/openmm/runner.py. For toy systems, returns 
    toy.pdb from the input stage's running directory. For real 
    systems, returns the input positions PDB used by that stage.
    """
    if model.using_toy():
        toy_pdb = os.path.join(running_dir, "toy.pdb")
        assert os.path.exists(toy_pdb), \
            f"toy.pdb not found in the input stage directory: {toy_pdb}"
        return toy_pdb
    else:
        assert input_stage.scope == "unpartitioned"
        # Get topology from the initial stage
        input_positions_filenames = md_engine.unpartitioned_settings\
            .get_positions_full_paths(
                #model.directory, input_stage.input_stage_index, 
                #input_stage.set_name, swarm_id=0)
                model.directory, 0, 
                "main", swarm_id=0)
        
        assert len(input_positions_filenames) > 0, \
            "No input positions/topology file found for DCD loading."
        return input_positions_filenames[0]

def _get_ligand_and_receptor_indices(
        model: "Seekr_model",
        stage: Stage_egress_pathway,
        ) -> tuple[list[int], list[int]]:
    """
    Get the indices of the ligand and receptor atoms.
    """
    ligand_indices = []
    receptor_indices = []
    md_engine = model.get_engine_by_scale_type("molecular_dynamics")
    assert md_engine.unpartitioned_settings is not None
    
    for selection in md_engine.unpartitioned_settings.selections:
        if selection.name == stage.ligand_selection_name:
            ligand_indices = selection.atom_indices
        elif selection.name == stage.receptor_selection_name:
            receptor_indices = selection.atom_indices
    assert len(ligand_indices) > 0, \
        f"No ligand indices found for selection {stage.ligand_selection_name}"
    assert len(receptor_indices) > 0, \
        f"No receptor indices found for selection {stage.receptor_selection_name}"
    return ligand_indices, receptor_indices

def _rmsd(
        ligand_indices: list[int],
        receptor_indices: list[int],
        positions: np.ndarray,
        ref_positions: np.ndarray,
        ) -> float:
    """
    Calculate the RMSD between ligand_indices of positions and ref_positions, 
    but aligned by receptor_indices.
    """
    
    pos_subset = []
    pos_subset_rmsd = []
    ref_subset = []
    ref_subset_rmsd = []
    for atom_index in receptor_indices:
        pos_subset.append(positions[atom_index])
        ref_subset.append(ref_positions[atom_index])
        
    for atom_index in ligand_indices:
        pos_subset_rmsd.append(positions[atom_index])
        ref_subset_rmsd.append(ref_positions[atom_index])

    pos_subset = np.array(pos_subset)
    ref_subset = np.array(ref_subset)
    pos_subset_rmsd = np.array(pos_subset_rmsd)
    ref_subset_rmsd = np.array(ref_subset_rmsd)
    
    avg_pos = np.array([float(np.mean(pos_subset[:,0])), 
                        float(np.mean(pos_subset[:,1])), 
                        float(np.mean(pos_subset[:,2]))])
    avg_ref = np.array([np.mean(ref_subset[:,0]), 
                        np.mean(ref_subset[:,1]), 
                        np.mean(ref_subset[:,2])])        
    pos_subset_centered = pos_subset - avg_pos
    ref_subset_centered = ref_subset - avg_ref
    pos_subset_rmsd_centered = pos_subset_rmsd - avg_pos
    ref_subset_rmsd_centered = ref_subset_rmsd - avg_ref
    rotation, rmsd_align = transform.Rotation.align_vectors(
        ref_subset_centered, pos_subset_centered)
    new_pos_subset_rmsd = rotation.apply(pos_subset_rmsd_centered)
    value = 0.0
    for i in range(len(ligand_indices)):
        increment = (new_pos_subset_rmsd[i,0] - ref_subset_rmsd_centered[i,0])**2 \
            + (new_pos_subset_rmsd[i,1] - ref_subset_rmsd_centered[i,1])**2 \
            + (new_pos_subset_rmsd[i,2] - ref_subset_rmsd_centered[i,2])**2
        value += increment

    rmsd = np.sqrt(value / len(ligand_indices))
    assert np.isfinite(rmsd), "Non-finite value detected."
    return rmsd

def _make_rmsd_matrix(
        traj: typing.Any,
        ligand_indices: list[int],
        receptor_indices: list[int],
        ) -> np.ndarray:
    """
    Make an RMSD matrix for the trajectory.
    """
    rmsd_matrix = np.zeros((len(traj), len(traj)))
    for i in range(len(traj)):
        for j in range(i, len(traj)):
            rmsd_matrix[i, j] = _rmsd(ligand_indices, receptor_indices, 
                                     traj.xyz[i], traj.xyz[j])
            rmsd_matrix[j, i] = rmsd_matrix[i, j]
    return rmsd_matrix

def _monte_carlo_equal_spacing(
        traj: typing.Any,
        rmsd_matrix: np.ndarray,
        ) -> typing.Any:
    """
    Perform a Monte Carlo simulation to find N equally spaced points
    in RMSD space along the trajectory.
    """
    N = 20
    num_monte_carlo_steps = 10000
    smoothed_indices = np.zeros(N, dtype=int)
    smoothed_indices[0] = 0
    smoothed_indices[-1] = len(traj)-1
    smoothed_indices[1:-1] = np.sort(np.random.choice(len(traj)-2, N-2, replace=False)+1)
    sum_rmsd = 0.0
    sum_rmsd_sq = 0.0
    for j in range(1,N):
        sum_rmsd += rmsd_matrix[smoothed_indices[j], smoothed_indices[j-1]]
        sum_rmsd_sq += rmsd_matrix[smoothed_indices[j], smoothed_indices[j-1]]**2

    avg_rmsd = sum_rmsd / (N-1)
    avg_rmsd_sq = sum_rmsd_sq / (N-1)
    variance = avg_rmsd_sq - (avg_rmsd**2)

    # Monte Carlo simulation to minimize variance between adjacent RMSD distances
    for i in range(num_monte_carlo_steps):
        new_smoothed_indices = np.copy(smoothed_indices)
        index = np.random.choice(N)
        if index == 0 or index == N-1:
            continue
        new_rough_traj_frame = np.random.choice(len(traj))
        if new_rough_traj_frame in new_smoothed_indices:
            # Don't allow duplicates
            continue
        new_smoothed_indices[index] = new_rough_traj_frame
        new_smoothed_indices = np.sort(new_smoothed_indices)
        
        sum_rmsd = 0.0
        sum_rmsd_sq = 0.0
        for j in range(1, N):
            sum_rmsd += rmsd_matrix[new_smoothed_indices[j], new_smoothed_indices[j-1]]
            sum_rmsd_sq += rmsd_matrix[new_smoothed_indices[j], new_smoothed_indices[j-1]]**2

        avg_rmsd = sum_rmsd / (N-1)
        avg_rmsd_sq = sum_rmsd_sq / (N-1)
        new_variance = avg_rmsd_sq - (avg_rmsd**2)
        
        if new_variance < variance:
            smoothed_indices = np.copy(new_smoothed_indices)
            variance = new_variance
            
    smoothed_trajectory = traj[smoothed_indices]
    return smoothed_trajectory

def _dijkstra_shortest_path(
        traj: typing.Any,
        rmsd_matrix: np.ndarray,
        ) -> typing.Any:
    """
    Perform a Dijkstra's algorithm to find the shortest path between the first frame
    (bound state) and the last frame (unbound state) in order to remove kinks/
    loops (kinetic traps).
    """
    G = nx.Graph()
    for i in range(len(traj)):
        G.add_node(i)
    for i in range(len(traj)):
        for j in range(i+1, len(traj)):
            G.add_edge(i, j, weight=rmsd_matrix[i, j]**2)
    shortest_path = nx.shortest_path(G, source=0, target=len(traj)-1, weight='weight')
    new_indices = np.array(shortest_path)
    new_traj = traj[new_indices]
    return new_traj

def _get_u2_values(
        traj: typing.Any,
        rmsd_matrix: np.ndarray,
        ) -> tuple[np.ndarray, np.ndarray]:
    """
    Calculate the number of points and total RMSD along a trajectory.
    """
    u = [0.0]
    rmsd_sum = 0.0
    for i in range(1, len(traj)):
        rmsd_sum += rmsd_matrix[i, i-1]
        u.append(rmsd_sum)
    u = np.array(u)
    num_points = int(rmsd_sum / RMSD_SPACING) + 1
    end_rmsd = (num_points-1) * RMSD_SPACING
    u2 = np.linspace(0, end_rmsd, num_points)
    assert np.isclose(u2[1] - u2[0], 0.15)
    return u, u2

def _interpolation(
        traj: typing.Any,
        ligand_indices: list[int],
        receptor_indices: list[int],
        u: np.ndarray,
        u2: np.ndarray,
        ) -> typing.Any:
    """
    Interpolate between frames in a trajectory to achieve 1.5 A spacing.
    """
    num_points = len(u2)
    new_positions = np.zeros((num_points, traj.xyz.shape[1], 3))
    for i in range(num_points):
        new_positions[i] = traj.xyz[0]

    for i, lig_index in enumerate(ligand_indices):
        old_points = traj.xyz[:, lig_index]
        tck, u1 = splprep(old_points.T, k=1, s=0, u=u)
        new_points = np.array(splev(u2, tck)).T
        new_positions[:, lig_index] = new_points

    for i, receptor_index in enumerate(receptor_indices):
        old_points = traj.xyz[:, receptor_index]
        tck, u1 = splprep(old_points.T, k=1, s=0, u=u)
        new_points = np.array(splev(u2, tck)).T
        new_positions[:, receptor_index] = new_points

    new_traj = mdtraj.Trajectory(new_positions, traj.topology)
    return new_traj

def _min_dist(
        ligand_indices: list[int],
        receptor_indices: list[int],
        positions: np.ndarray,
        ) -> float:
    """
    Calculate the minimum distance between the ligand and receptor atoms.
    """
    pos_subset_rec = []
    pos_subset_lig = []
    for atom_index in receptor_indices:
        pos_subset_rec.append(positions[atom_index])
        
    for atom_index in ligand_indices:
        pos_subset_lig.append(positions[atom_index])

    min_dist = 9e9
    for i in range(len(receptor_indices)):
        for j in range(len(ligand_indices)):
            dist = np.linalg.norm(pos_subset_lig[j] - pos_subset_rec[i])
            if dist < min_dist:
                min_dist = dist
        
    return min_dist
    
def _process_position_reporter(
        model: "Seekr_model",
        stage: Stage_egress_pathway,
        trajs: list,
        reporter: typing.Any,
        ) -> None:
    """
    For a position_reporter output, perform a Monte Carlo simulation to find
    N equally spaced points in RMSD space along the trajectory. Then perform
    a Dijkstra's algorithm to find the shortest path between the first frame
    (bound state) and the last frame (unbound state) in order to remove kinks/
    loops (kinetic traps). The rough trajectory is then smoothed to achieve
    as close to 1.5 A spacing as possible between adjacent frames.
    """
    stage_dir_name = stage.get_directory_name()
    if stage.scope == "unpartitioned":
        full_stage_dir = os.path.join(
            model.directory, base.UNPARTITIONED_DIRECTORY, 
            stage_dir_name)
    else:
        raise NotImplementedError(
            "egress_pathway_stage with partitioned scope is not yet implemented.")
    
    ligand_indices, receptor_indices = _get_ligand_and_receptor_indices(model, stage)
    counter = 0
    max_num_points = 0
    max_u2_value = 0.0
    u_list = []
    dijkstra_trajs = []
    included_indices = []
    for i, traj in enumerate(trajs):
        if traj.n_frames < MIN_TRAJ_FRAMES:
            continue # Skip trajectories that are too short.

        rmsd_matrix = _make_rmsd_matrix(traj, ligand_indices, receptor_indices)
        smooth_traj = _monte_carlo_equal_spacing(traj, rmsd_matrix)
        smooth_traj_name = os.path.join(full_stage_dir, 
                                        f"smoothed_{counter}.dcd")
        logging.info(f"Saving smoothed trajectory to {smooth_traj_name}")
        smooth_traj.save(smooth_traj_name)
        smooth_rmsd_matrix = _make_rmsd_matrix(
            smooth_traj, ligand_indices, receptor_indices)
        dijkstra_traj = _dijkstra_shortest_path(smooth_traj, smooth_rmsd_matrix)
        dijkstra_traj_name = os.path.join(full_stage_dir, 
                                          f"dijkstra_{i}.dcd")
        dijkstra_trajs.append(dijkstra_traj)
        logging.info(f"Saving Dijkstra trajectory to {dijkstra_traj_name}")
        dijkstra_traj.save(dijkstra_traj_name)
        dijkstra_rmsd_matrix = _make_rmsd_matrix(
            dijkstra_traj, ligand_indices, receptor_indices)
        u, u2 = _get_u2_values(dijkstra_traj, dijkstra_rmsd_matrix)
        u_list.append(u)
        if len(u2) > max_num_points:
            max_num_points = len(u2)
            max_u2_value = u2[-1]
        included_indices.append(i)
        counter += 1

    # Make extrapolated u2 values
    new_u2 = np.linspace(0.0, max_u2_value, max_num_points)
    summary_traj_list = []
    for i, index in enumerate(included_indices):
        dijkstra_traj = dijkstra_trajs[i]
        summary_traj = _interpolation(
            dijkstra_traj, ligand_indices, receptor_indices, u_list[i], new_u2)
        summary_traj_name = os.path.join(full_stage_dir, 
                                         f"summary_{index}.dcd")
        summary_traj_list.append(summary_traj)
        logging.info(f"Saving summary trajectory to {summary_traj_name}")
        # This code cut out any of the pathways exiting much beyond the surface
        #summary_traj[0:len(u_list[i])].save(summary_traj_name)
        # This code keeps the extended pathways flying out into the unbound state
        summary_traj.save(summary_traj_name)

    path_distance_matrix = np.zeros((len(included_indices), len(included_indices)))
    for i, index1 in enumerate(included_indices):
        summary_traj_i = summary_traj_list[i]
        for j, index2 in enumerate(included_indices):
            if i == j:
                continue
            summary_traj_j = summary_traj_list[j]
            for k in range(max_num_points):
                # See whether the mindist of both is greater than 0.3 nm.
                #  If so, skip. (Changing to 0.5 nm)
                min_dist_i = _min_dist(ligand_indices, receptor_indices, 
                                        summary_traj_i.xyz[k])
                min_dist_j = _min_dist(ligand_indices, receptor_indices,
                                        summary_traj_j.xyz[k])
                if min_dist_i > RMSD_PATH_CLUSTER_MINDIST and \
                        min_dist_j > RMSD_PATH_CLUSTER_MINDIST:
                    continue
                path_distance_matrix[i, j] += _rmsd(
                    ligand_indices, receptor_indices, summary_traj_i.xyz[k], 
                    summary_traj_j.xyz[k])

    # AffinityPropagation expects similarities (higher = more similar). Match
    # sklearn's Euclidean default: negative squared distance.
    path_affinity_matrix = -path_distance_matrix ** 2
    clustering = AffinityPropagation(
        affinity="precomputed",
        random_state=model.get_random_seed(),
    ).fit(path_affinity_matrix)
    cluster_sizes = np.bincount(clustering.labels_)
    cluster_sizes_sorted = np.argsort(cluster_sizes)[::-1]
    cluster_center_indices_sorted = []
    for i in cluster_sizes_sorted:
        cluster_center_indices_sorted.append(clustering.cluster_centers_indices_[i])
    for i, cluster_index in enumerate(cluster_center_indices_sorted):
        file_index = included_indices[cluster_index]
        centroid_filename = os.path.join(full_stage_dir, f"centroid_{i}_{file_index}.dcd")
        traj_filename = os.path.join(full_stage_dir, f"summary_{file_index}.dcd")
        shutil.copyfile(traj_filename, centroid_filename)

    logging.info(
        f"Assigned egress pathways {full_stage_dir}.")
    return

def is_finished_egress_pathway_logistic_with_message(
        model: "Seekr_model",
        stage: Stage_egress_pathway,
        use_json: bool = True
    ) -> tuple[bool, str | dict]:
    """
    A egress pathway logistic stage is finished when all of the new egress pathway
    directories have files in them.
    """
    all_have_files = True
    if use_json:
        all_have_files = stage.has_files(model, anchor_index=None)
        final_message = {
            "anchor_index": None,
            "has_files": all_have_files
        }
        return all_have_files, final_message
    else:
        message = ""
        if stage.has_files(model, anchor_index=None):
            message = "Structures generated."
        else:
            all_have_files = False
            message = "No structures yet generated."
        return all_have_files, message

def choose_next_simulation_egress_pathway_logistic(
        model: "Seekr_model",
        stage: Stage_egress_pathway,
        force_overwrite: bool,
        use_json: bool = True,
        silent: bool = False
        ) -> list:
    """
    Choose the next simulation to run for a swarm structures from reporter logistic stage.
    """
    run_instructions: list[typing.NamedTuple] = []
    completed, message \
        = is_finished_egress_pathway_logistic_with_message(
            model, stage, use_json)
    if not completed or force_overwrite:
        run_instruction = logistic_scale_base.Run_instruction_logistic(
            anchor_index=None,
            swarm_id=0,
            message=message
        )
        run_instructions.append(run_instruction)
    return run_instructions

def do_run_instruction_egress_pathway_logistic(
        model: "Seekr_model",
        force_overwrite: bool,
        stage_index: int,
        instruction: list
        ) -> None:
    """
    Read egress pathway output files from the predecessor MD stage,
    extract one representative frame per anchor, and write 
    starting structures for the next stage.
    """
    stage: Stage_egress_pathway \
        = model.stages[stage_index]
    input_stage = model.stages[stage.input_stage_index-1]
    md_engine = model.get_engine_by_scale_type("molecular_dynamics")
    assert md_engine is not None, "No MD engine found in model."
    assert stage.scope == "unpartitioned", \
        "This stage must be run from an unpartitioned scope"
    assert md_engine.unpartitioned_settings is not None, \
        "No unpartitioned settings found in MD engine."
    unpartitioned_settings = md_engine.unpartitioned_settings
    running_dir = os.path.join(
        model.directory, unpartitioned_settings.directory,
        input_stage.get_directory_name())
    basename = input_stage.sampling_method\
        .get_output_positions_basename_noext()
    assert stage.reporter_name != "", \
        "reporter_name must be set on the logistic stage."
    reporter = None
    for r in input_stage.reporters:
        if r.name == stage.reporter_name:
            reporter = r
            break
    assert reporter is not None, \
        f"Reporter '{stage.reporter_name}' not found in input stage "\
        f"'{input_stage.name}'. Check that reporter_name matches the "\
        "reporter's name."
    ext = reporter.get_reporter_extension()
    reporter_glob = os.path.join(
        running_dir, f"{basename}*.{reporter.name}.{ext}")
    reporter_files = sorted(glob.glob(reporter_glob))
    assert len(reporter_files) > 0, \
        f"No reporter output files found matching: {reporter_glob}"
    topology_file = _get_topology_file(model, md_engine, input_stage, running_dir)
    # Build (traj, json_data) pairs for every restart file that has a JSON sidecar.
    trajs = []
    # Concatenate the starting frame to the trajs!
    initial_stage_dir = os.path.join(
        model.directory, base.UNPARTITIONED_DIRECTORY, "stage_0_initial")
    initial_struct_basename = unpartitioned_settings.stage_parameters["initial"]["positions"]
    initial_struct_file = os.path.join(initial_stage_dir, f"{initial_struct_basename}.swarm_0.pdb")
    initial_traj = mdtraj.load(initial_struct_file, top=topology_file)
    for reporter_file in reporter_files:
        if ext == "pdb":
            traj = mdtraj.load(reporter_file)
        elif ext == "dcd":
            traj = mdtraj.load(reporter_file, top=topology_file)
        else:
            raise Exception(f"Unsupported reporter file extension: {ext}")
        combined_traj = mdtraj.join([initial_traj, traj])
        trajs.append(combined_traj)
    assert len(trajs) > 0, \
        f"No reporter files found matching: {reporter_glob}"
    _process_position_reporter(model, stage, trajs, reporter)
    return
