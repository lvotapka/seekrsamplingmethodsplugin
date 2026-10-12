from __future__ import annotations
"""
This logistic stage will take the outputs to a reporter and 
will write them as swarm starting structures for any stage
with this as a dependency.
"""

import os
import glob
import json
import logging
import typing
from shutil import rmtree

import numpy as np
import mdtraj
from attrs import define, field, validators, Factory

import seekr.modules.base as base
import seekr.modules.scales.logistic.base as logistic_scale_base

OUTPUT_POSITIONS_BASENAME_NOEXT = "position"

@define
class Stage_swarm_structures_from_reporter(
        logistic_scale_base.Stage_logistic):
    """
    This stage takes the output of a reporter within an input stage 
    and writes them as swarm starting structures for any stage
    with this as a dependency.

    If there's more than one reporter in the stage, assume the first one.
    """
    type: typing.Literal["swarm_structures_from_reporter_stage"] \
        = "swarm_structures_from_reporter_stage"
    scale_type: typing.Literal["logistic"] = "logistic"
    index: int = field(default=0,
                        validator=[validators.instance_of(int),
                                      validators.ge(0)])
    name: str = field(default="swarm_structures_from_reporter", 
                      validator=validators.instance_of(str))
    comment: str = field(default="This is a logistic stage.",
                         validator=validators.instance_of(str))
    output_positions_basename_noext: str = field(
        default="output_positions", 
        validator=validators.instance_of(str))
    output_positions_glob_list: list[str] = field(
        default=Factory(lambda: ["output_positions*.pdb", 
                                 "output_positions*.txt"]), 
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
    number_of_velocity_resamples: int = field(
        default=1,
        validator=[validators.instance_of(int),
                    validators.ge(1)])
    scope: str = field(default="partitioned",
                       validator=validators.in_(
                           ["partitioned", "unpartitioned"]))
    
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
        my_directory = self.get_directory_name()
        anchor_set_directory = os.path.join(
            model.directory, "anchor_sets", self.get_output_set_name())

        if anchor_index is None:
            stage_directory = self.get_directory_name()
            for stage_output_glob in stage_output_globs:
                output_files_glob = os.path.join(
                    model.directory, base.UNPARTITIONED_DIRECTORY, 
                    stage_directory, stage_output_glob)
                if len(glob.glob(output_files_glob)) > 0:
                    return True
        else: 
            for anchor in model.anchors:
                if anchor_index is not None \
                        and anchor.index != anchor_index:
                    continue
                for stage_parameter in anchor.stage_parameters.values():
                    stage_directory = stage_parameter["directory"]
                    if stage_directory == my_directory:
                        for stage_output_glob in stage_output_globs:
                            output_files_glob = os.path.join(
                                anchor_set_directory, anchor.directory, 
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
        anchor_set_directory = os.path.join(
            model.directory, "anchor_sets", self.get_output_set_name())
        for anchor in model.anchors:
            full_directory = os.path.join(
                anchor_set_directory, anchor.directory, my_directory)
            if os.path.exists(full_directory):
                rmtree(full_directory)
        md_engine = model.get_engine_by_scale_type("molecular_dynamics")
        assert md_engine.unpartitioned_settings is not None
        full_directory = os.path.join(
            model.directory, md_engine.unpartitioned_settings.directory, my_directory)
        if os.path.exists(full_directory):
            rmtree(full_directory)
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
class Input_stage_swarm_structures_from_reporter(
        logistic_scale_base.Input_stage_logistic):
    """
    This stage assigns the starting structures of the anchors of the model
    according to the output of a reporter stage.
    """
    type: typing.Literal["swarm_structures_from_reporter_stage"] \
        = "swarm_structures_from_reporter_stage"
    name: str = field(default="swarm_structures_from_reporter", 
                      validator=validators.instance_of(str))
    comment: str = field(
        default="Assign swarm starting structures from a reporter.",
        validator=validators.instance_of(str))
    # If these two are None, then assume previous/next stage.
    input_stage_name: str | None = field(default=None,
                                 validator=validators.optional(
                                     validators.instance_of(str)))
    set_name: typing.Optional[str] = field(
        default="main", 
        validator=validators.optional(validators.instance_of(str)))
    reporter_name: str = field(
        default="",
        validator=validators.instance_of(str))
    number_of_velocity_resamples: int = field(
        default=1,
        validator=[validators.instance_of(int),
                    validators.ge(1)])

    def get_scale_type(self) -> str:
        return "logistic"

    def make_model_stage(
            self, 
            index: int,
            scope: str,
            dummy_settings: typing.Any = None,
            dummy_directory: str = ""
            ) -> typing.Tuple[Stage_swarm_structures_from_reporter, str]:
        """
        Given this input stage object, create the model's Stage object.
        """
        stage = Stage_swarm_structures_from_reporter()
        stage.index = index
        stage.name = self.name
        stage.comment = self.comment
        stage.set_name = self.set_name
        stage.reporter_name = self.reporter_name
        stage.number_of_velocity_resamples \
            = self.number_of_velocity_resamples
        stage.scope = scope
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
        if input_stage.scope == "unpartitioned":
            # Get topology from the initial stage
            input_positions_filenames = md_engine.unpartitioned_settings\
                .get_positions_full_paths(
                    #model.directory, input_stage.input_stage_index, 
                    #input_stage.set_name, swarm_id=0)
                    model.directory, 0, 
                    "main", swarm_id=0)
        else:
            for anchor in model.anchors:
                if anchor.bulk_anchor:
                    continue
                stage_dir_name = input_stage.get_directory_name()
                anchor_set_directory = os.path.join(
                    model.directory, "anchor_sets", 
                    input_stage.get_output_set_name())
                full_stage_dir = os.path.join(
                    anchor_set_directory, stage_dir_name)
                basename_glob = input_stage.sampling_method\
                    .get_output_positions_basename_noext()+"*.pdb"
                input_positions_glob = os.path.join(
                    full_stage_dir, basename_glob)
                input_positions_filenames = glob.glob(
                    input_positions_glob)
        
        assert len(input_positions_filenames) > 0, \
            "No input positions/topology file found for DCD loading."
        return input_positions_filenames[0]

def _write_frame(
        traj: typing.Any,
        frame_num: int,
        output_path: str,
        using_toy: bool
        ) -> None:
    """
    Write a single frame from an mdtraj Trajectory to disk.
    Creates the parent directory if it does not already exist.
    Toy systems produce .txt (NumPy savetxt); real systems 
    produce .pdb.
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    if using_toy:
        np.savetxt(output_path, traj.xyz[frame_num])
    else:
        traj[frame_num].save_pdb(output_path)
    return

def _process_position_reporter(
        model: "Seekr_model",
        stage: Stage_swarm_structures_from_reporter,
        trajs: list,
        reporter: typing.Any,
        ) -> None:
    """
    For a position_reporter output, extract the frames and generate
    starting structures for the next stage.
    """
    output_ext = ".txt" if model.using_toy() else ".pdb"
    stage_dir_name = stage.get_directory_name()
    if stage.scope == "unpartitioned":
        full_stage_dir = os.path.join(
            model.directory, base.UNPARTITIONED_DIRECTORY, 
            stage_dir_name)
    else:
        anchor_set_directory = os.path.join(
                model.directory, "anchor_sets", 
                stage.get_output_set_name())
        full_stage_dir = os.path.join(
            anchor_set_directory, stage_dir_name)
    counter = 0
    for traj in trajs:
        for frame_num in range(traj.n_frames):
            if frame_num >= len(traj):
                continue
            for _ in range(stage.number_of_velocity_resamples):
                output_path = os.path.join(
                    full_stage_dir,
                    stage.output_positions_basename_noext + f".swarm_{counter}" \
                        + output_ext)
                _write_frame(traj, frame_num, output_path, 
                            model.using_toy())
                counter += 1

    logging.info(
        f"Assigned swarm structures {full_stage_dir}.")
    return

def is_finished_swarm_structures_from_reporter_logistic_with_message(
        model: "Seekr_model",
        stage: Stage_swarm_structures_from_reporter,
        use_json: bool = True
    ) -> tuple[bool, str | dict]:
    """
    A reassign anchors logistic stage is finished when all of the new anchor
    directories have files in them.
    """
    all_have_files = True
    if use_json:
        final_message: dict[int, dict] = {}
        if stage.scope == "partitioned":  
            for anchor in model.anchors:
                has_files = stage.has_files(model, anchor_index=anchor.index)
                progress = 0.0 if has_files else 1.0
                if not has_files:
                    all_have_files = False
                final_message[anchor.index] = {
                    "finished": has_files,
                    "progress": progress,
                    "num_swarms": 1,
                    "swarms": {"0":{
                        "finished": has_files,
                        "progress": progress,
                    }}
                }
            return all_have_files, final_message
        else:
            all_have_files = stage.has_files(model, anchor_index=None)
            progress = 0.0 if all_have_files else 1.0
            final_message = {
                "finished": all_have_files,
                "progress": progress,
                "num_swarms": 1,
                "swarms": {"0":{
                    "finished": all_have_files,
                    "progress": progress,
                }}
            }
            return all_have_files, final_message
    else:
        message = ""
        if stage.scope == "partitioned":
            for anchor in model.anchors:
                if not stage.has_files(model, anchor_index=anchor.index):
                    all_have_files = False
                    message += f"    Anchor {anchor.index} has no structures.\n"
            if all_have_files:
                message = "    All anchors assigned structures."
        else:
            if stage.has_files(model, anchor_index=None):
                message = "Structures generated."
            else:
                all_have_files = False
                message = "No structures yet generated."
        return all_have_files, message

def choose_next_simulation_swarm_structures_from_reporter_logistic(
        model: "Seekr_model",
        stage: Stage_swarm_structures_from_reporter,
        force_overwrite: bool,
        use_json: bool = True,
        silent: bool = False
        ) -> list:
    """
    Choose the next simulation to run for a swarm structures from reporter logistic stage.
    """
    run_instructions: list[typing.NamedTuple] = []
    completed, message \
        = is_finished_swarm_structures_from_reporter_logistic_with_message(
            model, stage, use_json)
    if not completed or force_overwrite:
        run_instruction = logistic_scale_base.Run_instruction_logistic(
            anchor_index=None,
            swarm_id=0,
            message=message
        )
        run_instructions.append(run_instruction)
    return run_instructions

def do_run_instruction_swarm_structures_from_reporter_logistic(
        model: "Seekr_model",
        force_overwrite: bool,
        stage_index: int,
        instruction: list
        ) -> None:
    """
    Read reporter output files from the predecessor MD stage,
    extract one representative frame per anchor, and write 
    starting structures for the next stage.
    """
    stage: Stage_swarm_structures_from_reporter \
        = model.stages[stage_index]
    input_stage = model.stages[stage.input_stage_index-1]
    md_engine = model.get_engine_by_scale_type("molecular_dynamics")
    assert md_engine is not None, "No MD engine found in model."
    if stage.scope == "unpartitioned":
        assert md_engine.unpartitioned_settings is not None, \
            "No unpartitioned settings found in MD engine."
        unpartitioned_settings = md_engine.unpartitioned_settings
        running_dir = os.path.join(
            model.directory, unpartitioned_settings.directory,
            input_stage.get_directory_name())
    else:
        raise NotImplementedError(
            "swarm_structures_from_reporter_stage with " \
            "partitioned scope is not yet implemented.")
    basename = input_stage.sampling_method\
        .get_output_positions_basename_noext()
    swarm_string = ".swarm_0"
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
        running_dir, f"{basename}*{swarm_string}.{reporter.name}.{ext}")
    reporter_files = sorted(glob.glob(reporter_glob))
    assert len(reporter_files) > 0, \
        f"No reporter output files found matching: {reporter_glob}"
    topology_file = _get_topology_file(model, md_engine, input_stage, running_dir)
    # Build (traj, json_data) pairs for every restart file that has a JSON sidecar.
    trajs = []
    for reporter_file in reporter_files:
        if os.path.getsize(reporter_file) == 0:
            continue
        try:
            if ext == "pdb":
                traj = mdtraj.load(reporter_file)
            elif ext == "dcd":
                traj = mdtraj.load(reporter_file, top=topology_file)
            else:
                raise Exception(f"Unsupported reporter file extension: {ext}")
            
        except OSError:
           continue
        trajs.append(traj)
    assert len(trajs) > 0, \
        f"No reporter files found matching: {reporter_glob}"
    _process_position_reporter(model, stage, trajs, reporter)
    return
