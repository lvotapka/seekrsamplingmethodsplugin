from __future__ import annotations
"""
modules/analyzers/random_accelerated_analyzer.py

Analyze the results of a Random Accelerated Molecular Dynamics (RAMD)
calculation.
"""

import os
import glob
import random
import statistics

from attrs import define, field, validators

import seekr.modules.base as base
import seekr.modules.structures as structures
import seekr.modules.analyzers.base as analyzers_base

DEFAULT_NUM_ERROR_SAMPLES = 200
BOOTSTRAP_FRACTION = 0.8

@define
class Ramd_analyzer:
    """
    Analyze the results of a Random Accelerated Molecular Dynamics 
    (RAMD) calculation.
    """
    model: structures.Seekr_model = field(
        validator=validators.instance_of(structures.Seekr_model))
    stage_index: int = field(
        validator=validators.instance_of(int),)
    lengths_steps: list[float] = field(
        factory=list, validator=validators.instance_of(list))
    tau_mean: float | None = field(
        default=None, validator=validators.optional(
        validators.instance_of(float)))
    tau_median: float | None = field(
        default=None, validator=validators.optional(
        validators.instance_of(float)))
    tau_err_median: float | None = field(
        default=None, validator=validators.optional(
        validators.instance_of(float)))
    tau_err_stdev: float | None = field(
        default=None, validator=validators.optional(
        validators.instance_of(float)))
    
    def extract_data(self) -> None:
        """
        Extract the necessary data from the RAMD calculation 
        output files.
        """
        stage = self.model.stages[self.stage_index-1]
        assert stage.scope == "unpartitioned"
        output_file_glob = os.path.join(
            self.model.directory, base.UNPARTITIONED_DIRECTORY,
            stage.get_directory_name(), "*.out")
        output_file_list = glob.glob(output_file_glob)
        output_file_list = base.order_files_numerically(
            output_file_list)
        step_counts = []
        for output_file in output_file_list:
            with open(output_file, "r") as f:
                for line in f:
                    if line.startswith("step_count:"):
                        length_steps = int(line.split()[-1])
                        step_counts.append(length_steps)
        self.lengths_steps = step_counts

    def process_data(self) -> None:
        """
        Process the extracted data to calculate any necessary 
        quantities for analysis.
        """
        if self.lengths_steps:
            steps_mean = statistics.mean(self.lengths_steps)
            self.tau_mean = steps_mean \
                * self.model.get_timestep_by_type("molecular_dynamics")
            steps_median = statistics.median(self.lengths_steps)
            self.tau_median = steps_median \
                * self.model.get_timestep_by_type("molecular_dynamics")
        
        return

    def error_analysis(
        self,
        error_analysis_method: str,
        num_error_samples: int | None = None,
        random_seed: int | None = None
    ) -> None:
        """
        Perform error analysis on the calculated quantities.
        """
        if len(self.lengths_steps) < 2:
            # Not enough data for error analysis
            self.tau_err_median = None
            self.tau_err_stdev = None
            return
        assert error_analysis_method in ["bootstrap"]
        if error_analysis_method == "bootstrap":
            fraction_count = int(len(self.lengths_steps) \
                                 * BOOTSTRAP_FRACTION)
            bootstrap_steps = []
            for i in range(num_error_samples or DEFAULT_NUM_ERROR_SAMPLES):
                if random_seed is not None:
                    random_seed_i = random_seed + i
                else:
                    random_seed_i = None
                random.seed(random_seed_i)
                #bootstrap_sample = random.sample(
                #    self.lengths_steps, fraction_count)
                bootstrap_sample = random.choices(
                    self.lengths_steps, k=len(self.lengths_steps))
                bootstrap_median = statistics.median(bootstrap_sample)
                bootstrap_steps.append(bootstrap_median)

            self.tau_err_median = statistics.median(bootstrap_steps) \
                * self.model.get_timestep_by_type("molecular_dynamics")
            self.tau_err_stdev = statistics.stdev(bootstrap_steps) \
                * self.model.get_timestep_by_type("molecular_dynamics")
        return
    
    def get_statistics_message(
            self,
            model: structures.Seekr_model,
            anchor_index_instruction: str | None = None,
            use_json: bool = False
            ) -> str | dict:
        """

        """
        if use_json:
            final_message = {"escape_events_counted": 
                             len(self.lengths_steps)}
        else:
            final_message = f"  escape events counted: "\
                f"{len(self.lengths_steps)}\n  "

        return final_message
    
    def print_results(self):
        """Print all results of the analysis calculation."""
        print("Printing results from RAMD SEEKR calculation")
        if self.tau_median is not None:
            print("Tau (ps):", 
                  analyzers_base.pretty_string_value_error(
                    self.tau_median, self.tau_err_stdev))
        
        return
    
    def make_image_directory(
            self, 
            image_directory: str
            ) -> None:
        """
        Create the directory and sub-directories for images to be 
        written. Does nothing at present...
        """
        #if image_directory != "" and not os.path.exists(image_directory):
        #    os.mkdir(image_directory)
        return
    
    def save_plots(
            self, 
            image_directory: str
            ) -> None:
        """
        Save a potentially useful series of plots of some quantities
        obtained during the analysis.
        
        Does nothing at present.
        """
        return
        

def analyze_random_accelerated(
        model: structures.Seekr_model,
        stage_index: int,
        skip_checks: bool = False,
        error_analysis_method: str = "bootstrap",
        only_extract_data: bool = False,
        num_error_samples: int | None = None
        ) -> Ramd_analyzer:
    """
    Analyze the results of a Random Accelerated Molecular Dynamics 
    (RAMD) calculation.
    
    Parameters
    ----------
    model : structures.Seekr_model
        The SEEKR model to analyze.
    
    Returns
    -------
    Ramd_analyzer
        An instance of Ramd_analyzer containing the analysis results.
    """
    analysis = Ramd_analyzer(model=model, stage_index=stage_index)
    analysis.extract_data()
    # TODO: implement some sort of checks here...
    #if not skip_checks:
    #    check_extraction(...)
    if only_extract_data:
        return analysis
    analysis.process_data()
    analysis.error_analysis(error_analysis_method, num_error_samples,
                            random_seed=model.get_random_seed())
    return analysis