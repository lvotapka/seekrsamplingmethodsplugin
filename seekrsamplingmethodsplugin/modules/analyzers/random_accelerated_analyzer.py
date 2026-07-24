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
import numpy as np

import seekr.modules.base as base
import seekr.modules.structures as structures
import seekr.modules.analyzers.base as analyzers_base

DEFAULT_NUM_ERROR_SAMPLES = 200
BOOTSTRAP_FRACTION = 0.8
DEFAULT_NUM_WINDOWS = 5

def bootstrap_ramd(
        lengths_steps: list[float],
        num_error_samples: int | None, 
        random_seed: int | None, 
        timestep: float,
        ) -> tuple[float | None, float | None]:
    bootstrap_mean_steps = []
    bootstrap_median_steps = []
    for i in range(num_error_samples or DEFAULT_NUM_ERROR_SAMPLES):
        if random_seed is not None:
            random_seed_i = random_seed + i
        else:
            random_seed_i = None
        random.seed(random_seed_i)
        bootstrap_sample = random.choices(
            lengths_steps, k=len(lengths_steps))
        bootstrap_mean = statistics.mean(bootstrap_sample)
        bootstrap_mean_steps.append(bootstrap_mean)
        bootstrap_median = statistics.median(bootstrap_sample)
        bootstrap_median_steps.append(bootstrap_median)

    tau_err_mean = statistics.mean(bootstrap_mean_steps) \
        * timestep
    tau_err_mean_stdev = statistics.stdev(bootstrap_mean_steps) \
        * timestep
    tau_err_median = statistics.median(bootstrap_median_steps) \
        * timestep
    tau_err_median_stdev = statistics.stdev(bootstrap_median_steps) \
        * timestep
    
    return tau_err_mean, tau_err_mean_stdev, tau_err_median, tau_err_median_stdev
    
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
    num_windows: int = field(
        default=DEFAULT_NUM_WINDOWS, validator=validators.instance_of(int))
    window_lengths_steps: list[float] = field(
        factory=list, validator=validators.instance_of(list))
    tau_mean: float | None = field(
        default=None, validator=validators.optional(
        validators.instance_of(float)))
    tau_mean_windows: list[float] = field(
        factory=list, validator=validators.instance_of(list))
    tau_median: float | None = field(
        default=None, validator=validators.optional(
        validators.instance_of(float)))
    tau_median_windows: list[float] = field(
        factory=list, validator=validators.instance_of(list))
    tau_err_mean: float | None = field(
        default=None, validator=validators.optional(
        validators.instance_of(float)))
    tau_err_mean_stdev: float | None = field(
        default=None, validator=validators.optional(
        validators.instance_of(float)))
    tau_err_mean_stdev_windows: list[float] = field(
        factory=list, validator=validators.instance_of(list))
    tau_err_mean_windows: list[float] = field(
        factory=list, validator=validators.instance_of(list))
    tau_err_median: float | None = field(
        default=None, validator=validators.optional(
        validators.instance_of(float)))
    tau_err_median_windows: list[float] = field(
        factory=list, validator=validators.instance_of(list))
    tau_err_median_stdev: float | None = field(
        default=None, validator=validators.optional(
        validators.instance_of(float)))
    tau_err_median_stdev_windows: list[float] = field(
        factory=list, validator=validators.instance_of(list))
    
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
            timestep = self.model.get_timestep_by_type("molecular_dynamics")
            steps_mean = statistics.mean(self.lengths_steps)
            self.tau_mean = steps_mean * timestep
            steps_median = statistics.median(self.lengths_steps)
            self.tau_median = steps_median * timestep
            
            # Now handle the windows - a subset of the data up to a step
            self.window_lengths_steps = np.linspace(
                0, len(self.lengths_steps), self.num_windows+1, 
                dtype=int)[1:].tolist()
            for window_length_steps in self.window_lengths_steps:
                window_steps = self.lengths_steps[:window_length_steps]
                window_steps_mean = statistics.mean(window_steps)
                self.tau_mean_windows.append(window_steps_mean * timestep)
                window_steps_median = statistics.median(window_steps)
                self.tau_median_windows.append(window_steps_median * timestep)
        
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
            self.tau_err_mean = None
            self.tau_err_mean_windows = [None] * len(self.window_lengths_steps)
            self.tau_err_median = None
            self.tau_err_median_windows = [None] * len(self.window_lengths_steps)
            self.tau_err_mean_stdev = None
            self.tau_err_mean_stdev_windows = [None] * len(self.window_lengths_steps)
            self.tau_err_median_stdev = None
            self.tau_err_median_stdev_windows = [None] * len(self.window_lengths_steps)
            return
        assert error_analysis_method in ["bootstrap"]
        if error_analysis_method == "bootstrap":
            timestep = self.model.get_timestep_by_type("molecular_dynamics")
            self.tau_err_mean, self.tau_err_mean_stdev, self.tau_err_median, \
                self.tau_err_median_stdev = bootstrap_ramd(
                    num_error_samples=num_error_samples,
                    random_seed=random_seed,
                    timestep=timestep,
                    lengths_steps=self.lengths_steps
                )
            # Now handle the windows - a subset of the data up to a step
            self.window_lengths_steps = np.linspace(
                0, len(self.lengths_steps), self.num_windows+1, 
                dtype=int)[1:].tolist()
            for window_length_steps in self.window_lengths_steps:
                window_steps = self.lengths_steps[:window_length_steps]
                tau_err_mean, tau_err_mean_stdev, tau_err_median, \
                    tau_err_median_stdev = bootstrap_ramd(
                        num_error_samples=num_error_samples,
                        random_seed=random_seed,
                        timestep=timestep,
                        lengths_steps=window_steps
                    )
                self.tau_err_mean_windows.append(tau_err_mean)
                self.tau_err_mean_stdev_windows.append(tau_err_mean_stdev)
                self.tau_err_median_windows.append(tau_err_median)
                self.tau_err_median_stdev_windows.append(tau_err_median_stdev)

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
        if self.tau_mean is not None:
            print("Tau (ps) - mean:",
                  analyzers_base.pretty_string_value_error(
                    self.tau_mean, self.tau_err_mean_stdev))
        if self.tau_median is not None:
            print("Tau (ps) - median:",
                  analyzers_base.pretty_string_value_error(
                    self.tau_median, self.tau_err_median_stdev))
        
        return
    
    def make_image_directory(
            self, 
            image_directory: str
            ) -> None:
        """
        Create the directory and sub-directories for images to be 
        written. 
        """
        if image_directory != "" and not os.path.exists(image_directory):
            os.mkdir(image_directory)
        return
    
    def save_plots(
            self, 
            image_directory: str
            ) -> None:
        """
        Save a potentially useful series of plots of some quantities
        obtained during the analysis. Plots convergence.
        """
        # Mean plots
        filename = os.path.join(image_directory, "tau_mean_convergence.png")
        tau_mean_conv = np.array(self.tau_mean_windows)
        tau_err_mean_stdev_windows = np.array(self.tau_err_mean_stdev_windows)
        ret = analyzers_base.plot_scalar_conv(
            tau_mean_conv, self.window_lengths_steps, filename,
            title="$\\tau_{mean}$ Convergence", label="\\tau_{mean} (ps)",
            conv_errors=tau_err_mean_stdev_windows, return_plot=False)
        
        # Median plots
        filename = os.path.join(image_directory, "tau_median_convergence.png")
        tau_median_conv = np.array(self.tau_median_windows)
        tau_err_median_stdev_windows = np.array(self.tau_err_median_stdev_windows)
        ret = analyzers_base.plot_scalar_conv(
            tau_median_conv, self.window_lengths_steps, filename, 
            title="$\\tau_{median}$ Convergence", label="\\tau_{median} (ps)", 
            conv_errors=tau_err_median_stdev_windows, return_plot=False)
        return
        
    def scramble_me(
            self,
            random_seed: int | None = None) -> "Ramd_analyzer":
        """
        Create essentially a copy of this object, but with the data
        scrambled using sampling with replacement, and then re-calculate 
        the relevant quantities for the copy. Needed for doing error
        analysis of a data-series computed quantities, like Spearman rho.
        """
        random.seed(random_seed)
        new_analyzer = Ramd_analyzer(model=self.model, stage_index=self.stage_index)
        bootstrap_sample = random.choices(
            self.window_lengths_steps, k=len(self.window_lengths_steps))
        new_analyzer.lengths_steps = bootstrap_sample
        new_analyzer.process_data()
        new_analyzer.error_analysis(error_analysis_method="bootstrap", 
                                    random_seed=random_seed)
        return new_analyzer


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