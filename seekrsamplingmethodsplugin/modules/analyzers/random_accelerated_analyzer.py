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
DEFAULT_NUM_WINDOWS = 5


def bootstrap_ramd(
        lengths_steps: list[float],
        num_error_samples: int | None,
        random_seed: int | None,
        timestep: float,
        ) -> tuple[float | None, float | None, float | None, float | None]:
    """
    Bootstrap the mean and median escape times.

    Returns
    -------
    tau_err_mean, tau_err_mean_stdev, tau_err_median, tau_err_median_stdev
        Mean/median of the bootstrap distribution of the mean/median, and the
        corresponding sample standard deviations (all in the same time units
        as ``timestep``). Returns Nones if fewer than 2 escape times.
    """
    if len(lengths_steps) < 2:
        return None, None, None, None

    bootstrap_mean_steps = []
    bootstrap_median_steps = []
    n_samples = num_error_samples or DEFAULT_NUM_ERROR_SAMPLES
    for i in range(n_samples):
        if random_seed is not None:
            random.seed(random_seed + i)
        bootstrap_sample = random.choices(
            lengths_steps, k=len(lengths_steps))
        bootstrap_mean_steps.append(statistics.mean(bootstrap_sample))
        bootstrap_median_steps.append(statistics.median(bootstrap_sample))

    tau_err_mean = statistics.mean(bootstrap_mean_steps) * timestep
    tau_err_mean_stdev = statistics.stdev(bootstrap_mean_steps) * timestep
    tau_err_median = statistics.median(bootstrap_median_steps) * timestep
    tau_err_median_stdev = statistics.stdev(bootstrap_median_steps) * timestep

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
        self.tau_mean_windows = []
        self.tau_median_windows = []
        self.window_lengths_steps = []
        if self.lengths_steps:
            timestep = self.model.get_timestep_by_type("molecular_dynamics")
            steps_mean = statistics.mean(self.lengths_steps)
            self.tau_mean = steps_mean * timestep
            steps_median = statistics.median(self.lengths_steps)
            self.tau_median = steps_median * timestep

            # Cumulative windows over the ordered escape-time list: the i-th
            # window uses the first n_i trajectories (file / swarm order).
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
        self.tau_err_mean_windows = []
        self.tau_err_mean_stdev_windows = []
        self.tau_err_median_windows = []
        self.tau_err_median_stdev_windows = []

        if len(self.lengths_steps) < 2:
            # Not enough data for error analysis
            self.tau_err_mean = None
            self.tau_err_median = None
            self.tau_err_mean_stdev = None
            self.tau_err_median_stdev = None
            n_win = len(self.window_lengths_steps)
            self.tau_err_mean_windows = [None] * n_win
            self.tau_err_median_windows = [None] * n_win
            self.tau_err_mean_stdev_windows = [None] * n_win
            self.tau_err_median_stdev_windows = [None] * n_win
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
            # Ensure window sizes exist even if process_data was not called.
            if not self.window_lengths_steps:
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
        analyzers_base.plot_scalar_conv(
            tau_mean_conv, self.window_lengths_steps, filename,
            title="$\\tau_{mean}$ Convergence", label="\\tau_{mean} (ps)",
            conv_errors=tau_err_mean_stdev_windows, return_plot=False)

        # Median plots
        filename = os.path.join(image_directory, "tau_median_convergence.png")
        tau_median_conv = np.array(self.tau_median_windows)
        tau_err_median_stdev_windows = np.array(self.tau_err_median_stdev_windows)
        analyzers_base.plot_scalar_conv(
            tau_median_conv, self.window_lengths_steps, filename,
            title="$\\tau_{median}$ Convergence", label="\\tau_{median} (ps)",
            conv_errors=tau_err_median_stdev_windows, return_plot=False)
        return

    def scramble_me(
            self,
            random_seed: int | None = None,
            n_trajectories: int | None = None,
            num_error_samples: int | None = None,
            skip_error_analysis: bool = True,
            ) -> "Ramd_analyzer":
        """
        Bootstrap-resample escape times and recompute tau statistics.

        Samples **with replacement** from ``lengths_steps`` (the escape
        times), not from ``window_lengths_steps`` (window sizes). Optional
        ``n_trajectories`` restricts the population to the first N escapes
        in file order, which is the correct uncertainty for a cumulative
        sampling window used in ranking-metric convergence plots.

        Parameters
        ----------
        skip_error_analysis :
            If True (default), only recompute point estimates. Nested
            bootstraps are expensive and usually unnecessary when this
            method is itself used inside an outer bootstrap over compounds.
        """
        random.seed(random_seed)
        new_analyzer = Ramd_analyzer(
            model=self.model,
            stage_index=self.stage_index,
            num_windows=self.num_windows,
        )
        if n_trajectories is None:
            population = list(self.lengths_steps)
        else:
            if n_trajectories < 1:
                raise ValueError(
                    "n_trajectories must be >= 1 when provided.")
            population = list(self.lengths_steps[:n_trajectories])
        if len(population) == 0:
            new_analyzer.lengths_steps = []
            return new_analyzer

        new_analyzer.lengths_steps = random.choices(
            population, k=len(population))
        new_analyzer.process_data()
        if (not skip_error_analysis) and len(new_analyzer.lengths_steps) >= 2:
            new_analyzer.error_analysis(
                error_analysis_method="bootstrap",
                num_error_samples=num_error_samples,
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
