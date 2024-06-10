from evolib.experiments.model_generation_experiment import ModelGenerationExperiment, FIT_CYCLE
from evolib.evolution.evolve_libraries import EvolvingLibraries, optimize
from evolib.integrate import gen_data_for, hyperparams_from_data
from evolib.reaction import Reaction

from typing import List
from dataclasses import dataclass, field
import os

import pandas as pd
import matplotlib.pyplot as plt


@dataclass(kw_only=True)
class RandomSearchExperiment(ModelGenerationExperiment):
  lib_size: int      # maximum number of reactions in a single system
  max_num_left: int  # maximum number of reactands
  max_num_right: int # maximum number of products
  pop_size: int      # number of libraries to evolve
  max_steps: int     # maximum number of steps to perform
  slow_thresh: float # cutoff threshold for slow reactions
  fixated_reactions: List[Reaction] = field(default_factory=lambda: [])

  def __post_init__(self):
    self.shortname = "rs"

  def experiment(self, mrep: int):
    opt = EvolvingLibraries(self.lib_size, self.max_num_left, self.max_num_right, self.ref_data_path, pop_size=self.pop_size, fixated_reactions=self.fixated_reactions, cpu_div=self.parallel_macroreps, rs_mode=True)
    res = optimize(opt, nsteps=self.max_steps, output_filename=os.path.join(self.experiment_path, f"rs_loss_rep-{mrep:03d}.csv"))
    opt.overall_fittest.clean_slow_reactions(self.slow_thresh)
    with open(os.path.join(self.experiment_path, f"rs_model_rep-{mrep:03d}.txt"), "w") as file:
      file.write(str(opt.overall_fittest))
    rs_data = gen_data_for(opt.overall_fittest, *hyperparams_from_data(opt.overall_fittest.ref_data_path))
    rs_data.to_csv(os.path.join(self.experiment_path, f"rs_model_data_rep-{mrep:03d}.csv"), index=False)

  def evaluate_convergence(self, experiment_path: str, fig=None, ax=None):
    if fig is None or ax is None:
      fig, ax = plt.subplots()
    data = [pd.read_csv(os.path.join(experiment_path, f"rs_loss_rep-{rep:03d}.csv")).drop(columns="model") for rep in range(self.num_macroreps)]
    mean_loss = pd.concat(data).groupby("step").overall_fittest_norm.mean()
    loss = [d.overall_fittest_norm for d in data]
    ax.plot(mean_loss, label="random search (mean)", color="blue", ls="--")
    for l in loss:
      ax.plot(l, label="", alpha=0.1, color="blue", ls="--")
    if fig is None or ax is None:
      outpath = os.path.join(experiment_path, "rs_loss.pdf")
      print("writing", outpath)
      fig.savefig(outpath)
    return fig, ax

  def evaluate_fit(self, experiment_path: str, fig=None, ax=None):
    if fig is None or ax is None:
      fig, ax = plt.subplots()
    data = [pd.read_csv(os.path.join(experiment_path, f"rs_model_data_rep-{rep:03d}.csv")) for rep in range(self.num_macroreps)]
    sim_data = pd.concat(data).groupby("time").mean().reset_index()
    ref_data = pd.read_csv(self.ref_data_path)
    ax.set_prop_cycle("color", FIT_CYCLE)
    ref_data.add_prefix("ref_").plot(ax=ax, x="ref_time", marker="+", ls="None", markersize=3)
    ax.set_prop_cycle("color", FIT_CYCLE)
    sim_data.add_prefix("rs_").plot(ax=ax, x="rs_time", ls="--", alpha=0.5)
    return fig, ax

  def evaluate_model(self, experiment_path: str, fig=None, ax=None):
    if fig is None or ax is None:
      fig, ax = plt.subplots(figsize=(1,1))
    data = [pd.read_csv(os.path.join(experiment_path, f"rs_loss_rep-{rep:03d}.csv")) for rep in range(self.num_macroreps)]
    data = pd.concat(data)
    model = data[data.overall_fittest == data.overall_fittest.max()].model.iloc[-1]
    from evolib.reaction_library import ReactionLibrary
    from evolib.reaction import Reaction
    model = eval(model)
    model.clean_slow_reactions(1e-3)
    ax.set_axis_off()
    ax.text(0.0, 0.0, model.to_latex())
    return fig, ax


