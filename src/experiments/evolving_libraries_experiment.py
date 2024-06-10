"""
Implementation of the experiment base class for testing evolving libraries.
Contains all the setup and evaluation code.

MIT License

Copyright (c) 2024 Justin Kreikemeyer

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

from evolib.experiments.model_generation_experiment import ModelGenerationExperiment, FIT_CYCLE
from evolib.evolution.evolve_libraries import EvolvingLibraries, optimize
from evolib.integrate import gen_data_for, hyperparams_from_data
from evolib.reaction import Reaction
from evolib.wnt import wnt_penalty, species_subgroups

from dataclasses import dataclass, field
from typing import List
import os

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


@dataclass(kw_only=True)
class EvolvingLibrariesExperiment(ModelGenerationExperiment):
  lib_size: int      # maximum number of reactions in a single system
  max_num_left: int  # maximum number of reactands
  max_num_right: int # maximum number of products
  pop_size: int      # number of libraries to evolve
  max_steps: int     # maximum number of steps to perform
  slow_thresh: float # cutoff threshold for slow reactions
  co_prob: float = 0.75
  co_points: int = 4
  mut_prob: float = 0.25
  n_parents: int = 3
  fixated_reactions: List[Reaction] = field(default_factory=lambda: [])
  wnt_model_constraints: bool = False

  def __post_init__(self):
    self.shortname = "evolib"

  def experiment(self, mrep: int):
    extra_args = dict()
    if self.wnt_model_constraints:
      extra_args = dict(subgroups=species_subgroups, penalty_func=wnt_penalty)
    opt = EvolvingLibraries(
       self.lib_size,
       self.max_num_left,
       self.max_num_right,
       self.ref_data_path,
       pop_size=self.pop_size,
       co_prob=self.co_prob,
       mut_prob=self.mut_prob,
       co_points=self.co_points,
       n_parents=self.n_parents,
       fixated_reactions=self.fixated_reactions,
       cpu_div=self.parallel_macroreps,
       **extra_args
    )
    res = optimize(opt, nsteps=self.max_steps, output_filename=os.path.join(self.experiment_path, f"evolib_loss_rep-{mrep:03d}.csv"))
    opt.overall_fittest.clean_slow_reactions(self.slow_thresh)
    with open(os.path.join(self.experiment_path, f"evolib_model_rep-{mrep:03d}.txt"), "w") as file:
      file.write(str(opt.overall_fittest))
    evolib_data = gen_data_for(opt.overall_fittest, *hyperparams_from_data(opt.overall_fittest.ref_data_path))
    evolib_data.to_csv(os.path.join(self.experiment_path, f"evolib_model_data_rep-{mrep:03d}.csv"), index=False)

  def evaluate_convergence(self, experiment_path: str, fig=None, ax=None):
    if fig is None or ax is None:
      fig, ax = plt.subplots()
    data = [pd.read_csv(os.path.join(experiment_path, f"evolib_loss_rep-{rep:03d}.csv")).drop(columns="model") for rep in range(self.num_macroreps)]
    mean_loss = pd.concat(data).groupby("step").overall_fittest_norm.mean()
    loss = [d.overall_fittest_norm for d in data]
    ax.plot(mean_loss, label="evolib (mean)", color="black")
    for l in loss:
      ax.plot(l, label="", alpha=0.1, color="black", ls="-")
    if fig is None or ax is None:
      outpath = os.path.join(experiment_path, "evolib_loss.pdf")
      print("writing", outpath)
      fig.savefig(outpath)
    return fig, ax

  def evaluate_fit(self, experiment_path: str, fig=None, ax=None):
    if fig is None or ax is None:
      fig, ax = plt.subplots()
    data = [pd.read_csv(os.path.join(experiment_path, f"evolib_model_data_rep-{rep:03d}.csv")) for rep in range(self.num_macroreps) if os.path.exists(os.path.join(experiment_path, f"evolib_model_data_rep-{rep:03d}.csv"))]
    loss_data = [pd.read_csv(os.path.join(experiment_path, f"evolib_loss_rep-{rep:03d}.csv")).iloc[-1]["overall_fittest"] for rep in range(self.num_macroreps) if os.path.exists(os.path.join(experiment_path, f"evolib_model_data_rep-{rep:03d}.csv"))]
    sim_data = data[np.argmax(loss_data)] 
    ref_data = pd.read_csv(self.ref_data_path)
    ax.set_prop_cycle("color", FIT_CYCLE)
    ref_data.add_prefix("ref ").iloc[::2].plot(ax=ax, x="ref time", marker="+", ls="None", markersize=2, alpha=0.5)
    ax.set_prop_cycle("color", FIT_CYCLE)
    sim_data.add_prefix("evolib ").plot(ax=ax, x="evolib time", ls="-", alpha=1.0)
    return fig, ax

  def evaluate_model(self, experiment_path: str, fig=None, ax=None):
    if fig is None or ax is None:
      fig, ax = plt.subplots(figsize=(0.5,0.5))
    data = [pd.read_csv(os.path.join(experiment_path, f"evolib_loss_rep-{rep:03d}.csv")) for rep in range(self.num_macroreps)]
    data = pd.concat(data)
    model = data[data.overall_fittest == data.overall_fittest.max()].model.iloc[-1]
    from evolib.reaction_library import ReactionLibrary
    from evolib.reaction import Reaction
    model = eval(model)
    original_length = len(model.reactions)
    model.clean_slow_reactions(self.slow_thresh)
    n_more = original_length - len(model.reactions)
    ax.set_axis_off()
    from evolib.wnt import wnt_model
    ax.text(0.0, 0.0, model.to_latex(wnt_model) + (f"\\\\(+ {n_more} more)" if n_more > 0 else ""))
    return fig, ax


