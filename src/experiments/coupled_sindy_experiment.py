"""
Implementation of the experiment base class for testing coupled SINDy.
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
from evolib.reaction_library import ReactionLibrary, species_from_csv
from evolib.reaction import Reaction
from evolib.sindy.coupled_sindy import CoupledSINDy
from evolib.reaction_enumerator import ReactionEnumerator
from evolib.integrate import gen_data_for, hyperparams_from_data

from dataclasses import dataclass, field
from typing import List
import os
import csv

import pandas as pd
import matplotlib.pyplot as plt


@dataclass(kw_only=True)
class CoupledSindyExperiment(ModelGenerationExperiment):
  max_num_left: int  # maximum number of reactands
  max_num_right: int # maximum number of products
  slow_thresh: float 
  method: str = "nnls"
  fixated_reactions: List[Reaction] = field(default_factory=lambda: [])

  def __post_init__(self):
    self.shortname = "csindy"

  def experiment(self, mrep: int):
    num_species, species_names = species_from_csv(self.ref_data_path)
    reactions = list(ReactionEnumerator(num_species, self.max_num_left, self.max_num_right, 1, reaction_blacklist=self.fixated_reactions, species_names=species_names).generator())
    library = ReactionLibrary(reactions, num_species, ref_data_path=self.ref_data_path, species_names=species_names, fixated_reactions=self.fixated_reactions)
    csindy = CoupledSINDy(library, verbose=True)
    res = csindy.optimize(method="nnls")
    model = csindy.get_sparse_reaction_library(res, thresh=self.slow_thresh)
    with open(os.path.join(self.experiment_path, f"csindy_model_rep-{mrep:03d}.txt"), "w") as file:
      file.write(str(model))
    with open(os.path.join(self.experiment_path, f"csindy_loss_rep-{mrep:03d}.csv"), "w") as file:
      writer = csv.writer(file)
      writer.writerow("step,r2_score,norm,model".split(","))
      writer.writerow([0,res.r2_score,res.norm,model.__repr__()])
    csindy_data = gen_data_for(model, *hyperparams_from_data(model.ref_data_path))
    csindy_data.to_csv(os.path.join(self.experiment_path, f"csindy_model_data_rep-{mrep:03d}.csv"), index=False)

  def evaluate_convergence(self, experiment_path: str, fig=None, ax=None):
    if fig is None or ax is None:
      fig, ax = plt.subplots()
    data = pd.read_csv(os.path.join(experiment_path, f"csindy_loss_rep-000.csv"))
    ax.hlines(data.norm, ax.get_xlim()[0], ax.get_xlim()[1], label="c-SINDy", color="red")
    if fig is None or ax is None:
      outpath = os.path.join(experiment_path, "csindy_loss.pdf")
      print("writing", outpath)
      fig.savefig(outpath)
    return fig, ax

  def evaluate_fit(self, experiment_path: str, fig=None, ax=None):
    if fig is None or ax is None:
      fig, ax = plt.subplots()
    if os.path.exists(os.path.join(experiment_path, f"csindy_model_data_rep-000.csv")):
      data = pd.read_csv(os.path.join(experiment_path, f"csindy_model_data_rep-000.csv"))
      ax.set_prop_cycle("color", FIT_CYCLE)
      data.add_prefix("c-sindy ").plot(ax=ax, x="c-sindy time", ls=":", alpha=0.5)
    return fig, ax

  def evaluate_model(self, experiment_path: str, fig=None, ax=None):
    if fig is None or ax is None:
      fig, ax = plt.subplots(figsize=(1,1))
    data = pd.read_csv(os.path.join(experiment_path, f"csindy_loss_rep-000.csv"))
    model = data[data.r2_score == data.r2_score.max()].model.iloc[-1]
    from evolib.reaction_library import ReactionLibrary
    from evolib.reaction import Reaction
    model = eval(model)
    original_length = len(model.reactions)
    model.clean_slow_reactions(self.slow_thresh)
    n_more = original_length - len(model.reactions)
    ax.set_axis_off()
    if len(model.reactions) > 50:
      ax.text(0.0, 0.0, "incomprehensible") 
    else:
      ax.text(0.0, 0.0, model.to_latex() + (f"\\\\(+ {n_more} more)" if n_more > 0 else ""))
    return fig, ax

