"""
Base class for all experiment types. Handles a lot of common functionality
like systematically choosing a filepath and multiprocessing.

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

from dataclasses import dataclass
from abc import abstractmethod
from typing import List
from concurrent.futures import ProcessPoolExecutor
import os
from datetime import datetime

EXPERIMENTS_BASEPATH = "experiment_results"
ALL_COLORS = ["#00843D","#78BE20","#772583","#981D97","#FF6600","#FFA300","#00566e","#009CA6"]
FIT_CYCLE = ALL_COLORS[::2]

@dataclass(kw_only=True)
class ModelGenerationExperiment:
  # name of the experiment (determines the folder name)
  name: str
  # path to the reference data to fit
  ref_data_path: str

  # number of times the complete optimization process is repeated
  num_macroreps: int = 1
  # number of macroreps to run in parallel
  parallel_macroreps: int = 1

  def _prepare(self):
    """Prepare the experiment folder."""
    # path is base + subclass prefix + name + datetime
    if not hasattr(self, "shortname"):
        print("Error: Experiment subclass must specify shortname in __post_init__.")
        exit(-1)
    self.experiment_path = os.path.join(EXPERIMENTS_BASEPATH, self.shortname, self.name, str(datetime.now().isoformat(timespec="minutes")).replace(":", "_"))
    if os.path.exists(self.experiment_path):
      print("Error: Wait a minute, experiment path already exists. Not overwriting.")
      exit(-1)
    else:
      os.makedirs(self.experiment_path)
    with open(os.path.join(self.experiment_path, "hyperparameters.py"), "w") as file:
      file.write("""
from evolib.reaction import Reaction
from evolib.experiments.coupled_sindy_experiment import CoupledSindyExperiment
from evolib.experiments.evolving_libraries_experiment import EvolvingLibrariesExperiment
from evolib.experiments.random_search_experiment import RandomSearchExperiment
""")
      file.write("exp = ")
      file.write(str(self))
 
  # overload this with the experiment-specific logic
  @abstractmethod
  def experiment(self, macrorep_nr: int):
    pass

  # overload these methods with the experiment-specific logic
  @abstractmethod
  def evaluate_convergence(self, experiment_path: str, fig=None, ax=None):
    pass

  @abstractmethod
  def evaluate_fit(self, experiment_path: str, fig=None, ax=None):
    pass

  @abstractmethod
  def evaluate_model(self, experiment_path: str, fig=None, ax=None):
    pass

  def run(self):
    """Run the experiment as specified."""
    self._prepare()
    if self.parallel_macroreps > 1:
      pool = ProcessPoolExecutor(self.parallel_macroreps)
      pool.map(self.experiment, range(self.num_macroreps))
      pool.shutdown(wait=True)
    else:
      for mrep in range(self.num_macroreps):
        self.experiment(mrep)
    #self.evaluate(self.experiment_path)

