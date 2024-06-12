"""
List of experiments used in the paper. The corresponding variables have to
start with "exp_" in order to be discovered for running them later.

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

from evolib.experiments.coupled_sindy_experiment import CoupledSindyExperiment
from evolib.experiments.evolving_libraries_experiment import EvolvingLibrariesExperiment
from evolib.experiments.random_search_experiment import RandomSearchExperiment
from evolib.wnt import fixated_reactions as wnt_fixated_reactions

# sir

exp_sir_gp = EvolvingLibrariesExperiment(
  name="sir",
  ref_data_path="data/sir_100.csv",
  num_macroreps=10,
  parallel_macroreps=10,
  lib_size=2,     # ground truth model has 2 reactions
  max_num_left=2,
  max_num_right=3,
  pop_size=100,
  max_steps=100,
  slow_thresh=1e-6,
  co_prob=0.8,
  mut_prob=0.2,
  co_points=1,
  n_parents=10
)

exp_sir_rs = RandomSearchExperiment(
  name="sir",
  ref_data_path="data/sir_100.csv",
  num_macroreps=10,
  parallel_macroreps=10,
  lib_size=2,     # ground truth model has 2 reactions
  max_num_left=2,
  max_num_right=3,
  pop_size=100,
  max_steps=100,
  slow_thresh=1e-6 
)

exp_sir_csdiny = CoupledSindyExperiment(
  name="sir",
  ref_data_path="data/sir_100.csv",
  num_macroreps=1,
  parallel_macroreps=1,
  max_num_left=2,
  max_num_right=3,
  slow_thresh=1e-6 
)

# predatorprey

exp_predatorprey_gp = EvolvingLibrariesExperiment(
  name="predatorprey",
  ref_data_path="data/predatorprey_100.csv",
  num_macroreps=10,
  parallel_macroreps=10,
  lib_size=3,     # ground truth model has 3 reactions
  max_num_left=2,
  max_num_right=3,
  pop_size=100,
  max_steps=100,
  slow_thresh=1e-6,
  co_prob=0.8,
  mut_prob=0.2,
  co_points=1,
  n_parents=10
)

exp_predatorprey_rs = RandomSearchExperiment(
  name="predatorprey",
  ref_data_path="data/predatorprey_100.csv",
  num_macroreps=10,
  parallel_macroreps=10,
  lib_size=3,     # ground truth model has 3 reactions
  max_num_left=2,
  max_num_right=3,
  pop_size=100,
  max_steps=100,
  slow_thresh=1e-6
)

exp_predatorprey_csdiny = CoupledSindyExperiment(
  name="predatorprey",
  ref_data_path="data/predatorprey_100.csv",
  num_macroreps=1,
  parallel_macroreps=1,
  max_num_left=2,
  max_num_right=3,
  slow_thresh=1e-6
)

# wnt from scratch

exp_wnt_gp = EvolvingLibrariesExperiment(
  name="wnt",
  ref_data_path="data/wnt_0.csv",
  num_macroreps=10,
  parallel_macroreps=10,
  lib_size=60,     # ground truth model has 43 reactions
  max_num_left=2,
  max_num_right=2,
  pop_size=200,
  max_steps=10000,
  slow_thresh=1e-6,  # smallest rate in ground truth model is 1e-4
  co_prob=0.8,
  mut_prob=0.2,
  co_points=5,
  n_parents=10
)

exp_wnt_rs = RandomSearchExperiment(
  name="wnt",
  ref_data_path="data/wnt_0.csv",
  num_macroreps=10,
  parallel_macroreps=10,
  lib_size=60,     # ground truth model has 43 reactions
  max_num_left=2,
  max_num_right=2,
  pop_size=200,
  max_steps=10000,
  slow_thresh=1e-6  # smallest rate in ground truth model is 1e-4
)

exp_wnt_csindy = CoupledSindyExperiment(
  name="wnt",
  ref_data_path="data/wnt_0.csv",
  num_macroreps=1,
  parallel_macroreps=1,
  max_num_left=2,
  max_num_right=2,
  slow_thresh=1e-6 
)

# wnt extension

exp_wnt_ext_gp = EvolvingLibrariesExperiment(
  name="wnt_ext",
  ref_data_path="data/wnt_0.csv",
  num_macroreps=10,
  parallel_macroreps=10,
  lib_size=20,     # ground truth model has 43 reactions
  max_num_left=2,
  max_num_right=2,
  pop_size=200,
  max_steps=10000,
  slow_thresh=1e-6,  # smallest rate in ground truth model is 1e-4
  co_prob=0.8,
  mut_prob=0.2,
  co_points=5,
  n_parents=10,
  fixated_reactions=wnt_fixated_reactions
)

exp_wnt_ext_rs = RandomSearchExperiment(
  name="wnt_ext",
  ref_data_path="data/wnt_0.csv",
  num_macroreps=10,
  parallel_macroreps=10,
  lib_size=20,     # ground truth model has 43 reactions
  max_num_left=2,
  max_num_right=2,
  pop_size=200,
  max_steps=10000,
  slow_thresh=1e-6,  # smallest rate in ground truth model is 1e-4
  fixated_reactions=wnt_fixated_reactions
)

exp_wnt_ext_csindy = CoupledSindyExperiment(
  name="wnt_ext",
  ref_data_path="data/wnt_0.csv",
  num_macroreps=1,
  parallel_macroreps=1,
  max_num_left=2,
  max_num_right=2,
  slow_thresh=1e-6,
  fixated_reactions=wnt_fixated_reactions
)

# constrained wnt

exp_wnt_constrained_gp = EvolvingLibrariesExperiment(
  name="wnt_constrained",
  ref_data_path="data/wnt_0.csv",
  num_macroreps=10,
  parallel_macroreps=10,
  lib_size=60,     # ground truth model has 43 reactions
  max_num_left=2,
  max_num_right=2,
  pop_size=200,
  max_steps=10000,
  slow_thresh=1e-6,  # smallest rate in ground truth model is 1e-4
  co_prob=0.8,
  mut_prob=0.2,
  co_points=5,
  n_parents=10,
  wnt_model_constraints=True
)

exp_wnt_ext_constrained_gp = EvolvingLibrariesExperiment(
  name="wnt_ext_constrained",
  ref_data_path="data/wnt_0.csv",
  num_macroreps=10,
  parallel_macroreps=10,
  lib_size=20,     # ground truth model has 43 reactions
  max_num_left=2,
  max_num_right=2,
  pop_size=200,
  max_steps=10000,
  slow_thresh=1e-6,  # smallest rate in ground truth model is 1e-4
  co_prob=0.8,
  mut_prob=0.2,
  co_points=5,
  n_parents=10,
  fixated_reactions=wnt_fixated_reactions,
  wnt_model_constraints=True
)

if __name__ == "__main__":
    # run all experiments defined above
    to_run = []
    locals_copy = locals().copy()
    for local in locals_copy:
      if local.startswith("exp"):
        to_run.append(locals_copy[local])
    for e in to_run:
        print("Now running experiment", e)
        e.run()

