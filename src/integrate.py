"""
Functions to perform numerical integration for reaction networks.

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

from evolib.reaction import Reaction
from evolib.reaction_library import ReactionLibrary
from evolib.sindy.util import prepare_sindy_data

import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import pandas as pd
import numpy as np
from tqdm import tqdm


def gen_data_for(rs: ReactionLibrary, init_val: np.array, t_0=0.0, t_end=1.0, n_points=100):
  """Returns a pandas dataset of species concentrations over time."""
  if len(rs.reactions) > 0 and rs.reactions[0].species_names != []:
    sim_data = pd.DataFrame(columns=["time"] + rs.reactions[0].species_names)
  else:
    sim_data = pd.DataFrame(columns=["time"] + list(map(lambda x: f"S{x}", range(0, rs.num_species))))
  solver = "LSODA"
  if n_points > 0:
    print("Integrating from", t_0, "to", t_end, "with", solver, "measuring at", n_points, "points...")
    with tqdm(total=n_points, unit="%") as pbar:
      res = solve_ivp(rs.apply_odes, (t_0, t_end), init_val, solver, np.linspace(t_0, t_end, n_points), args=[pbar, [t_0, (t_end-t_0)/n_points]])#, rtol=1e-3, atol=1e-6)
  else:
    res = solve_ivp(rs.apply_odes, (t_0, t_end), init_val, solver, None)#, rtol=1e-3, atol=1e-6)
  if not res.success:
    print("ERROR")
    print(res)
    exit(-1)
  sim_data.time = res.t
  for idx, col in enumerate(sim_data.columns[1:]):
    sim_data[col] = res.y[idx]
  return sim_data

def hyperparams_from_data(path):
  """Helper to retrieve some hyperparameters from an existing dataset."""
  data = pd.read_csv(path)
  init_val = data[[col for col in data.columns if col != "time"]].iloc[0]
  t_0 = data.time.iloc[0]
  t_end = data.time.iloc[-1]
  npoints = len(data.time) - 1
  return init_val, t_0, t_end, npoints

def plot_sim_trace_of(rs: ReactionLibrary, n_points=100, t_end=None):
  """Plots simulated concentration trajectories over time."""
  ref_data = pd.read_csv(rs.ref_data_path)
  init_val, t_0, _t_end, _ = hyperparams_from_data(rs.ref_data_path)
  if t_end is None:
    t_end = _t_end
  sim_data = gen_data_for(rs, init_val, t_0, t_end, n_points)
  fig, ax = plt.subplots()
  ref_data.add_prefix("ref_").plot(ax=ax, x="ref_time", marker="o", lw=0)
  plt.gca().set_prop_cycle(None)
  sim_data.add_prefix("sim_").plot(ax=ax, x="sim_time", ls="--")
  plt.xlabel("time")
  plt.ylabel("amount")
  #plt.show()

def plot_system_comparison(rs1: ReactionLibrary, rs2: ReactionLibrary, ref_data_path: str, n_points=100, rs1_label="sim1", rs2_label="sim2"):
  """Plots simulated concentration trajectories over time for two systems for comparison."""
  init_val_1, t_0_1, t_end_1, _ = hyperparams_from_data(ref_data_path)
  init_val_2, t_0_2, t_end_2, _ = hyperparams_from_data(ref_data_path)
  if rs1.ref_data_path != "":
    init_val_1, t_0_1, t_end_1, _ = hyperparams_from_data(rs1.ref_data_path)
  if rs2.ref_data_path != "":
    init_val_2, t_0_2, t_end_2, _ = hyperparams_from_data(rs2.ref_data_path)
  ref_data = pd.read_csv(ref_data_path)
  sim1_data = gen_data_for(rs1, init_val_1, t_0_1, t_end_1, n_points)
  sim2_data = gen_data_for(rs2, init_val_2, t_0_2, t_end_2, n_points)
  fig, ax = plt.subplots()
  ref_data.add_prefix("ref_").plot(ax=ax, x="ref_time", marker="o", lw=0)
  plt.gca().set_prop_cycle(None)
  sim1_data.add_prefix(rs1_label+"_").plot(ax=ax, x=f"{rs1_label}_time", ls="--")
  plt.gca().set_prop_cycle(None)
  sim2_data.add_prefix(rs2_label+"_").plot(ax=ax, x=f"{rs2_label}_time", ls="-")
  return fig, ax

def plot_der_trace_of(rs: ReactionLibrary, n_points=100, t_end=None):
  """Plots derivatives over time."""
  data, times = prepare_sindy_data(rs.ref_data_path)
  data_deriv = np.gradient(data, times, axis=0, edge_order=2)
  ref_data = pd.read_csv(rs.ref_data_path)
  ref_data[ref_data.columns[1:]] = data_deriv
  sim_deriv = np.array([sum(r.get_deriv_terms(d) * r.rate for r in rs.reactions) for d in data])
  sim_data = ref_data.copy()
  sim_data[sim_data.columns[1:]] = sim_deriv
  fig, ax = plt.subplots()
  ref_data.add_prefix("ref_").plot(ax=ax, x="ref_time", marker="o", lw=0)
  plt.gca().set_prop_cycle(None)
  sim_data.add_prefix("sim_").plot(ax=ax, x="sim_time", ls="--")
  plt.xlabel("time")
  plt.ylabel(r"$\frac{d amount}{d t}$")
  #plt.show()

 
if __name__ == "__main__":
  # just for testing / requires previous execution of gen_data.py
  plot_sim_trace_of(
    ReactionLibrary([
      Reaction([0], [0, 0], 15, 2),
      Reaction([1], [], 12, 2),
      Reaction([0, 1], [1, 1], 0.01, 2),
    ], 2, "data/predatorprey_10.csv"),
    n_points=1000
  )

