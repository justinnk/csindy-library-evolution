"""
Class and helper functions for evolvling libraries.
Applies a genetic algorithm to search for optimal sublabraries to use
for sparse regression with coupled-SINDy.

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
from evolib.reaction_library import ReactionLibrary, species_from_csv, get_model_intersection
from evolib.reaction_enumerator import ReactionEnumerator
from evolib.integrate import gen_data_for, plot_system_comparison, hyperparams_from_data, plot_sim_trace_of, plot_der_trace_of
from evolib.sindy.coupled_sindy import CoupledSINDy, optimize_coupled_sindy
from evolib.wnt import fixated_reactions as wnt_fixated_reactions
from evolib.wnt import species_subgroups, wnt_penalty
from evolib.model_groundtruths import sir, predatorprey, wnt_model

from typing import List
import os
from multiprocessing import Pool
from copy import deepcopy
from random import choices
import signal
import matplotlib.pyplot as plt
from collections import defaultdict
import csv

import pandas as pd


# from: https://stackoverflow.com/questions/18499497/how-to-process-sigterm-signal-gracefully
class GracefulKiller:
  kill_now = False
  def __init__(self):
    signal.signal(signal.SIGINT, self.exit_gracefully)
    signal.signal(signal.SIGTERM, self.exit_gracefully)

  def exit_gracefully(self, signum, frame):
    self.kill_now = True

  def reset(self):
    self.kill_now = False


class EvolvingLibraries:
  """
  Implementation of a genetic algorithm to evolve libraries of reactions
  in order to find the optimal subsample to be used with sparse regression
  in the coupled SINDy approach.
  """
  def __init__(self,
               # Reaction System
               num_reactions: int,              # (maximum) number of reactions per library
               max_num_left: int,               # maximum number of species on the left side of reactions
               max_num_right: int,              # maximum number of species on the right side of reactions
               ref_data_path: str,              # path to the reference data set
               *,
               fixated_reactions: list = None,  # list of known reactions that are included in the model
               subgroups: list = [],            # list of subgroups in which the species interact
               # GP
               pop_size: int = 10,                   # number of candidate systems in the GP
               co_prob: float = 0.2,                # crossover probability
               mut_prob: float = 0.8,               # mutation probability
               n_parents: int = 3,                   # number of (best) parent soltions picked for reproduction
               co_points: int = 4,                   # number of reactions swapped by crossover
               constraints: List["Constraint"] = [], # List of constraints on the model structure
               penalty_func = None,                 # function taking a ReactionSystem and returning an integer
               verbose: bool = False,
               rs_mode: bool = False,           # perform random search instead of evolution
               cpu_div: int = 1                 # by how much to divide the available cpu cores (used for parallel experiments)
  ):
    # reaction system
    self.num_reactions = num_reactions
    self.max_num_left = max_num_left
    self.max_num_right = max_num_right
    self.ref_data_path = ref_data_path
    self.num_species, self.species_names = species_from_csv(ref_data_path)
    self.fixated_reactions = fixated_reactions if fixated_reactions is not None else []
    # gp 
    self.verbose = verbose
    self.pop_size = pop_size
    self.co_prob = co_prob
    self.mut_prob = mut_prob
    if n_parents > pop_size and not rs_mode:
      self.log(f"Error: n_parents ({n_parents}) >= pop_size ({pop_size})")
      exit(-1)
    self.n_parents = n_parents
    if co_points >= num_reactions and not rs_mode:
      self.log(f"Error: co_points ({co_points}) >= num_reactions ({num_reactions})")
      exit(-1)
    self.co_points = co_points
    self.population = []
    self.constraints = constraints
    self.penalty_func = penalty_func
    self.rs_mode = rs_mode
    self.pool = Pool(os.cpu_count() // cpu_div)
    # internal state
    self._enumerator = ReactionEnumerator(self.num_species, self.max_num_left, self.max_num_right, 1, self.constraints, shuffle=True, reaction_blacklist=self.fixated_reactions, species_names=self.species_names, subgroups=subgroups)
    self.log(f"Library size: {self._enumerator.get_number()}")
    self.fittest = None # fittest in the current gen
    self.overall_fittest = None # overall fittest
    self.overall_fittest_value = float("-inf") # fitness of overall fittest
    self.overall_fittest_norm = float("inf") # fitness of overall fittest
    self.no_improvement_steps = 0 # steps without improvement of the fittest
    self.init_population()
    assert self.pop_size > 3

  def log(self, *msg: str):
    if self.verbose: print("[EvolvingLibraries]", *msg)

  def print_population(self):
    for idx, member in enumerate(self.population):
      print("="*10, idx, "="*10)
      print(member)

  def init_population(self):
    """Initialize population by repeatedly mutating empty reactions."""
    # random initial population
    self.population = []
    for idx in range(self.pop_size):
      self.population.append(
        #ReactionLibrary.from_ref_data([self._enumerator.get_random() for _ in range(self.num_reactions)], self.ref_data_path, self._enumerator)
        ReactionLibrary([self._enumerator.get_random() for _ in range(self.num_reactions)], self.num_species, ref_data_path=self.ref_data_path, species_names=self.species_names, enumerator=self._enumerator, fixated_reactions=self.fixated_reactions)
      )

  @classmethod
  def get_fitness(cls, *args) -> float: #member: ReactionLibrary, penalty_func) -> float:
    """Fit model to data as close as possible and return loss."""
    member, penalty_func = args[0]
    #self.log("fitting library...")
    res = CoupledSINDy(member).optimize()
    penalty = 0
    if penalty_func != None:
      penalty = penalty_func(member)
    for idx, r in enumerate(member.non_fixated_reactions):
      r.rate = res.coef[idx]
    return res.r2_score - penalty, res.norm, member

  def step(self, gen: int) -> float:
    # determine fitness
    self.log("fitness evaluation...")
    fitness_list = self.pool.map(self.get_fitness, zip(self.population, [self.penalty_func for _ in range(self.pop_size)]))
    #fitness_list = [self.get_fitness([m, self.penalty_func]) for m in self.population]
    # debugging
    self.log("Population Gen", gen)
    self.log("="*10)
    if self.verbose: self.print_population()
    # sort by fitness
    fitness_sorted = sorted(fitness_list, key=lambda x: x[0], reverse=True)
    self.log(fitness_sorted)
    self.log(f"Fittest individual in gen {gen} is {fitness_sorted[0][2]} with {fitness_sorted[0][0]}")
    # update (overall) fittest
    self.fittest = deepcopy(fitness_sorted[0][2])
    if fitness_sorted[0][0] > self.overall_fittest_value:
      self.overall_fittest_value = fitness_sorted[0][0]
      self.overall_fittest_norm = fitness_sorted[0][1]
      self.overall_fittest = deepcopy(self.fittest)

    # if random search model is enabled, don't apply genetic operators
    # and continue with a new random generation (just copying over the fittest)
    if self.rs_mode:
      self.init_population()
      self.population[0] = deepcopy(self.overall_fittest)
    else:
      new_pop = []
      # copy fittest (elitism)
      new_pop.append(deepcopy(fitness_sorted[0][2])) 
      # mutate or cross over the fittest parents to fill new population
      parents = [fitness_sorted[idx][2] for idx in range(0, self.n_parents)]
      for _ in range(self.pop_size - 1):
        repro_choice = choices([0, 1], weights=[self.co_prob, self.mut_prob], k=1)[0]
        # crossover
        if repro_choice == 0:
          member = deepcopy(choices(parents)[0])
          member2 = deepcopy(choices(parents)[0])
          for _ in range(self.co_points):
            member.crossover(member2)
        # mutation
        elif repro_choice == 1:
          member = deepcopy(choices(parents)[0])
          member.mutate(self.constraints)
        else:
          raise Exception("Invalid repro operation.")
        # reset reaction constants
        for reac in member.non_fixated_reactions: 
          reac.rate = 0.0
        # add new member
        new_pop.append(member)
      self.population = new_pop
    return fitness_sorted[0][0], self.overall_fittest_value, self.overall_fittest_norm, self.overall_fittest

def optimize(opt: EvolvingLibraries, fitness_thresh: float = 0.999999999, nsteps: int = 100, output_filename: str = "gp_loss.csv", killer: GracefulKiller = None):
  """Wrapper to perform optimization using opt and writing the progress to a csv file."""
  print("Initial Population")
  print("="*10)
  opt.print_population()
  if os.path.dirname(output_filename) != "" and not os.path.exists(os.path.dirname(output_filename)):
    os.makedirs(os.path.dirname(output_filename))
  with open(output_filename, "w") as file:
    writer = csv.writer(file)
    writer.writerow("step curr_fittest overall_fittest overall_fittest_norm model".split(" "))
    for s in range(nsteps):
      curr_fittest, overall_fittest, overall_fittest_norm, overall_fittest_model = opt.step(s)
      print(f"{s},{curr_fittest},{overall_fittest},{overall_fittest_norm}")
      writer.writerow([s, curr_fittest, overall_fittest, overall_fittest_norm, overall_fittest_model.__repr__()])
      file.flush()
      if curr_fittest >= fitness_thresh: break
      if killer is not None and killer.kill_now: break
  print("="*10)
  print("Fittest in last gen:", opt.fittest.print_model())
  print("Fittest overall:", opt.overall_fittest.print_model())


if __name__ == "__main__":
  # for easy testing
  ref_model = "SIR" # <-- change this
  # (if you choose WNT, it may take quite a while to get results)

  if ref_model == "SIR":
    ref_data_path = "data/sir_0.csv"
    gt_model = sir
  elif ref_model == "PP":
    ref_data_path = "data/predatorprey_0.csv"
    gt_model = predatorprey
  elif ref_model == "WNT":
    ref_data_path = "data/wnt_0.csv"
    gt_model = wnt_model
  else:
    print("Unknown reference.")
    exit(-1)

  killer = GracefulKiller()
  # shared arguments for optimizer <--- change these to your liking/use case!
  #opt_args = dict(num_reactions=10, max_num_left=2, max_num_right=2, ref_data_path=ref_data_path, pop_size=1000, co_points=2, mut_prob=0.2, co_prob=0.8, n_parents=5, fixated_reactions=wnt_fixated_reactions, cpu_div=2, penalty_func=wnt_penalty, subgroups=species_subgroups)
  opt_args = dict(num_reactions=2, max_num_left=2, max_num_right=3, ref_data_path=ref_data_path, pop_size=20, co_points=1, mut_prob=0.2, co_prob=0.8, n_parents=3, verbose=False)
  opt = EvolvingLibraries(**opt_args)
  res = optimize(opt, killer=killer, nsteps=200)
  ref_data = pd.read_csv(ref_data_path)
  #plot_sim_trace_of(opt.overall_fittest)
  print("GP Model")
  opt.overall_fittest.clean_slow_reactions(0)
  print(opt.overall_fittest)
  print("Ground Truth")
  print(gt_model)
  print("Intersection")
  print(get_model_intersection(gt_model, opt.overall_fittest))
  plot_sim_trace_of(opt.overall_fittest)
  plt.show()
  plot_der_trace_of(opt.overall_fittest)
  plt.show()
  gp_data = gen_data_for(opt.overall_fittest, *hyperparams_from_data(opt.overall_fittest.ref_data_path))
  gp_data.to_csv("gp_model_data.csv", index=False)

  killer.reset()
  
  # random search using the same opt_args
  opt_rs = EvolvingLibraries(**opt_args, rs_mode=True)
  res = optimize(opt_rs, killer=killer, nsteps=200)
  #plot_sim_trace_of(opt.overall_fittest)
  #plt.show()
  ref_data = pd.read_csv(ref_data_path)
  print("RS Model")
  opt_rs.overall_fittest.clean_slow_reactions(0)
  print(opt_rs.overall_fittest)
  rs_data = gen_data_for(opt_rs.overall_fittest, *hyperparams_from_data(opt_rs.overall_fittest.ref_data_path))
  rs_data.to_csv("rs_model_data.csv", index=False)

  # coupled SINDy
  reactions = list(opt._enumerator.generator())
  library = ReactionLibrary.from_ref_data(reactions, ref_data_path)
  res, model = optimize_coupled_sindy(library)
  print("coupled sindy r2/norm:", res.r2_score, "/", res.norm)
  print("cSINDy Model")
  print(model)
  print("Ground Truth")
  print(gt_model)
  print("Intersection")
  print(get_model_intersection(gt_model, model))
  #plot_sim_trace_of(model)
  #plt.show()
  sindy_data = gen_data_for(model, *hyperparams_from_data(model.ref_data_path))
  sindy_data.to_csv("csindy_model_data.csv", index=False)

  # some evaluation...
  print("RS r2/norm:", opt_rs.overall_fittest_value, "/", opt_rs.overall_fittest_norm)
  print("GP r2/norm:", opt.overall_fittest_value, "/", opt.overall_fittest_norm)
  print("SINDy r2/norm:", res.r2_score, "/", res.norm)
  print("GP r2/SINDy r2:", opt.overall_fittest_value / res.r2_score)
  fig, ax = plot_system_comparison(opt.overall_fittest, model, ref_data_path, rs1_label="GP", rs2_label="c-SINDy")
  diff_gp = ((ref_data - gp_data) ** 2).mean().drop("time")
  print(diff_gp, "\nsum", diff_gp.sum())
  diff_sindy = ((ref_data - sindy_data) ** 2).mean().drop("time")
  print(diff_sindy, "\nsum", diff_sindy.sum())
  plt.show()

