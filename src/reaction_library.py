"""
Class to represent a library of biochemical reactions and some helper methods.

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


from evolib.reaction import Reaction, full_hist
from evolib.reaction_enumerator import ReactionEnumerator

from random import randint, choice
from typing import List, Union
from copy import deepcopy
from math import gcd
from collections import defaultdict
from ast import literal_eval

import numpy as np


def species_from_csv(path: str):
  """Return (number of species, list of names) from the csv file at path."""
  with open(path, "r") as file:
    header = file.readline().rstrip()
    species_names = header.split(",")[1:]
  return len(species_names), species_names

def get_model_intersection(model1, model2):
  """Return a ReactionLibrary with all common reactions of model1 and model2, disregarding rates and fixated reactions."""
  inter_model = deepcopy(model1)
  inter_model.fixated_reactions = []
  inter_model.non_fixated_reactions = []
  for reac1 in model1.non_fixated_reactions:
    for reac2 in model2.reactions:
      if reac1.equal_structure(reac2):
        inter_model.non_fixated_reactions.append(deepcopy(reac1))
  return inter_model


class ReactionLibrary:
  """A library of biochemical reactions, some of which can be marked as immutable during learning."""
  def __init__(self, non_fixated_reactions: List[Reaction], num_species, *, ref_data_path: str = "", species_names=[], enumerator=None, fixated_reactions=[]):
    self.fixated_reactions = fixated_reactions
    self.non_fixated_reactions = non_fixated_reactions
    self.ref_data_path = ref_data_path
    self.enumerator = enumerator
    self.num_species = num_species
    self.species_names = species_names
    for reaction in self.reactions:
      reaction.species_names = species_names
    self.clean_duplicate_reactions()

  @classmethod
  def from_ref_data(cls, reactions: List[Reaction], ref_data_path: str, enumerator=None):
    """Initializes the library for re-discovering a specific reference dataset."""
    if isinstance(ref_data_path, str) and ref_data_path != "":
      num_species, species_names = species_from_csv(ref_data_path)
      return cls(reactions, num_species, ref_data_path=ref_data_path, species_names=species_names, enumerator=enumerator)
    else:
      print("Error: ref_data_path should be a non-empty string!")
      exit(-1)

  @property
  def reactions(self):
    """All reactions within the library (fixated or not)."""
    if self.fixated_reactions is None:
      return self.non_fixated_reactions
    return self.fixated_reactions + self.non_fixated_reactions

  def __str__(self) -> str:
    out = ""
    if self.fixated_reactions is not None:
      for reaction in self.fixated_reactions:
        out += str(reaction) + " (fix)\n"
    for reaction in self.non_fixated_reactions:
      out += str(reaction) + "\n"
    return out

  def __repr__(self) -> str:
    return f"ReactionLibrary({self.non_fixated_reactions}, {self.num_species}, ref_data_path={self.ref_data_path.__repr__()}, species_names={self.species_names}, fixated_reactions={self.fixated_reactions})"

  def to_latex(self, gt_model=None) -> str:
    """Returns a LaTeX math formulation of the library as string."""
    out = "$\\begin{aligned}"
    for reaction in self.fixated_reactions:
      out += reaction.to_latex(gt_model) + r"\\[-0.2cm]"
    if len(self.fixated_reactions) > 0:
      out += r"\cline{1-3}"
    for reaction in self.non_fixated_reactions:
      out += reaction.to_latex(gt_model) + r"\\[-0.2cm]"
    out += "\\end{aligned}$"
    return out

  def print_model(self) -> str:
    """Same as __str__ but skipping 0-rate reactions."""
    out = ""
    if self.fixated_reactions is not None:
      for reaction in self.fixated_reactions:
        if reaction.rate > 0.0:
          out += str(reaction) + " (fix)\n"
    for reaction in self.non_fixated_reactions:
      if reaction.rate > 0.0:
        out += str(reaction) + "\n"
    return out

  def clean_duplicate_reactions(self):
    """Clean duplicate reactions and fill up empty slots again."""
    len_before = len(self.non_fixated_reactions)
    # remove duplicates
    self.non_fixated_reactions = list(set(self.non_fixated_reactions))
    # fill back up (if there is an enumerator)
    if self.enumerator is not None and len(self.non_fixated_reactions) < len_before:
      while len(self.non_fixated_reactions) < len_before:
        self.non_fixated_reactions.append(self.enumerator.get_random())
      # recurse in case we added more reactions which could be duplicates
      self.clean_duplicate_reactions()

  def clean_slow_reactions(self, thresh: float):
    """Remove all reactions with a rate below thres (inline/destructive)."""
    non_fixated_reactions = []
    for r in self.non_fixated_reactions:
      if r.rate > thresh:
        non_fixated_reactions.append(r)
    self.non_fixated_reactions = non_fixated_reactions

  def crossover(self, other: "ReactionLibrary") -> None:
    """Swap one reaction, chosen uniformly at random, between this and other."""
    # determine places for crossover
    place_self = randint(0, len(self.non_fixated_reactions) - 1)
    place_other = randint(0, len(other.non_fixated_reactions) - 1)
    # swap reactions
    temp = other.non_fixated_reactions[place_other]
    other.non_fixated_reactions[place_other] = self.non_fixated_reactions[place_self]
    self.non_fixated_reactions[place_self] = temp
    # remove possible duplicates
    self.clean_duplicate_reactions()
    other.clean_duplicate_reactions()

  def mutate(self, constraints = []):
    """Swap one random reaction chosen uniformly at ranom with a new random one."""
    # determine place for mutation
    place = randint(0, len(self.non_fixated_reactions) - 1)
    #self.reactions[place].mutate()
    self.non_fixated_reactions[place] = self.enumerator.get_random()
    self.clean_duplicate_reactions()

  def apply_odes(self, t, state, pbar=None, pbar_state=None):
    """Return the gradient in the phase space at the point 'state'."""
    if pbar is not None and pbar_state is not None:
      last_t, dt = pbar_state
      n = (t - last_t)/dt
      pbar.update(n)
      pbar_state[0] = last_t + dt * n

    deriv = np.zeros_like(state)
    for reaction in self.reactions:
      deriv += reaction.get_deriv_terms(state) * reaction.rate
    return deriv


