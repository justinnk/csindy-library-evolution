"""
Class to enumerate all possible reactions within some
constraints on, e.g., the number of reactants.

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

from random import shuffle, choice
from typing import List, Generator
from itertools import product, combinations_with_replacement 
from math import gcd

from evolib.reaction import Reaction


class ReactionEnumerator:
  """
  Enumerates possible reactions based on number of species and maximum allowed number of reactants and products.
  Reactions that have no effect, e.g. " -> " or "A -> A", are omitted.
  Uses caching to avoid repeated enumeration, so reusing a single instance is desirable.
  You will probably want to use either the `generator()` method, which returns a Python generator or
  the `get_random()` method, which returns a single randomly chosen reaction.
  """
  def __init__(
      self,
      num_species: int,
      max_num_left: int = 2,
      max_num_right: int = 2,
      max_stoichiometry: int = 99,
      constraints: List = [],
      shuffle: bool = True,
      species_names: List[str] = [],
      reaction_blacklist: List[Reaction] = [],
      subgroups: List = []
    ): 
    self.num_species = num_species
    self.max_num_left = max_num_left
    self.max_num_right = max_num_right
    self.max_stoichiometry = max_stoichiometry
    self.constraints = constraints
    self.shuffle = shuffle
    self.species_names = species_names
    self.reaction_blacklist = [(r.reactands, r.products) for r in reaction_blacklist]
    self.subgroups = subgroups
    # internal
    self._cache = {}


  def generator(self) -> Generator[Reaction, None, None]:
    """Returns a generator object that delivers a (possibly shuffled) list of all possible reactions."""
    # views a reaction c0 S0 + c1 S1 + ... -> d0 S0 + d1 S1 + ... @ r0 as a list [c0, c1, ..., d0, d1, ...] of coefficients
    left = self._get_one_side(self.max_num_left)
    right = self._get_one_side(self.max_num_right)
    # the possible rules are all combinations of the above...
    res = map(lambda x: (list(x[0]), list(x[1])), product(left, right))
    # ...minus the empty reaction and reactions that have no effect on the state
    res = filter(lambda x: self._filter_nochange(*x), res)
    res = filter(lambda x: self._filter_blacklisted(*x), res)
    res = filter(lambda x: self._filter_stoichiometry(*x), res)
    if len(self.subgroups) > 0:
      res = filter(lambda x: self._filter_subgroups(*x), res)
    if self.shuffle:
      res = list(res)
      shuffle(res)
    # for the final result, we map the stoichiometric coefficients c and d to reaction instances
    return map(lambda x: self._to_reaction(*x), res)

  def get_random(self) -> Reaction:
    """Returns a random reaction from the set of all possible reactions."""
    left = self._get_one_side(self.max_num_left)
    right = self._get_one_side(self.max_num_right)
    l_chosen = choice(left)
    r_chosen = choice(right)
    while not self._filter_nochange(l_chosen, r_chosen) or\
          not self._filter_blacklisted(l_chosen, r_chosen) or\
          not self._filter_stoichiometry(l_chosen, r_chosen) or\
          (len(self.subgroups) > 0 and not self._filter_subgroups(l_chosen, r_chosen)):
      l_chosen = choice(left)
      r_chosen = choice(right)
    # for the final result, we map the stoichiometric coefficients c and d to reaction instances
    return self._to_reaction(l_chosen, r_chosen)

  def get_number(self, consider_blacklisted: bool = False) -> int:
    """Returns the number of all possible (sensible) reactions."""
    left = len(self._get_one_side(self.max_num_left))
    right = len(self._get_one_side(self.max_num_right))
    if not consider_blacklisted:
      return left * right - left
    return left * right - left - len(self.reaction_blacklist)

  def get_lhs_number(self) -> int:
    """Returns the number of all possible combinations of 0 to max_num_left species."""
    left = len(self._get_one_side(self.max_num_left))
    return left

  def get_rhs_number(self) -> int:
    """Returns the number of all possible combinations of 0 to max_num_right species."""
    right = len(self._get_one_side(self.max_num_right))
    return right

  def _get_one_side(self, max_order: int = 2) -> Generator[List, None, None]:
    # (views a reaction c0 S0 + c1 S1 + ... -> d0 S0 + d1 S1 + ... @ r0 as a list [c0, c1, ..., d0, d1, ...] of coefficients)
    # to get all n-ary interactions on one side of the reaction (c's or d's), we find all strings of c's or d's whose sum is <= max_order
    # if there is no list cached for the required max order, recalculate
    if max_order not in self._cache:
      self._cache[max_order] = []
      # enumerate all combinations of 1, 2, ..., max_order species id's
      for n in range(1, max_order + 1):
        self._cache[max_order].extend(list(combinations_with_replacement(range(self.num_species), n)))
      # don't forget empty left/right side
      self._cache[max_order].append(())
    # return the cached list of partial reactions
    return self._cache[max_order]

  def _filter_subgroups(self, left, right) -> bool:
    for idx1 in range(len(left)):
      for idx2 in range(idx1 + 1, len(left)):
        for sg in self.subgroups:
          if (left[idx1] in sg and left[idx2] not in sg) or\
             (left[idx1] not in sg and left[idx2] in sg):
            return False
    for idx1 in range(len(right)):
      for idx2 in range(idx1 + 1, len(right)):
        for sg in self.subgroups:
          if (right[idx1] in sg and right[idx2] not in sg) or\
             (right[idx1] not in sg and right[idx2] in sg):
            return False
    return True

  def _filter_stoichiometry(self, left, right) -> bool:
    hist_l = {species: len([1 for _s in left if _s == species]) for species in range(self.num_species)}
    hist_r = {species: len([1 for _s in right if _s == species]) for species in range(self.num_species)}
    changes = [hist_r[r] - hist_l[l] for l, r in zip(hist_l, hist_r)]
    #print(left, right, changes, any([abs(c) > self.max_stoichiometry for c in changes]))
    return not any([abs(c) > self.max_stoichiometry for c in changes])

  def _filter_blacklisted(self, left, right):
    return (list(left), list(right)) not in self.reaction_blacklist

  def _filter_nochange(self, left: list, right: list):
    return left != right # Note: this depends on the order in left and right being the same, could add sorted to be safe.

  def _to_reaction(self, left: list, right: list):
    return Reaction(list(left), list(right), 0.0, self.num_species, species_names=self.species_names)

