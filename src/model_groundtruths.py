"""
Contains definitions for the benchmarked models using the evolib API 
and some helper functions for fixating reactions while learning and
constructing a complete library.

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
from evolib.reaction_enumerator import ReactionEnumerator
from evolib.wnt import wnt_model

from typing import Union, List


def library_for(model: ReactionLibrary, ref_data_path: str, max_num_left: int, max_num_right: int, *, shuffle: bool=False, avoid_colinear: bool=False, verbose: bool=False) -> ReactionLibrary:
  """Generate a complete library of possible reactions for the given model assuming the contained number of species."""
  reac_enum = ReactionEnumerator(num_species=model.num_species, max_num_left=max_num_left, max_num_right=max_num_right, max_stoichiometry=1, constraints=[], shuffle=False, species_names=model.species_names)
  if verbose:
    print("[library_for]", "Library generation.")
  reactions = list(reac_enum.generator())
  if verbose:
    print("[library_for]", "Length of library:", len(reactions), "(excl. colinear), ", reac_enum.get_number(), "(incl. colinear)")
  library = ReactionLibrary(reactions, model.num_species, ref_data_path=ref_data_path, species_names=reac_enum.species_names, fixated_reactions=[])
  if verbose:
    print("[library_for]", "Final length of library:", len(library.reactions), "(", len(library.fixated_reactions), "of them fixated)")
  return library


def fixate_reactions(library: ReactionLibrary, reactions: Union[List[Reaction],List[int]], invert: bool = False) -> ReactionLibrary:
  """Fixate list of reactions (based on reaction or index in list). If invert=True, fixate all but those reactions."""
  if len(reactions) == 0:
    return library
  non_fixated_reactions = []
  for idx, reac in enumerate(library.non_fixated_reactions):
    temp = ((isinstance(reactions[0], int) and idx in reactions) or (isinstance(reactions[0], Reaction) and reac in reactions))
    if (not invert and temp) or (invert and not temp):
      library.fixated_reactions.append(reac)
    else:
      non_fixated_reactions.append(reac)
  library.non_fixated_reactions = non_fixated_reactions
  return library

def get_latex(model: ReactionLibrary, model_name: str):
  """Output a reaction system as pdf."""
  import matplotlib.pyplot as plt
  import scienceplots
  plt.style.use(["science", "ieee"])
  fig, ax = plt.subplots(figsize=(1, 1))
  ax.set_axis_off()
  ax.text(0, 0, model.to_latex())
  return fix, ax #fig.savefig(f"{name}_model_gt.pdf")

# compartmental susceptible-infected-recovered model of disease dynamics
sir = ReactionLibrary([
  Reaction([0, 1], [1, 1], 0.02, 3),
  Reaction([1], [2], 5.0, 3)
], 3, species_names=["S", "I", "R"])

# lotka-volterra dynamics of predator and prey relation
predatorprey = ReactionLibrary([
  Reaction([0, 1], [1, 1], 0.01, 2),
  Reaction([1], [], 8.0, 2),
  Reaction([0], [0, 0], 10.0, 2)
], 2, species_names=["S", "W"])

# wnt model from Staehlke 2020 (model defined in wnt_model.py)
wnt = wnt_model



