from evolib.model_groundtruths import *
from evolib.wnt import wnt_model, fiete_fitted_reactions
from copy import deepcopy

import unittest

class TestModelGroundtruth(unittest.TestCase):

    def test_library_for(self):
      self.assertTrue(True)

    def test_fixate_reactions(self):
      library = fixate_reactions(deepcopy(sir), [1])
      self.assertTrue(len(library.fixated_reactions) == 1 and len(library.non_fixated_reactions) == 1, "Length of new library wrong.")
      self.assertTrue(library.fixated_reactions == [sir.reactions[1]] and library.non_fixated_reactions == [sir.reactions[0]], "Wrong terms fixated.")

      library = fixate_reactions(deepcopy(sir), [1], invert=True)
      self.assertTrue(len(library.fixated_reactions) == 1 and len(library.non_fixated_reactions) == 1, "Length of new library wrong.")
      self.assertTrue(library.fixated_reactions == [sir.reactions[0]] and library.non_fixated_reactions == [sir.reactions[1]], "Wrong terms fixated.")

      library1 = fixate_reactions(deepcopy(wnt_model), fiete_fitted_reactions)
      library2 = fixate_reactions(deepcopy(wnt_model), fiete_fitted_reactions, invert=True)
      self.assertTrue(library1.fixated_reactions == library2.non_fixated_reactions, "Inversion failed.")
      self.assertTrue(library1.non_fixated_reactions == library2.fixated_reactions, "Inverseion failed.")


if __name__ == '__main__':
    unittest.main()
