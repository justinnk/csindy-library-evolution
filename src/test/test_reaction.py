from evolib.model_groundtruths import *
from evolib.reaction import Reaction
from copy import deepcopy

import unittest

class TestModelGroundtruth(unittest.TestCase):

    def test_library_for(self):
      self.assertTrue(True)

    def test_eq(self):
      base = Reaction([0, 1], [1, 1], 1.0, 3, species_names=list("ABC"))
      different_rate = Reaction([0, 1], [1, 1], 0.0, 3, species_names=list("ABC"))
      different_names = Reaction([0, 1], [1, 1], 1.0, 3, species_names=list("SIR"))
      different_reactants = Reaction([0, 1, 1], [1, 1], 1.0, 3, species_names=list("SIR"))
      different_products = Reaction([0, 1], [0, 1], 1.0, 3, species_names=list("SIR"))
      different_num_species = Reaction([0, 1], [0, 1], 1.0, 2, species_names=list("SIR"))

      self.assertFalse(base == different_rate, "Reactions with different rate are not equivalent.")
      self.assertTrue(base == different_names, "Species names are not important.")
      self.assertFalse(base == different_reactants, "Reactants must be equal.")
      self.assertFalse(base == different_products, "Products must be equal.")
      self.assertFalse(base == different_num_species, "Number of species must be equal.")


if __name__ == '__main__':
    unittest.main()
