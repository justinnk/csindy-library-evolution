"""
A wrapper for the classical SINDy implementation.
Not used in the paper, included here only for the curious :)
Does not produce reaction-style outputs and is just used
to test fitting the data in a typical way. It seems that the
results are also generally not interpretable as reactions 
(the same functional terms appear with different coefficients
in the resulting reactions).
"""

from evolib.sindy.util import prepare_sindy_data

from pysindy import SINDy
from pysindy.feature_library import PolynomialLibrary
from pysindy.optimizers import EnsembleOptimizer, STLSQ
import pandas as pd
from sklearn.linear_model import Lasso, LinearRegression, LassoCV

_str_to_opt = {
    "stlsq": STLSQ,
    "lasso": Lasso,
    "lsq": LinearRegression,
}

class MySINDy:
  def __init__(self, ref_data_path: str, optimizer: str, **opt_kwargs):
    self.library = PolynomialLibrary(degree=2)
    self.sindy = SINDy(optimizer=_str_to_opt[optimizer](**opt_kwargs, fit_intercept=False), feature_library=self.library)
    self.data, self.times = prepare_sindy_data(ref_data_path)

  def optimize(self):
    self.sindy.fit(self.data, t=self.times)

  def score(self):
    return self.sindy.score(self.data, t=self.times)

class EnsembleSINDy:
  def __init__(self, ref_data_path: str, optimizer: str, **opt_kwargs):
    self.library = PolynomialLibrary(degree=2)
    self.sindy = SINDy(optimizer=_str_to_opt[optimizer](**opt_kwargs, fit_intercept=False), feature_library=self.library)
    self.data, self.times = prepare_sindy_data(ref_data_path)

  def optimize(self, **kw_args):
    self.sindy.fit(self.data, t=self.times, library_ensemble=True, **kw_args)

  def score(self):
    return self.sindy.score(self.data, t=self.times)


if __name__ == "__main__":
  sindy = MySINDy("data/sir_100.csv", "lasso")
  #sindy = EnsembleSINDy("../data/wnt_5000.csv", "stlsq", alpha=0.1, threshold=1e-5)
  sindy.optimize()
  #sindy.optimize(n_models=100, n_candidates_to_drop=10)
  print(sindy.score())
  print(sindy.sindy.equations(precision=4))
  sindy.sindy.print(["S", "I", "R"])

