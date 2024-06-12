"""
Implementation of coupled SINDy based on the publication:

[1] Burrage, P.M., Weerasinghe, H.N., Burrage, K.: Using a library of chemical re-
actions to fit systems of ordinary differential equations to agent-based models: a
machine learning approach. Numerical Algorithms (1 2024). https://doi.org/10.
1007/s11075-023-01737-0

The underlying Sparse Identification of Non-linear Dynamics (SINDy) was proposed in:

[2] Brunton, S.L., Proctor, J.L., Kutz, J.N.: Discovering governing equations from
data by sparse identification of nonlinear dynamical systems. Proceedings of the
National Academy of Sciences 113, 3932–3937 (4 2016). https://doi.org/10.1073/
pnas.1517384113


License for the below implementation:

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
from evolib.integrate import plot_sim_trace_of, plot_system_comparison
from evolib.sindy.util import prepare_sindy_data
from evolib.model_groundtruths import library_for, sir, predatorprey, wnt

from copy import deepcopy

from scipy.optimize import nnls, lsq_linear
from sklearn.linear_model import Lasso, LinearRegression
from sklearn.metrics import r2_score
import numpy as np
import matplotlib.pyplot as plt
from pysindy.optimizers.stlsq import STLSQ
from pysindy.optimizers.base import EnsembleOptimizer
from evolib.sindy.wrapped_optimizer import WrappedOptimizer
from tqdm import tqdm


class Res:
  """Class to wrap the results from the sparse regression produced by different optimizers."""
  def __init__(self, coef: np.array, data_deriv: np.array, design_deriv: np.array, fixed_coef: np.array = []):
    self.coef = coef
    self.fixed_coef = fixed_coef
    self.r2_score = r2_score(data_deriv, design_deriv @ coef)
    self.norm = np.linalg.norm(data_deriv - design_deriv @ coef)

  def __repr__(self):
    return f"Res({self.coef}, {self.r2_score}, {self.norm})"


class CoupledSINDy:
  """
  Implementation of the coupled SINDy extension of SINDy using the evolib API.
  """
  def __init__(self, library: ReactionLibrary, verbose: bool = False, consider_rate_bounds: bool = False):
    """
    library: list of reactions (ReactionLibrary) to choose from for regression.
    verbose: print additional information.
    consider_rate_bounds: recognize the bounds for the rate specified in each reaction.
    """
    # load file, split concentration data and time
    data, times = prepare_sindy_data(library.ref_data_path)
    self.data = data
    # calculate numerical derivative and put into required matrix form
    self.data_deriv = np.gradient(data, times, axis=0, edge_order=2)
    self.data_deriv = np.transpose(self.data_deriv).reshape(-1)
    self.times = times
    self.library = library
    self.verbose = verbose
    self.consider_rate_bounds = consider_rate_bounds
    if self.verbose:
      self.log("length of library:", len(self.library.reactions), "reactions")

  def log(self, *msg):
    if self.verbose: print("[CoupledSINDy]", *msg)

  def optimize(self, method="nnls"):
    """
    method: The method to use in [nnls, stlsq, lasso, lib_ensemble, time_ensemble]
    Returns: Res(rate constants, r2 score, norm)
    If uncertain, leave the default method "nnls". It was it was found to
    perform far better than all the others tested (the above plus fista and 
    scipy.optimize.least_squares).
    """
    self.log("applying function.")
    progress_bar = tqdm if self.verbose else lambda x: x # only show progress bar if verbose
    theta = np.array([])
    # when there are fixated reactions, subtract their contribution from the design matrix before continuing
    if self.library.fixated_reactions is not None and len(self.library.fixated_reactions) > 0:
      theta_fix = np.transpose(np.stack([r.apply_to(self.data) for r in progress_bar(self.library.fixated_reactions)]))
      theta_fix = theta_fix.reshape(-1, theta_fix.shape[-1])
      coeff_fix = np.array([r.rate for r in self.library.fixated_reactions])
      theta = np.transpose(np.stack([r.apply_to(self.data) for r in progress_bar(self.library.non_fixated_reactions)]))
      self.data_deriv = self.data_deriv - theta_fix @ coeff_fix
    # otherwise use the "standard" design matrix
    else:
      theta = np.transpose(np.stack([r.apply_to(self.data) for r in progress_bar(self.library.reactions)])) # library of reactions
    theta = theta.reshape(-1, theta.shape[-1])
    # collect bounds on the rate constants (default 0, inf) (if enabled)
    # this also forces the method to be constrained_lsq
    bounds = (0, np.inf)
    if self.consider_rate_bounds:
      method = "constrained_lsq"
      lb = np.zeros(len(self.library.non_fixated_reactions))
      ub = np.ones(len(self.library.non_fixated_reactions)) * float("inf")
      for idx, reac in enumerate(self.library.non_fixated_reactions):
        lb[idx] = reac.rate_bounds[0]
        ub[idx] = reac.rate_bounds[1]
      bounds = (lb, ub)
      self.log("rate constants are bounded to:", bounds)
    self.log("starting regression.")
    self.log("shape of design matrix (theta):", theta.shape)
    self.log("shape of data matrix:", self.data_deriv.shape)
    # apply a solver to the constructed problem
    # generally, nnls (the default) works far better than all others
    if method == "nnls": # non-negative least squares
      res = nnls(theta, self.data_deriv, maxiter=1e8)
      res = Res(res[0], self.data_deriv, theta)
      #res = LinearRegression(positive=True, fit_intercept=False).fit(theta, self.data_deriv)
      #res = Res(res.coef_, self.data_deriv, theta)
    elif method == "stlsq": # sequential thresholded least squares (from (py-)SINDy)
      res = STLSQ(threshold=1e-4, max_iter=1000).fit(theta, self.data_deriv)
      res = Res(res.coef_[0], self.data_deriv, theta)
    elif method == "lasso": # the LASSO
      lasso = Lasso(alpha=10.3, max_iter=1000000, tol=1e-12, fit_intercept=False)
      res = lasso.fit(theta, self.data_deriv)
      res = Res(res.coef_, self.data_deriv, theta)
    elif method == "constrained_lsq": # constrained least squares
      res = lsq_linear(theta, self.data_deriv, bounds=bounds)
      res = Res(res.x, self.data_deriv, theta)
    elif method == "lib_ensemble": # library ensembling from PySINDy
      res = EnsembleOptimizer(WrappedOptimizer(LinearRegression(positive=True, fit_intercept=False)), library_ensemble=True, n_models=100, n_candidates_to_drop=2).fit(theta, self.data_deriv)
      res = Res(res.coef_[0], self.data_deriv, theta)
    elif method == "time_ensemble": # time ensembling from PySINDy
      res = EnsembleOptimizer(WrappedOptimizer(LinearRegression(positive=True, fit_intercept=False)), bagging=True, n_models=100, n_subset=10).fit(theta, self.data_deriv)
      res = Res(res.coef_[0], self.data_deriv, theta)
    else:
      raise Exception(f"Method {method} not supported. Choose one of nnls, stlsq, lasso, constrained_lsq, lib_ensemble, time_ensemble!")
    self.log("done.\n", res)
    return res

  def get_sparse_reaction_library(self, res: Res, thresh=None):
    """Get the model determined via regression."""
    library = deepcopy(self.library)
    reactions = []
    if thresh is not None:
      res.coef[res.coef < thresh] = 0.0
    active_reaction_idxs = [idx for idx, x in enumerate(res.coef) if x != 0.0]
    for active_reaction_idx in active_reaction_idxs:
      library.non_fixated_reactions[active_reaction_idx].rate = res.coef[active_reaction_idx]
      reactions += [library.non_fixated_reactions[active_reaction_idx]]
    library.non_fixated_reactions = reactions + library.fixated_reactions
    library.fixated_reactions = []
    return library


def optimize_coupled_sindy(library, method="nnls", thresh=None, verbose=False):
  """Convenience function to apply coupled sindy to library and return resulting model."""
  csindy = CoupledSINDy(library, verbose=verbose)
  res = csindy.optimize(method=method)
  return res, csindy.get_sparse_reaction_library(res, thresh=thresh)


if __name__ == "__main__":
  # for quick testing (run from inside the src/ folder so the data path is correct!)
  from evolib.reaction_enumerator import ReactionEnumerator
  from pprint import pprint 

  gt_model = sir # <--- change this...
  ref_data_path = "data/sir_100.csv" # <--- ...and this
  library = library_for(gt_model, ref_data_path, 2, 3, verbose=True)

  method = "nnls"
  res, model = optimize_coupled_sindy(library, method, verbose=True)
  print("coeff. of determination:", res.r2_score)
  print("residuals:", res.norm)
  print(model)
  plot_sim_trace_of(model, n_points=500)
  plot_system_comparison(model, gt_model, ref_data_path)
  plt.show()

