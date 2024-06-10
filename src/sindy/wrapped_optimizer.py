"""
(unchanged) WrappedOptimizer class from PySINDy package (seems unavailable in packaged
version). Included here only to test coupled SINDy with ensembling defined in PySINDy.
License is as follows:

MIT License

Copyright (c) for portions of project PySINDy are held by Markus Quade, 2019 as
part of project sparsereg. All other copyright for project PySINDy are held by
Brian de Silva and Kathleen Champion 2019.

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


import numpy as np
from sklearn.multioutput import MultiOutputRegressor

from pysindy.optimizers.base import BaseOptimizer

COEF_THRESHOLD = 1e-14


class WrappedOptimizer(BaseOptimizer):
    """Wrapper class for generic optimizers/sparse regression methods

    Enables single target regressors (i.e. those whose predictions are
    1-dimensional) to perform multi target regression (i.e. predictions
    are 2-dimensional).  Also allows unbiasing & normalization for
    optimizers that would otherwise not include it.

    Args:
        optimizer: wrapped optimizer/sparse regression method

    Parameters
    ----------
    optimizer: estimator object
        The optimizer/sparse regressor to be wrapped, implementing ``fit`` and
        ``predict``. ``optimizer`` should also have the attribute ``coef_``.
        Any optimizer that supports a ``fit_intercept`` argument should
        be initialized to False.

    """

    def __init__(self, optimizer, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.optimizer = MultiOutputRegressor(optimizer)

    def _reduce(self, x, y):
        coef_shape = (y.shape[1], x.shape[1])
        self.coef_ = np.zeros(coef_shape)
        self.ind_ = np.ones(coef_shape)
        self.optimizer.fit(x, y)
        coef_list = [
            np.reshape(est.coef_, (-1, coef_shape[1]))
            for est in self.optimizer.estimators_
        ]
        self.coef_ = np.concatenate(coef_list, axis=0)
        self.ind_ = np.abs(self.coef_) > COEF_THRESHOLD
        self.intercept_ = 0.0
        return self

    def predict(self, x):
        return self.optimizer.predict(x)

    @property
    def complexity(self):
        return np.count_nonzero(self.coef_) + np.count_nonzero(self.intercept_)
