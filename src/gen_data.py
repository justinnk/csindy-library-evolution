"""
Script to generates synthetic datasets from the model descriptions.

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

from evolib.integrate import gen_data_for
from evolib.model_groundtruths import sir, predatorprey, wnt
from evolib.wnt import wnt_init

import matplotlib.pyplot as plt


def plot(data, title=""):
  data.plot(x="time", marker="o", lw=0, title=title)


if __name__ == "__main__":
  for n_points in [0, 10, 100, 1000, 5000]:
    sir_data = gen_data_for(sir, [1980.0, 20.0, 0.0], n_points=n_points)
    plot(sir_data, title=f"SIR {n_points}")
    sir_data.to_csv(f"data/sir_{n_points}.csv", index=False)

    predatorprey_data = gen_data_for(predatorprey, [1000.0, 20.0], t_end=2.0, n_points=n_points)
    plot(predatorprey_data, title=f"predator prey {n_points}")
    predatorprey_data.to_csv(f"data/predatorprey_{n_points}.csv", index=False)

    wnt_data = gen_data_for(wnt, wnt_init, t_end=1440, n_points=n_points)
    wnt_data = wnt_data[(wnt_data.time > 15) & (wnt_data.time < 500)]
    plot(wnt_data, title=f"wnt pathway {n_points}")
    wnt_data.to_csv(f"data/wnt_{n_points}.csv", index=False)

    #plt.show()


