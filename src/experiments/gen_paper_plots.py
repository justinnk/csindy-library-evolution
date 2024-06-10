"""
Script to generate the plots from the experiment data as presented in the paper.

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


if __name__ == "__main__":
  from evolib.experiments.experiment_definitions import *
  from evolib.experiments.model_generation_experiment import EXPERIMENTS_BASEPATH
  from glob import glob
  import os
  import sys
  import matplotlib.pyplot as plt
  import scienceplots
  plt.style.use(["science", "ieee"])

  FIGURES = "figures"
  DATES = dict(
    sir=dict(
      evolib="2024-04-28T20_57",
      csindy="2024-04-28T20_57",
      rs="2024-04-28T20_57",
    ),
    predatorprey=dict(
      evolib="2024-05-05T09_47",#"2024-04-28T20_50",
      csindy="2024-04-28T20_50",
      rs="2024-04-28T20_50",
    ),
    wnt=dict(
      evolib="2024-04-28T21_02",
      csindy="2024-05-02T10_52",
      rs="2024-04-30T10_47",
    ),
    wnt_ext=dict(
      evolib="2024-05-02T23_49",#"2024-04-28T20_30",
      csindy="2024-04-30T17_45",
      rs="2024-05-05T00_29"#"2024-04-29T18_49",
    ),
    wnt_constrained=dict(
      evolib="2024-05-01T10_06",
    ),
    wnt_ext_constrained=dict(
      evolib="2024-05-02T15_16",#"2024-05-01T13_55",
    )
  )

  def evaluate(what: str, path: str, fig=None, ax=None):
    print(path)
    sys.path.append(path)
    from hyperparameters import exp
    #try:
    fig, ax = getattr(exp, f"evaluate_{what}")(path, fig, ax)
    #except Exception as ex:
    #  print(ex)
    #  exit(-1)
    del sys.modules["hyperparameters"]
    sys.path.remove(path)
    return fig, ax

  # model

  for model in DATES:
    for method in DATES[model]:
      fig, ax = evaluate("model", os.path.join(EXPERIMENTS_BASEPATH, method, model, DATES[model][method]))
      fig.tight_layout()
      fig.savefig(f"{FIGURES}/{model}_model_{method}.pdf")

  # fit

  fig, ax = evaluate("fit", os.path.join(EXPERIMENTS_BASEPATH, "evolib", "sir", DATES["sir"]["evolib"]))
  fig, ax = evaluate("fit", os.path.join(EXPERIMENTS_BASEPATH, "csindy", "sir", DATES["sir"]["csindy"]), fig, ax)
  ax.set_xlabel("time [units]")
  ax.set_ylabel("amount")
  ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  #ax.set_title("SIR Model")
  fig.set_size_inches(3, 1.8)
  fig.tight_layout()
  fig.subplots_adjust(right=1.0)
  fig.savefig(f"{FIGURES}/sir_fit.pdf")

  fig, ax = evaluate("fit", os.path.join(EXPERIMENTS_BASEPATH, "evolib", "predatorprey", DATES["predatorprey"]["evolib"]))
  fig, ax = evaluate("fit", os.path.join(EXPERIMENTS_BASEPATH, "csindy", "predatorprey", DATES["predatorprey"]["csindy"]), fig, ax)
  ax.set_xlabel("time [units]")
  ax.set_ylabel("amount")
  ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  #ax.set_title("Predator-Prey Model")
  fig.set_size_inches(3, 1.8)
  fig.tight_layout()
  fig.subplots_adjust(right=1.0)
  fig.savefig(f"{FIGURES}/predatorprey_fit.pdf")

  
  fig, ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True, figsize=(6,4))

  _ax = ax[0,0]
  fig, _ax = evaluate("fit", os.path.join(EXPERIMENTS_BASEPATH, "evolib", "wnt", DATES["wnt"]["evolib"]), fig, _ax)
  #fig, ax = evaluate("fit", os.path.join(EXPERIMENTS_BASEPATH, "csindy", "wnt", DATES["wnt"]["csindy"]), fig, ax)
  _ax.set_xlabel("time [min]")
  _ax.set_ylabel("amount")
  _ax.set_yscale("log")
  _ax.set_ylim((1e-6, 2e5))
  _ax.get_legend().remove()
  _ax.set_title("From Scratch")
  #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  #fig.set_size_inches(3, 1.8)
  #fig.tight_layout()
  #fig.savefig(f"{FIGURES}/wnt_fit.pdf")

  _ax = ax[1,0]
  fig, _ax = evaluate("fit", os.path.join(EXPERIMENTS_BASEPATH, "evolib", "wnt_ext", DATES["wnt_ext"]["evolib"]), fig, _ax)
  #fig, ax = evaluate("fit", os.path.join(EXPERIMENTS_BASEPATH, "csindy", "wnt_ext", DATES["wnt_ext"]["csindy"]), fig, _ax)
  _ax.set_xlabel("time [min]")
  _ax.set_ylabel("amount")
  _ax.set_yscale("log")
  _ax.set_ylim((1e-6, 2e5))
  _ax.set_title("Extension")
  _ax.get_legend().remove()
  #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  #fig.set_size_inches(3, 1.8)
  #fig.tight_layout()
  #fig.subplots_adjust(right=1.0)
  #fig.savefig(f"{FIGURES}/wnt_ext_fit.pdf")

  _ax = ax[0,1]
  fig, _ax = evaluate("fit", os.path.join(EXPERIMENTS_BASEPATH, "evolib", "wnt_constrained", DATES["wnt_constrained"]["evolib"]), fig, _ax)
  _ax.set_xlabel("time [min]")
  _ax.set_ylabel("amount")
  _ax.set_yscale("log")
  _ax.set_ylim((1e-6, 2e5))
  _ax.get_legend().remove()
  _ax.set_title("From Scratch (Constrained)")
  #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  #fig.set_size_inches(3, 1.8)
  #fig.tight_layout()
  #fig.subplots_adjust(right=1.0)
  #fig.savefig(f"{FIGURES}/wnt_constrained_fit.pdf")

  _ax = ax[1,1]
  fig, _ax = evaluate("fit", os.path.join(EXPERIMENTS_BASEPATH, "evolib", "wnt_ext_constrained", DATES["wnt_ext_constrained"]["evolib"]), fig, _ax)
  _ax.set_xlabel("time [min]")
  _ax.set_ylabel("amount")
  _ax.set_yscale("log")
  _ax.set_ylim((1e-6, 2e5))
  _ax.get_legend().remove()
  _ax.set_title("Extension (Constrained)")
  #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  #fig.set_size_inches(3, 1.8)
  #fig.tight_layout()
  #fig.subplots_adjust(right=1.0)
  #fig.savefig(f"{FIGURES}/wnt_ext_constrained_fit.pdf")

  #fig.legend(*ax[0,0].get_legend_handles_labels(), loc="upper center", bbox_to_anchor=(0.5, 0.03), ncol=3)
  fig.tight_layout()
  fig.savefig(f"{FIGURES}/wnt_all_fit.pdf")

  # convergence

  fig, ax = evaluate("convergence", os.path.join(EXPERIMENTS_BASEPATH, "evolib", "sir", DATES["sir"]["evolib"]))
  fig, ax = evaluate("convergence", os.path.join(EXPERIMENTS_BASEPATH, "csindy", "sir", DATES["sir"]["csindy"]), fig, ax)
  fig, ax = evaluate("convergence", os.path.join(EXPERIMENTS_BASEPATH, "rs", "sir", DATES["sir"]["rs"]), fig, ax)
  ax.set_yscale("log")
  ax.set_xlim((-0.5, 60))
  ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
  ax.set_xlabel("Steps")
  ax.set_ylabel("$\ell^2$ Norm")
  ax.set_title("SIR Model")
  fig.tight_layout()
  fig.subplots_adjust(bottom=0.2)
  fig.savefig(f"{FIGURES}/sir_convergence.pdf")

  fig, ax = evaluate("convergence", os.path.join(EXPERIMENTS_BASEPATH, "evolib", "predatorprey", DATES["predatorprey"]["evolib"]))
  fig, ax = evaluate("convergence", os.path.join(EXPERIMENTS_BASEPATH, "csindy", "predatorprey", DATES["predatorprey"]["csindy"]), fig, ax)
  fig, ax = evaluate("convergence", os.path.join(EXPERIMENTS_BASEPATH, "rs", "predatorprey", DATES["predatorprey"]["rs"]), fig, ax)
  ax.set_yscale("log")
  ax.set_xlim((-0.5, 60))
  ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
  ax.set_xlabel("Steps")
  ax.set_ylabel("$\ell^2$ Norm")
  ax.set_title("Predator-Prey Model")
  fig.tight_layout()
  fig.subplots_adjust(bottom=0.2)
  fig.savefig(f"{FIGURES}/predatorprey_convergence.pdf")


  fig, ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True, figsize=(5,3))

  _ax = ax[0,0]
  fig, _ax = evaluate("convergence", os.path.join(EXPERIMENTS_BASEPATH, "evolib", "wnt", DATES["wnt"]["evolib"]), fig, _ax)
  fig, _ax = evaluate("convergence", os.path.join(EXPERIMENTS_BASEPATH, "csindy", "wnt", DATES["wnt"]["csindy"]), fig, _ax)
  fig, _ax = evaluate("convergence", os.path.join(EXPERIMENTS_BASEPATH, "rs", "wnt", DATES["wnt"]["rs"]), fig, _ax)
  _ax.set_yscale("log")
  #_ax.legend()
  _ax.set_xlabel("Steps")
  _ax.set_ylabel("$\ell^2$ Norm")
  _ax.set_title("From Scratch")
  _ax.grid(True)
  #fig.tight_layout()
  #fig.savefig(f"{FIGURES}/wnt_convergence.pdf")

  _ax = ax[1,0]
  fig, _ax = evaluate("convergence", os.path.join(EXPERIMENTS_BASEPATH, "evolib", "wnt_ext", DATES["wnt_ext"]["evolib"]), fig, _ax)
  fig, _ax = evaluate("convergence", os.path.join(EXPERIMENTS_BASEPATH, "csindy", "wnt_ext", DATES["wnt_ext"]["csindy"]), fig, _ax)
  fig, _ax = evaluate("convergence", os.path.join(EXPERIMENTS_BASEPATH, "rs", "wnt_ext", DATES["wnt_ext"]["rs"]), fig, _ax)
  _ax.set_yscale("log")
  #_ax.legend()
  _ax.set_xlabel("Steps")
  _ax.set_ylabel("$\ell^2$ Norm")
  _ax.set_title("Extension")
  _ax.grid(True)
  #fig.tight_layout()
  #fig.savefig(f"{FIGURES}/wnt_ext_convergence.pdf")

  _ax = ax[0,1]
  fig, _ax = evaluate("convergence", os.path.join(EXPERIMENTS_BASEPATH, "evolib", "wnt_constrained", DATES["wnt_constrained"]["evolib"]), fig, _ax)
  _ax.set_yscale("log")
  #_ax.legend()
  _ax.set_xlabel("Steps")
  #_ax.set_ylabel("$\ell^2$ Norm")
  _ax.set_title("From Scratch (Constrained)")
  _ax.grid(True)
  #fig.tight_layout()
  #fig.savefig(f"{FIGURES}/wnt_constrained_convergence.pdf")

  _ax = ax[1,1]
  fig, _ax = evaluate("convergence", os.path.join(EXPERIMENTS_BASEPATH, "evolib", "wnt_ext_constrained", DATES["wnt_ext_constrained"]["evolib"]), fig, _ax)
  _ax.set_yscale("log")
  #_ax.legend()
  _ax.set_xlabel("Steps")
  #_ax.set_ylabel("$\ell^2$ Norm")
  _ax.set_title("Extension (Constrained)")
  _ax.grid(True)
  #fig.tight_layout()
  #fig.savefig(f"{FIGURES}/wnt_ext_constrained_convergence.pdf")

  #fig.subplots_adjust(bottom=0.3, wspace=0.33)
  #ax[1,1].legend(*ax[0,0].get_legend_handles_labels(), loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=False, shadow=False, ncol=3)
  fig.legend(*ax[0,0].get_legend_handles_labels(), loc="upper center", bbox_to_anchor=(0.5, 0.03), ncol=3)
  fig.tight_layout()
  fig.savefig(f"{FIGURES}/wnt_all_convergence.pdf")


