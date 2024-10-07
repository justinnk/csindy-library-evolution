[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11654440.svg)](https://doi.org/10.5281/zenodo.11654440)
# Evolving Libraries (evolib)
 
Code artifacts for the paper "Discovering Biochemical Reaction Models by Evolving Libraries" to appear at the CMSB'24.

## :cd: Setup

### Required Software

- Linux operating system; tested on Fedora 38 Workstation and Ubuntu 22.04.4 LTS.
  - Should generally also work on Windows but some adaptations may be required
- Python; tested on version 3.11.9 (Fedora) and 3.10.12 (Ubuntu)
- The dependencies listed in `requirements.txt`
  - can be installed automatically, see the next section

### Installation Guide
(for Linux)

1. Clone the repository:
```shell
git clone <url>
```
2. Create a virtual environment and install the dependencies
```shell
# depending on your linux distribution, you
# may use either "python3" or "python"
python3 -m venv .venv                                    
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```
3. (locally) install the evolib python package in editable mode
```shell
cd src
pip install -e .
```
(if this fails, you may need to start over from step 1, but this time update `setuptools` in step 2 with `pip install --upgrade setuptools`)

You should now have a virtual environment stored in the folder `.venv` with all the dependencies and the module `evolib` (this repository) installed.

## :bar_chart: Reproduce Results

These are the steps to reproduce the figures from the paper.
The steps assume that you are inside the virtual environment created in the installation guide above.
Running all experiments may take around 35h, depending on your hardware.

0. Navigate to the source directory: `cd src`
1. Generate the datasets: `python experiments/gen_data.py`.
2. Run the experiments: `python experiments/run_experiments.py`. This will run the experiments one after another and store the data in `src/experiment_results/<method>/<model>/<timestamp>/`. The default settings assume a machine with at least 10 cpu cores. This can be adjusted for each experiment by changing the `parallel_macroreps` parameter.
3. Open the file `gen_paper_plots.py` and fill in the correct folder names (timestamps) at the top. You can use `ls experiment_results/*/*` in the `src` directory to get a good overview.
4. Plot the results: `python gen_paper_plots.py`

The plots from the paper (and some additional plots) should now be placed in the `figures` folder.
 
> **Note:** The genetic algorithm is a stochastic optimization procedure. Hence, the results, especially for the concrete models shown for Wnt might not exactly coincide with those shown in the paper. However, they should follow a similar trend.

> **Note:** In step 2, there may be some errors from LSODA regarding a failiure of integration. This is expected for the csindy/Wnt experiments.
> If this stops the execution of further experiments, comment out the experiments already finished (you can easily check this with something like `ls experiment_results/*/*`) in `run_experiments.py` and run the script again.
> This will continue the experiments where they were left off.
> Generally, it is recommended to store the standard output in a file for later analysis of where a possible failure may have been using e.g. `(unbuffer python experiments/run_experiments.py 2>&1) | tee experiment_log.txt`.

## :file_folder: Overview of Contents

The following table provides an overview over the contents of the `src/` directory. At the top of each file, there is also a brief explanation regarding its purpose.

| Folder/file              | Content/Purpose                                                                     |
| ------:                  | :--------                                                                           |
| `reaction.py`            | Representation of a single biochemical reaction.                                    |
| `reaction_library.py`    | Representation of a library of biochemical reactions.                               |
| `reaction_enumerator.py` | Functions to list a complete library of reactions within some provided constraints. |
| `model_groundtruths.py`  | Benchmark model definitions.                                                        |
| `wnt.py`                 | Definition of the WNT model (in its own file for better overview).                  |
| `gen_data.py`            | Script to generate synthetic reference data from the ground truth models.           |
| `integrate.py`           | Some functions for simulating a reaction system and collecting trajectory data.     |
| `evolution/*`            | Python code implementing the genetic algorithm and associated helpers.              |
| `sindy/*`                | Python code implementing sparse regression with (coupled-)SINDy.                    |
| `experiments/*`          | Definitions of the experiments and plotting scripts run for the paper.              |
| `experiment_results/`    | Output directory for experimental data (empty by default).                          |
| `data/`                  | Output directory for synthetic reference data (empty by default).                   |
| `figures/`               | Output directory for figures.                                                       |

## :balance_scale: License

This project is licensed under the MIT License contained in `LICENSE`, unless indicated otherwise at the top of a file.

## :page_facing_up: Cite

```
Kreikemeyer, J.N., Burrage, K., Uhrmacher, A.M. (2024). Discovering Biochemical Reaction Models by Evolving Libraries.
In: Gori, R., Milazzo, P., Tribastone, M. (eds) Computational Methods in Systems Biology. CMSB 2024.
Lecture Notes in Computer Science, vol 14971. Springer, Cham.
https://doi.org/10.1007/978-3-031-71671-3_10
```
