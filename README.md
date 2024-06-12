# Evolving Libraries (evolib)
 
Code artifacts for the paper "Discovering Biochemical Reaction Models by Evolving Libraries" to appear at the CMSB'24.

## :cd: Setup

### Required Software

- Linux; tested on Fedora 38 Workstation.
  - Should generally also work on Windows but some adaptations may be required
- Python; tested on version 3.11.9
- The dependencies listed in `requirements.txt`
  - can be installed automatically, see the next section

### Installation
(for Linux)

1. Clone the repository:
```
git clone <url>
```
2. Create a virtual environment and install the dependencies
```
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```
3. (locally) install the evolib python package in live mode
```
cd src
pip install -e .
```
(if this fails, you may need to start again from step 1, but this time update setuptools in step 2 with `pip install --upgrade setuptools`)

You should now have a virtual environment with all the dependencies and the module `evolib` (this repository) installed.

## :bar_chart: Reproduce Results

0. Navigate to the source directory: `cd src`
1. Generate the datasets: `python gen_data.py`.
2. Run the experiments: `python run_experiments.py`. This will run the experiments one after another and store the data in `src/experiment_results/<method>/<model>/<timestamp>/`.
3. Open the file `gen_paper_plots.py` and fill in the correct folder names (timestamps) at the top. You can use `ls experiment_results/*/*` in the `src` directory to get a good overview.
4. Plot the results: `python gen_paper_plots.py`
 
> **Note:** The genetic algorithm is a stochastic optimization procedure. Hence, the results, especially for the concrete models shown for Wnt might not exactly coincide with those shown in the paper. However, they should follow a similar trend.

> **Note:** In step 2, there may be some errors from LSODA regarding a failiure of integration. This is expected.
> If this stops the execution of further experiments, comment out the experiments already finished (you can easily check this with something like `ls experiment_results/*/*`) in `run_experiments.py` and run the script again.
> This will continue the experiments where they were left off.
> Generally, it is recommended to store the standard output in a file for later analysis of where a possible failure may have been using e.g. `(unbuffer python experiments/run_experiments.py 2>&1) | tee experiment_log.txt`.

## :file_folder: Overview of Contents

The following table provides an overview over the contents of the `src/` directory.

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

This project is liecensed under the MIT License contained in `LICENSE`, unless indicated otherwise at the top of a file.

## :page_facing_up: Cite

```
Justin N. Kreikemeyer, Kevin Burrage, Adelinde M. Uhrmacher. "Discovering 
Biochemical Reaction Models by Evolving Libraries", 2024. To appear in 
proceedings of the 22nd International Conference on Computational Methods
in Systems Biology, CMSB 2024.
```
