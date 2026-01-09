# Chemputer and Chemputation -- A Universal Chemical Compound Synthesis Machine

**Authors:** Leroy Cronin, Sebastian Pagel, Abhishek Sharma  
**Affiliation:** [University of Glasgow](https://www.chem.gla.ac.uk/cronin/)

## Abstract

Chemputation reframes synthesis as the programmable execution of reaction code on a universally re-configurable hardware graph. Here we prove that a chemputer equipped with a finite, but extensible, set of reagents, catalysts and process conditions, together with a chempiler that maps reaction graphs onto hardware, is universal: it can generate any stable, isolable molecule in finite time and in analytically detectable quantity, provided real-time error correction keeps the per-step fidelity above the threshold set by the molecule's assembly index. The proof is constructed by casting the platform as a Chemical Synthesis Turing Machine (CSTM). The CSTM formalism supplies (i) an eight-tuple state definition that unifies reagents, process variables (including catalysts) and tape operations; (ii) the Universal Chemputation Principle; and (iii) a dynamic-error-correction routine ensuring fault tolerant execution. Linking this framework to assembly theory strengthens the definition of a molecule by demanding practical synthesizability and error correction becomes a prerequisite for universality. We validate the abstraction against >100 XDL programs executed on a modular chemputer rigs spanning single step to multi-step routes. Mapping each procedure onto CSTM shows that the cumulative number of unit operations grows linearly with synthetic depth. Together, these results elevate chemical synthesis to the status of a general computation: algorithms written in XDL are compiled to hardware, executed with closed-loop correction, and produce verifiable molecular outputs. By formalising chemistry in this way, the chemputer offers a path to shareable, executable chemical code, interoperable hardware ecosystems, and ultimately a searchable, provable atlas of chemical space. 

<p align="center">
  <img src="https://github.com/user-attachments/assets/56f085b8-2dac-405b-90fa-28556664d904" />
</p>

## Data Availability

The dataset provided in the `data/` directory has been adapted from the following publication:

> *Science* **2022**, *377*. DOI: [10.1126/science.abo0058](https://www.science.org/doi/10.1126/science.abo0058)


## Notebooks

The `notebooks/` directory contains Jupyter notebooks for demonstrating the usage of the codebase and reproducing figures.

- `Figure6.py`: Python script to generate Figure 6.
- `ChemicalTuringMachine.nb`: Mathematica notebook for the simulation of a chemical turing machine
- `reaction_display.ipynb`: Examples of visualizing chemical reactions and molecules using `rxnframe`.
- `Figure7.ipynb`: Code to generate subfigures of Figure 7.

## Installation

This project uses [uv](https://docs.astral.sh/uv/getting-started/installation/) for dependency management. To install the dependencies, run:

```bash
uv sync
```

**Note:** To run the Mathematica notebooks (`.nb` files), you will need a valid installation of Wolfram Mathematica.

## Citation

If you use this code or data, please cite the following:

> [Chemputer and Chemputation -- A Universal Chemical Compound Synthesis Machine](https://arxiv.org/abs/2408.09171)

