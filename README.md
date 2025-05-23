# Unified Hybrid Modeling Appraoches of Microbial and Multicellular Systems

## Overview

This repository provides a modular and extensible framework for designing and simulating **unified hybrid models** of complex biological systems. These models are built as **Petri net-based hypernodes**, which integrate both mechanistic (ODE-based) and constraint-based (FBA) dynamics.

The system supports dynamic modeling of:

* Microbial ecosystems
* Host-microbiota interactions
* Multicellular tissue environments
* Synthetic and hybrid biological constructs

## Hypernode Assembly

Hypernode refers to a high-level computational construct that encapsulates the interactions and internal dynamics of multiple biological units—referred to as biounits—within a unified hybrid Petri net formalism. Each hypernode functions as a composite Petri net graph where:

Nodes represent biological entities such as metabolites, cell populations, or signaling intermediates.

Transitions model both continuous (ODE-governed) and discrete (event-driven or LP-constrained) processes.

Edges encode causal dependencies or stoichiometric relationships.

Hypernodes enable:

The simultaneous simulation of constraint-based models (e.g. flux balance analysis) and fine-grained mechanistic models (e.g. kinetic pathways).

Integration of heterogeneous biological components—such as microbes and host cell types—within a consistent mathematical structure.

Encapsulation of intra- and inter-agent dynamics, allowing scalable and interpretable simulation of ecosystems, tissues, or synthetic constructs.

This abstraction supports multiscale systems biology by bridging metabolic, cellular, and population-level phenomena.

## Key Feature: Unified Agent-Based Modeling

A central innovation of this framework is the introduction of a **biounits/** directory. This directory serves as a general-purpose container for any biological unit involved in the simulation—regardless of whether it is a microbe, a cell type, a synthetic agent, or any composite entity.

### Why `biounits/`?

* **Generalization**: It unifies the treatment of cell types and organisms, making the framework applicable to both unicellular and multicellular systems.
* **Scalability**: Each biounit is self-contained, enabling modular expansion, reuse, or replacement.
* **Integration-ready**: All biounits are compatible with the unified Petri net formalism, and their behavior is defined through standard formats (.mat, .csv).

### Examples of Biounits

| Type                  | Example                       |
| --------------------- | ----------------------------- |
| Microbial species     | `Escherichia_coli_SE11`       |
| Host cell types       | `Intestinal_epithelial_cell`  |
| Engineered constructs | `Synthetic_butyrate_producer` |
| Tissue compartments   | `Colon_crypt_base_cell`       |

## Simulation Workflow

1. **Organism/Cell Configuration**:

   * Define metabolic models (.mat files).
   * Provide `metabolites_metadata.csv` and `reactions_metadata.csv`.

2. **Hypernode Assembly**:

   * Connect multiple biounits into a coherent Petri net.
   * Specify population dynamics (duplication, death, starvation).

3. **Validation & Rendering**:

   * Validate arcs using `validate_pnpro()`.
   * Optionally visualize or refine layout.

4. **Simulation Execution**:

   * Run hybrid dynamic simulations with both ODEs and FBA.

## Applications

This framework has been designed with broad applicability in mind, including but not limited to:

* Functional modeling of **gut microbiota**.
* Simulating **host-pathogen interactions**.
* Investigating **immune cell dynamics** within a tissue niche.
* Testing **synthetic biology circuits** in microbial consortia.

## Project Philosophy

* **Unified Formalism**: All biological processes—discrete or continuous, metabolic or cellular—are encoded as hybrid Petri nets.
* **Modular Construction**: Each `biounit` is a plug-and-play module.
* **Reproducibility**: Every simulation component is version-controlled and traceable.

## Getting Started

* Configure your `biounits/` directory with appropriate subdirectories.
* Modify `config_minimal_doublet.yaml` for species setup.
* Use `main_doublet.R` to generate, validate, and simulate your system.

---

For more advanced use cases, including multi-organ or spatiotemporal modeling, please refer to the documentation and case studies in the `/demo` folder.
