
# UnifiedGreatMod: Hybrid Modeling of Microbial and Multicellular Systems

## Overview

This repository houses the **UnifiedGreatMod** framework, a novel, modular, and extensible computational platform for designing and simulating **unified hybrid models** of complex biological systems.Building on Extended Stochastic Petri Net formalisms, 'UnifiedGreatMod' integrates fine-grained **(Ordinary Differential Equation - ODE-based)** dynamics with coarse-grained **constraint-based (Flux Balance Analysis - FBA)** mechanistic modelling approaches.


'UnifiedGreatMod' addresses the inherent multi-scale complexity of biological systems by providing a mathematically rigorous and computationally efficient platform capable of handling interactions across different levels of biological organization. The framework supports dynamic modeling and analysis of:

* **Microbial Ecosystems**: From simplified communities to complex microbiotas. [cite: 9]
* **Host-Microbiota Interactions**: Deciphering the dynamic interplay between host and microbial components. [cite: 9]
* **Multicellular Tissue Environments**: Capturing emergent behaviors in complex cellular landscapes.
* **Synthetic and Hybrid Biological Constructs**: For bioengineering and system design applications.

## Core Concepts

### Hypernode Assembly

A **Hypernode** is a high-level computational construct within `UnifiedGreatMod` that encapsulates the interactions and internal dynamics of multiple biological units (referred to as **biounits**) within a unified, *hybrid, and harmonized multi-low-level formalism* based on **Extended Stochastic Petri Nets (ESPNs)**. Each hypernode functions as a composite ESPN graph where:

* **Nodes** (places $P$) represent biological entities such as molecular species (metabolites, signaling intermediates) or macroscopic quantities (cell populations, biomass concentrations). These can hold `tokens` representing the modeled entities of the system, defining its `marking`.
* **Transitions** ($T$) model both continuous (Ordinary Differential Equation - ODE-governed) and discrete (event-driven or Linear Programming - LP-constrained) processes or events, such as biochemical reactions, population dynamics, or metabolite exchange, including automated Flux Balance Analysis (FBA) resolution. These transitions are categorized into Mass Action (`Tma`) or general (`Tg`), with `Tg` further subdivided to include FBA-driven rates (`TgFBA`).
* **Edges** (input $I$ and output $O$ functions) encode causal dependencies, stoichiometric relationships, and flux multiplicities, defining how transitions alter the marking of places within the composite graph.

This abstraction enables:

* The simultaneous simulation of **constraint-based models**, and **fine-grained mechanistic models**, capturing time-dependent behaviors.
* The `harmonization` of heterogeneous biological components, such as distinct microbial species (each potentially represented by its Genome-scale Metabolic Model - GEM) and host cell types, within a consistent mathematical and computational structure. This is facilitated by dynamically coupled `FBA[]` commands within ESPN transitions, linking FBA flux solutions to ODE rates, and vice-versa.
* Encapsulation of intra- and inter-entities dynamics, allowing scalable and interpretable simulation of complex biological `ecosystems`, `tissues`, or `synthetic constructs` that exhibit emergent properties from multi-level interactions.

### Unified Biounit Handling

The **`biounits/`** directory serves as a general-purpose container for any biological unit involved in the simulation—regardless of whether it is a microbe, a cell type, a synthetic agent, or any composite entity.

* **Generalization**: Treatment unification of diverse biological components, making the framework applicable to both unicellular and multicellular systems.
* **Scalability**: Each biounit is self-contained, enabling modular expansion, reuse, or replacement.
* **Integration-Ready**: All biounits are compatible with the unified ESPN formalism, and their behavior is defined through standard formats (metabolic models, metadata and unified model configuration).

#### Examples of Unified Models (Hypernodes)

| Type | Example | Description | Biounit(s) |
| :---------------------- | :------------------------------- | :---------------------------------------------------------- |:---------------------------------------------------------- |
| Microbial cells | *E. coli* Growth | ODE-GEMs based model of a microorganism capturing dynamic metabolic responses to environmental changes. | *Escherichia coli* Str K-12 MG1655 |
| Microbial community | SIHUMIx (Simplified Human Intestinal Microbiota) [cite: 6357, 752] | ODE-GEMs model heterogeneous for a population microbial assemblies. | *Anaerostipes caccae* DSM 14662; *Bacteroides thetaiotaomicron* VPI 5482; *Bifidobacterium longum* NCC2705; *Blautia producta* DSM 2950; *Clostridium butyricum* DSM 10702; *Clostridium ramosum* VPI 0427 DSM 1402; *Escherichia coli* Str K-12 MG1655; *Lactobacillus plantarum* subsp. *plantarum* ATCC 14917 |
| Host cell types | Intestinal epithelial cell | ODE-GEMs based model of a human cell, particularly for studying host-pathogen interactions and nutrient absorption. | Human Intestinal Epithelial Cells (IECs) |
| Tissue | Health/Cancer multicell-type tissue | ODE-GEMs model for a tissue region or heterogeneous cellular states, often applied to study altered metabolic pathways in disease contexts. | Multiple human cell types within a tissue (e.g., cancer cells, stromal cells, immune cells) |

## Key Features & Methodological Advancements

'UnifiedGreatMod' integrates several advanced computational techniques to achieve its holistic modeling capabilities:

* **Automated FBA Workflow**: Automates the entire FBA process, from reading `.mat` files and converting them to optimized sparse `.txt` formats, to dynamically updating flux bounds during simulations based on time-varying conditions.
* **Sparse Matrix Optimization**: Leverages sparse matrix representation for metabolic networks, drastically reducing memory consumption and improving computational efficiency for large-scale GEMs.
* **Biomass Parameterization**: Allows for explicit parameterization of bacterial biomass (min, mean, max) directly within the FBA setup, enabling more realistic growth dynamics. [cite: 6]
* **Extended Stochastic Petri Nets (ESPNs)**: Utilizes ESPN as the core meta-formalism to represent complex interactions, allowing for custom rate functions and flexible integration of heterogeneous modeling paradigms.
* **Coloured Petri Nets (CPNs) for Community Modeling**: Introduces "colors" (e.g., GEMs ID, trophic layers, spatial coordinates) within the ESPN formalism to manage multiple interacting biounits in a single model. This enables the simulation of complex community structures and inter-species metabolic exchanges.
* **Modular R Interface (`epimod`, `epimod_FBAfunctions_enhanced`)**: Offers a user-friendly R interface to orchestrate model generation, parameterization, simulation, and analysis, promoting reproducibility and extensibility.

## Simulation Workflow

The `UnifiedGreatMod` framework orchestrates a multi-step simulation workflow, from initial data preparation to advanced model analysis.

1.  **Biounit Metadata Creation**:
    * For each unique biological unit (biounit) or species:
        * Prepare or acquire its Genome-scale Metabolic Model (GEM) in `.mat` format, placed in `models/<model_label>/data/raw_mats/`.
        * Create `metabolites_metadata.csv` and `reactions_metadata.csv` files, located in `models/<model_label>/data/biounits/<biounit_name>/`. These files detail the properties and annotations for each metabolite and reaction within the biounit's GEM.

2.  **Biounit Configuration**:
    * Define species-specific parameters (e.g., biomass `max`, `mean`, `min`, population `starv`, `dup`, `death` rates) and overall model configuration in `models/<model_label>/config/config_<model_label>.yaml`.
    * Specify global boundary conditions for the integrated metabolic network in `models/<model_label>/config/boundary_conditions.json`.

3.  **Model Generation & Preprocessing**:
    * Execute R scripts (e.g., `src/R/utils/run_model_test.R` or custom `setup_models.R`) to automate preprocessing steps:
        * Convert GEMs from `.mat` to an optimized sparse `.txt` format for each biounit, leveraging `epimod_FBAfunctions_enhanced::FBA4Greatmod.generation()` and `writeFBAfile()` for efficiency and reduced memory footprint.
        * Generate or load the primary Extended Stochastic Petri Net (ESPN) or Coloured Petri Net (CPN) `.PNPRO` file for the community model (e.g., `models/<model_label>/pnpro/<model_label>.PNPRO`). This PNPRO model utilizes `FBA[]` commands within transitions to link dynamically to the `.txt` metabolic models.
        * Prepare model-specific C++ functions by programmatically injecting parameters into generic templates from `src/cpp/templates/general_functions_template.cpp`.

4.  **Solver Compilation**:
    * Use `epimod::model_generation()` to compile the `.PNPRO` file and the associated C++ code into a high-performance executable solver (`.solver` file). This step transforms the high-level PN formalism into low-level mathematical processes (e.g., ODEs).

5.  **Simulation Execution**:
    * Run hybrid dynamic simulations using `epimod::model_analysis()`. This function takes the compiled `.solver` file, initial conditions, and simulation parameters (e.g., final time `f_time`, sampling interval `s_time`, ODE solver `solver_type`) to generate time-course trajectories.The simulation dynamically resolves coupled ODEs and FBA problems based on defined synchronization criteria.

6.  **Sensitivity Analysis**:
    * Prior to full-scale simulations, perform sensitivity analysis using `epimod::sensitivity_analysis()`. This identifies key model parameters that most significantly influence model outcomes, allowing for focused parameter estimation and hypothesis testing.

## Applications & Case Studies

This framework is designed for broad applicability and has been applied to a range of complex biological systems:

* **Minimal Microbiota (EcCb)**: Functional modeling of simplified microbial communities (e.g., *Escherichia coli* and *Clostridium butyricum* doublet), investigating cross-feeding and competition.
* **Simplified Human Intestinal Microbiota (SIHUMIx)**: Community-scale modeling of more complex gut microbial consortia, focusing on metabolic outputs like short-chain fatty acid (SCFA) production and inter-species interactions.

## Project Philosophy

* **Unified Formalism**: All biological processes—discrete or continuous, metabolic or cellular—are encoded as hybrid Petri nets.
* **Modular Construction**: Each `biounit` is a plug-and-play module.
* **Reproducibility**: Every simulation component is version-controlled and traceable, adhering to best practices for computational biology.

## Getting Started

To get started with UnifiedGreatMod:

1.  **Clone this repository:**
    ```bash
    git clone [https://github.com/your_username/UGM_Microbiota.git](https://github.com/<your_username>/UGM_Microbiota.git)
    cd UGM_Microbiota
    ```
2.  **Install Dependencies:**
    * **R Packages**: Ensure `epimod` and `yaml`, `jsonlite` are installed.
        ```R
        install.packages(c("epimod", "yaml", "jsonlite"))
        # Install your local enhanced epimod_FBAfunctions_enhanced package
        install.packages("src/R/epimod_FBAfunctions_enhanced", repos = NULL, type = "source")
        ```
    * **External Tools**: Make sure `GreatSPN` (for PNPRO creation/editing) and `GLPK` (GNU Linear Programming Kit) are installed and accessible in your system's PATH. Refer to the `GreatSPN` documentation for specific installation instructions[cite: 3].
3.  **Configure a Model:**
    * Navigate to `models/EcCb/`.
    * Review `config/config_minimal_doublet.yaml` and `config/boundary_conditions.json`.
4.  **Run a Test Simulation:**
    * Open R, set your working directory to the repository root.
    * Execute the test script for the EcCb model:
        ```R
        source("src/R/utils/run_model_test.R")
        run_model_test(model_label = "EcCb") # Or "minimal_doublet" depending on how you structure subfolders
        ```