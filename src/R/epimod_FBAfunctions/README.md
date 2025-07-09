
# Using _UnifiedGreatMod_ to integrate constrained-based models into an unique framework

To integrate a constrained-based models into _GreatMod_ framework, we developed the "epimod_FBAfunctions" R package. It offers essential functions to read and modify COBRA MAT files, translating metabolic network representations into a format suitable for the epimod functions based on the GLPK linear programming solver other integrate Sobol's variance-based Sensitivity Analysis (SA) to classify reactions in a FBA-model based on their impact on model outcomes.

It is an R package for supporting [GreatMod](https://qbioturin.github.io/epimod/) for the integration of metabolic models, represented as a Flux Balance Analysis (FBA) problem, into a dynamic model, represented as a set of ordinary differential equations (ODEs). Additionally, It can be used for exploring the sensitivity of input parameters in genome-wide models represented as a Flux Balance Analysis (FBA) problem, using Sobol variance-based Sensitivity Analysis as global strategy.

# GreatMod, FBA and Sobol variance-based Sensitivity Analysis

## FBA General Functions

_UnifiedGreatMod_ (https://qbioturin.github.io/epimod/) is a general modelling framework to simulate biological complex systems using an
intuitive graphical interface for model construction and automatic derivation of the low-level mathematical processes characterising
the system dynamics. Using **epimod_FBAfunctions**, gene expression data and dietary inputs can be integrated to set flux constraints. Our computational pipeline allows for (i) computing differential reaction activities from transcriptomics to predict how changes in metabolic enzyme expression affect metabolic fluxes and (ii) constraining boundary reactions based on metabolite concentrations from dietary inputs. The critical aspect is selecting the integration method, implemented as described by Kashaf et al. (2017), which involves evaluating Gene-Protein-Reaction rules to derive reaction expression values from gene expression data.

1. **class_generation** Flux Balance Analysis model generation for GreatMod
2. **FBAgreatmodeClass** generates the FBA_greatmod class in which all the elements needed for defining the FBA model to pass to GreatMod are saved. It takes in input a Matlab file
3. **readMat** Reads matlab model and imports matlab cobra (ONLY) models into the R environment
4. **FBA4Greatmod.generation** generates the *FBA_greatmod* class in which all the elements needed for defining the FBA model to pass to GreatMod are saved. It takes in input a Matlab file.

## FBA_greatmod class methods

- **initialize** initializes an FBA_greatmod object with provided stoichiometric matrix, bounds, objective function, and optionally reaction and metabolite names.
- **setObjFun** sets the objective function in an FBA model to maximize the flux of a specified reaction.
- **getExchangesR** identifies and retrieves the exchange reactions in an FBA model.
- **getConstraints** retrieves the flux constraints for a specific reaction in an FBA model.
- **writeFBAfile** writes the FBA model's details to a file suitable for GLPK solver input.
- **setDiet** updates the dietary constraints in an FBA model based on input from a file or data frame.
- **setConstraints**  updates the flux constraints for a specific reaction in an FBA model.

## Global Sensitivity Analysis Functions

Also, **epimod_FBAfunctions** allows for Global SA. It is an analysis workflow ideal for any genome-wide constraint-based model combined with the random sampling scheme of parameters. The sampling method is a variance-based global SA that exploits Monte Carlo simulations to estimate the sensitivity indices of input parameters. The sampling involves generating a set of samples using Sobol's Quasi Random Numbers sequences using a base of two to form successively finer uniform partitions of the unit interval, reordering the coordinates in each dimension, and then computing the model's output for each sample. Global SA is used to stratify FBA model reactions according to their contribution to model outcomes. 

We applied Sobol's variance-based SA to calculate how much of the variance of output variables is explained by a given input and the interaction repercussion among inputs.

4. **ParallelFBA_sens** solves parallel Flux Balance Analysis optimizations given the parameter set generated with saltelli sampling.
5. **SA_FBA_R** implememts Sobol's variance-based SA Global Sensitivity Analysis using R.

## Installing epimod_FBAfunctions

To install **epimod_FBAfunctions** you can use use **devtools**:

```
install.packages("devtools")
library(devtools)
install_github("https://github.com/qBioTurin/epimod_FBAfunctions", ref="master")
```

## Repositories

The following list is a selection of project developed with **epimod** and **epimod_FBAfunctions**, providing both the necessary files and explanations to perform the analysis. 

### Transcriptional data onto Escherichia coli genome-scale metabolic model (iML1515) growing on different regime of carbon feeding

* [E.coli carbon regimes](https://github.com/qBioTurin/Ec_coli_modelling): _epimod_ integrates Flux Balance Analysis (FBA) and Petri Net dynamical models based on Ordinary Differential Equations (ODEs) into a unified framework. This repository includes an example demonstrating its capability to simulate multi-state systems by solving the FBA optimization problem with constraints on metabolite rate changes at specific times. _epimod_ can simulate _E. coli_'s metabolic output under various transcriptional programs and nutrient conditions, providing insights into metabolic reprogramming and metabolic network design.

### Modeling host and drug responses to Clostridium difficile infection
* [CDI](https://github.com/qBioTurin/EpiCell_CDifficile): We combined the ODE-based dynamical model with the Genome-Scale Metabolic Model (GSMM) into a unified multiscale hybrid model. _C. difficile_ strain CD196 GSMM was originally coded in the AGORA/MATLAB (1.03 version) format (Magnusdottir et al.,2016) and then elaborated using the R functions collected in [epimod_FBAfunctions](https://github.com/qBioTurin/epimod_FBAfunctions).
