# Use Case 2 - DrugLogics Pipeline

This folder contains all the scripts, notebooks, and figures that were used and generated for the original publication of NeKo.

To utilize this folder effectively, you will need both the `druglogics-synergy` and `synergy-tutorial` repositories, available at the following links:

- [druglogics-synergy](https://github.com/druglogics/druglogics-synergy)
- [synergy-tutorial](https://github.com/druglogics/synergy-tutorial)

The folder is organized as follows:

## usecase_neko_ags

- `20170321_HGNC-nodelist.xlsx`: An Excel file containing a list of nodes (genes/proteins) with their corresponding HGNC (HUGO Gene Nomenclature Committee) identifiers.
- `cascade_net_refined.sif`: A refined SIF file representing the cascade network used in the analysis.
- `hgnc_conversion_key_for_comparison_cascade_1_vs_neko_topology.xlsx`: An Excel file providing a conversion key to compare the topology between Cascade1.0 and NeKo networks, using HGNC identifiers.
- `logic_model_network_cascade.sif`: A SIF file representing the logic model network specifically for the cascade analysis.
- `logic_model_network.sif`: A SIF file representing the logic model network used built with NeKo.
- `neko_config.tab`: A tab-separated configuration file for the NeKo analysis, specifying parameters and settings.
- `neko_drugpanel.tab`: A tab-separated file containing information about the drug panel used in the NeKo analysis.
- `neko_modeloutputs.tab`: A tab-separated file listing the output results from the NeKo model simulations.
- `neko_perturbations.tab`: A tab-separated file detailing the perturbations applied during the NeKo analysis.
- `neko_training_data.tab`: A tab-separated file containing the training data used for the NeKo model.
- `signor_db.tsv`: A tab-separated values file containing data from the Signor database.
- `Signor_filtered.tsv`: A filtered tab-separated values file from the Signor database.

## usecase_neko_ags
- Folder containing all the generated figures from the analysis and cytoscape

## Comparison_cascade_neko.ipynb

- Notebook containing the analysis of the network built with NeKo and its comparison with the one from Cascade 1.0

## use_case_ags_synergy_network_building.ipynb

- Notebook containing the code to build the network using NeKo

## neko_data.tab

- A tab-separated file with the values from the druglogics pipeline to predict the effect of drug synergies 

## cytoscape_session_use_case_2.cys

- Cytoscape session containing the topological analysis of both the Cascade network and the NeKo one.

## use_case_2_neko.R

- R script for the analysis of the NeKo network.