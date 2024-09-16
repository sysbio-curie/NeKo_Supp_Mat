# Use Case 1 - Medulloblastoma Sub-group Network

This folder contains all the scripts, notebooks, and figures used and generated for the original publication of NeKo.

The folder is organized as follows:

## Data
This folder contains the data used for the medulloblastoma analysis. The data is organized as follows:

- `41586_2016_BFnature16546_MOESM106_ESM.xlsx`: Supplementary material from a 2016 Nature publication related to the analysis.
- `antoine_discrete_data.csv`: Discrete data used in the analysis, processed by Antoine.
- `antoine_gene_type.csv`: Information about different gene types categorized by Antoine.
- `data_MB_genes_Northcott2017.xlsx`: Data on medulloblastoma genes from the Northcott 2017 study.
- `gene_per_subgroup_correlation.csv`: Correlation data for genes categorized per subgroup.
- `gene_per_subgroup_correlation_sorted.csv`: A sorted version of the gene per subgroup correlation data.
- `gene_per_subgroup_correlation_sorted_melt.csv`: A 'melted' version of the sorted gene per subgroup correlation data.
- `gene_per_subgroup_unique_correlation.csv`: Unique correlation data for genes per subgroup.
- `NIHMS781143-supplement-ST5.xlsx`: Supplemental material related to the NIHMS781143 record.
- `Northcott_Lin_genes.csv`: Gene expression data specific to the medulloblastoma subgroup as outlined in the Northcott and Lin studies.
- `script_nice_table.R`: An R script for creating well-formatted tables.
- `SIGNOR_complexes.csv`: Data on protein complexes from the SIGNOR database.
- `signor_db.tsv`: Whole SIGNOR database.

## Huri Analysis
This directory contains the following subfolders and their respective files:

### Data

1. `HI-union_translated.csv`: Translated version of the HI-union dataset.
2. `HI-union.tsv`: Tab-separated values file of the HI-union dataset.
3. `HURI_full_data_trasnlated.csv`: Translated version of the complete HURI dataset.
4. `HuRI_translated.csv`: Translated version of the HuRI dataset.
5. `HuRI.tsv`: Tab-separated values file of the HuRI dataset.
6. `legacy`: Contains older versions of the datasets.
7. `Lit-BM.tsv`: Tab-separated values file of the Lit-BM dataset.
8. `Test_space_screens-19_translated.csv`: Translated version of the Test space screens-19 dataset.
9. `Test_space_screens-19.tsv`: Tab-separated values file of the Test space screens-19 dataset.

### Figures

1. `group3_HURI.pdf`: PDF network figure for group 3 in the HURI analysis.
2. `group4_HURI.pdf`: PDF network figure for group 4 in the HURI analysis.
3. `SSH_HURI.pdf`: PDF network figure for SSH in the HURI analysis.
4. `WNT_HURI.pdf`: PDF network figure for WNT in the HURI analysis.

### Notebooks

1. `All_groups.ipynb`: Notebook summarizing the full group analysis in HURI.
2. `G3_HURI.ipynb`: Notebook for group 3 analysis in HURI.
3. `G4_HURI.ipynb`: Notebook for group 4 analysis in HURI.
4. `groups_comparison.ipynb`: Notebook comparing various groups within the HURI dataset.
5. `Huri_pre_processing.ipynb`: Notebook for pre-processing HURI data.
6. `SSH_HURI.ipynb`: Notebook for SSH analysis in HURI.
7. `WNT_HURI.ipynb`: Notebook for WNT analysis in HURI.

### SIF Files

1. `all_groups_HURI.sif`: SIF file for all groups in the HURI analysis.
2. `group3_HURI.sif`: SIF file for group 3 in the HURI analysis.
3. `group4_HURI.sif`: SIF file for group 4 in the HURI analysis.
4. `SHH_HURI.sif`: SIF file for SHH in the HURI analysis.
5. `WNT_HURI.sif`: SIF file for WNT in the HURI analysis.

## Northcott et al.
This directory contains the following subfolders and their respective files:

### Cytoscape Analysis

1. `Analysis_Lin_groups.cys`: Cytoscape session file focusing on Lin groups.
2. `Analysis_Northcott_groups.cys`: Cytoscape session file focusing on Northcott groups.

### Figures

1. **Group 3**: Contains figures related to Group 3.
   - `Group3_omnipath_Northcott.pdf`: PDF network figure for group 3 using Omnipath.
   - `Group3_signor_Northcott.pdf`: PDF network figure for group 4 using Signor.
2. **Group 4**: Contains figures related to Group 4.
   - `Group4_omnipath_Northcott.pdf`: PDF network figure for group 4 using Omnipath.
   - `Group4_signor_northcott.pdf`: PDF network figure for group 4 using Signor.
3. **SHH**: Contains figures related to SHH group.
   - `SHH_omnipath_Northcott.pdf`: PDF network figure for SHH group using Omnipath.
   - `SHH_signor_northcott.pdf`: PDF network figure for SHH group using Signor.
4. **Venn Diagrams**: Contains Venn diagrams.
   - **Group 3**: Venn diagrams related to Group 3.
      - `Venn_2_Group3_data_omnipath.png`: Venn diagram comparing Group 3 data from Omnipath.
      - `Venn_2_Group3_data_signor.png`: Venn diagram comparing Group 3 data from Signor.
      - `Venn_2_Group3_signor_omnipath.png`: Venn diagram comparing Group 3 from Signor and Omnipath.
   - **Group 4**: Venn diagrams related to Group 4.
      - `Venn_2_Group4_data_omnipath.png`: Venn diagram comparing Group 4 data from Omnipath.
      - `Venn_2_Group4_data_signor.png`: Venn diagram comparing Group 4 data from Signor.
      - `Venn_2_Group4_signor_omnipath.png`: Venn diagram comparing Group 4 from Signor and Omnipath.
   - **SHH**: Venn diagrams related to SHH.
      - `Venn_2_SHH_data_omnipath.png`: Venn diagram comparing SHH data from Omnipath.
      - `Venn_2_SHH_data_signor.png`: Venn diagram comparing SHH data from Signor.
      - `Venn_2_SHH_signor_omnipath.png`: Venn diagram comparing SHH from Signor and Omnipath.
   - **WNT**: Venn diagrams related to WNT.
      - `Venn_2_wnt_data_omnipath.png`: Venn diagram comparing WNT data from Omnipath.
      - `Venn_2_wnt_data_signor.png`: Venn diagram comparing WNT data from Signor.
      - `Venn_2_wnt_signor_omnipath.png`: Venn diagram comparing WNT from Signor and Omnipath.
5. **WNT**: Contains figures related to WNT.
   - `WNT_omnipath_Northcott_bfs.pdf`: PDF network figure for WNT group using Omnipath and the BFS algorithm.
   - `WNT_omnipath_northcott.pdf`: PDF network figure for WNT group using Omnipath.
   - `WNT_signor_northcott.pdf`: PDF network figure for WNT group using Signor.

### Northcott Analysis
- `Northcott_analysis.ipynb`: Jupyter Notebook for the Northcott analysis.

### Notebooks

1. `comparison_group3_group4.ipynb`: Notebook comparing Group 3 and Group 4.
2. `comparison_group3_signor_omnipath.ipynb`: Notebook comparing Group 3 from Signor and Omnipath.
3. `comparison_group4_signor_omnipath.ipynb`: Notebook comparing Group 4 from Signor and Omnipath.
4. `comparison_shh_signor_omnipath.ipynb`: Notebook comparing SHH from Signor and Omnipath.
5. `comparison_wnt_signor_omnipath.ipynb`: Notebook comparing WNT from Signor and Omnipath.
6. **Group 3**: Contains notebooks and data specific to Group 3.
   - `fill_the_net.ipynb`: Notebook for filling the network data.
   - `Group3_generic.ipynb`: Generic notebook for Group 3 analysis.
   - `Group3_signor.ipynb`: Notebook for Group 3 Signor analysis.
   - `net_signor_g3.pdf`: PDF network figure for group 3.
7. **Group 4**: Contains notebooks and data specific to Group 4.
   - `fill_the_net.ipynb`: Notebook for filling the network data.
   - `Group4_generic.ipynb`: Generic notebook for Group 4 analysis.
   - `Group4_signor.ipynb`: Notebook for Group 4 Signor analysis.
   - `net_signor.pdf`: PDF network figure for group 4.
8. `groups_comparison.ipynb`: Notebook comparing various groups.
9. `Northcott_Lin_genes`: Contains data related to Northcott Lin genes.
10. **SHH Group**: Contains notebooks and data specific to the SHH group.
   - `fill_the_net.ipynb`: Notebook for filling the network data.
   - `net_signor.pdf`: PDF network figure for SHH group.
   - `SHH_generic.ipynb`: Generic notebook for SHH analysis.
   - `SHH_signor.ipynb`: Notebook for SHH Signor analysis.
11. **WNT Group**: Contains notebooks and data specific to the WNT group.
   - `fill_the_net.ipynb`: Notebook for filling the network data.
   - `net_signor.pdf`: PDF network figure for WNT group.
   - `pypath_log`: Log file related to pypath operations.
   - `WNT_generic.ipynb`: Generic notebook for WNT analysis.

### SIF Files

1. **Group 3**: Contains SIF files related to Group 3.
   - `Group3_omnipath_Northcott.sif`: SIF file for Group 3 Omnipath Northcott analysis.
   - `group3_signor_Northcott.sif`: SIF file for Group 3 Signor Northcott analysis.
   - `Signor_with_complexes_Northcott.sif`: SIF file including complexes for Group 3 Signor Northcott analysis.
2. **Group 4**: Contains SIF files related to Group 4.
   - `Group4_omnipath_Northcott.sif`: SIF file for Group 4 Omnipath Northcott analysis.
   - `Group4_signor_northcott_radial.sif`: Radial SIF file for Group 4 Signor Northcott analysis.
   - `Group4_signor_northcott.sif`: SIF file for Group 4 Signor Northcott analysis.
   - `Signor_with_complexes_Northcott.sif`: SIF file including complexes for Group 4 Signor Northcott analysis.
3. **SHH**: Contains SIF files related to SHH.
   - `SHH_omnipath_Northcott.sif`: SIF file for SHH Omnipath Northcott analysis.
   - `SHH_signor_northcott.sif`: SIF file for SHH Signor Northcott analysis.
   - `Signor_with_complexes_Northcott.sif`: SIF file including complexes for SHH Signor Northcott analysis.
4. **WNT**: Contains SIF files related to WNT.
   - `Signor_with_complexes_Northcott.sif`: SIF file including complexes for WNT Signor Northcott analysis.
   - `WNT_omnipath_northcott.sif`: SIF file for WNT Omnipath Northcott analysis.
   - `WNT_signor_Northcott.sif`: SIF file for WNT Signor Northcott analysis.

## How to Navigate

- **Cytoscape Analysis**: This folder contains Cytoscape session files for different groups and analyses.
- **Figures**: This folder contains subfolders for different groups (Group 3, Group 4, SHH, WNT) and figures related to the analysis.
- **Notebooks**: Contains various Jupyter Notebooks used for performing the analysis. It includes comparative analysis among different groups, detailed analysis within groups, and logs related to the process.
- **SIF Files**: Contains SIF (Simple Interaction Format) files which are used to represent the interactions within each group for various analyses.