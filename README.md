# Synaptic-Proteome-Diversity
This repository includes the code develop for the scientifc article entiteled 'Synaptic proteome diversity is primarily driven by gene regulation of glutamate receptors and their regulatory proteins'. This work is freely available at the preprint repository bioRxiv, which can be downloaded from the following link: https://www.biorxiv.org/content/10.1101/2024.04.04.588090v1

Article Abstract,

Electrophysiological features of excitatory synapses vary widely throughout the brain, granting neuronal circuits the ability to decode and store diverse patterns of information. Synapses formed by the same neurons have similar electrophysiological characteristics, belonging to the same type. However, these are generally confined to microscopic brain regions, precluding their proteomic analysis. This has greatly limited our ability to investigate the molecular basis of synaptic physiology. Here we introduce a procedure to characterise the proteome of individual synaptic types. We reveal a remarkable proteomic diversity among the synaptic types of the trisynaptic circuit. Differentially expressed proteins participate in well-known synaptic processes, controlling the signalling pathways preferentially used among diverse synapses. Noteworthy, all synaptic types differentially express proteins directly involved in the function of glutamate receptors. Moreover, neuron-specific gene expression programs would participate in their regulation. Indeed, genes coding for these proteins exhibit such distinct expression profiles between neuronal types that they greatly contribute to their classification. Our data is an important resource for exploring the molecular mechanisms behind electrophysiological properties of different hippocampal synaptic types. Our combined analysis of proteomics and transcriptomics data uncovers a previously unrecognised neuron-specific transcriptomic control of synaptic proteome diversity, directed towards the regulation of glutamate receptors and their regulatory proteins.

Authors and Affiliations,

Rita Reig-Viader1,2*, Diego del Castillo-Berges1,2*, Albert Burgas-Pau1,2,3,4*, Daniel Arco-Alonso1*, David Ramos-Vicente1,5, Carlos Sindreu6, Àlex Bayés1,2

1.	Molecular Physiology of the Synapse Laboratory, Institut de Recerca Sant Pau (IR Sant Pau), Barcelona, Spain.
2.	Universitat Autònoma de Barcelona (UAB), Bellaterra (Cerdanyola del Vallès), Spain.
3.	Unitat mixta d’Investigació IRTA-UAB en Sanitat Animal. Centre de Recerca en Sanitat Animal (CReSA), Bellaterra (Cerdanyola del Vallès), Spain.  
4.	IRTA. Programa de Sanitat Animal. Centre de Recerca en Sanitat Animal (CreSA), Bellaterra, 08193, Catalonia. Spain.
5.	Neurodegenerative Diseases Research Group, Vall d’Hebron Research Institute, Centre for Networked Biomedical Research on Neurodegenerative Diseases (CIBERNED), Barcelona, Catalonia, Spain.
6.	Institut d'Investigacions Biomèdiques August Pi i Sunyer (IDIBAPS), Barcelona, Spain.

* Equal correspondance
   
Corresponding author: Àlex Bayés, Molecular Physiology of the Synapse Laboratory, Institut de Recerca Sant Pau (IR Sant Pau), Barcelona, Spain., C/Sant Quintí, 77-79, 08041 Barcelona, Spain.
Email: abayesp@santpau.cat 

# Information on files in this repository:

In this repository we have deposited files relevant to the above-mentioned manuscript. This code allows to perform all bioinformatic analysis performed in the manuscript. 

Short description on the porpouse of each file:

Source_Data_1_Iteration_Classes.R: R script to iterate the statistical analysis performed with Seurat to identify genes differentially expressed between neuronal classes.

Source_Data_2_Iteration_Types.R: R script to iterate the statistical analysis performed with Seurat to identify genes differentially expressed between neuronal types.

Source_Data_3_Analysis_Classes.R: R script to generate data tables and graphs for genes differentially expressed between neuronal classes.  

Source_Data_4_Analysis_Types_CA3.R: R script to generate data tables and graphs for genes differentially expressed between neuronal Types in the CA3 Class. In the test_set folder you can find the dataset of DE genes between CA3 neurons that will allow to test this code. 

Source_Data_5_Split_Types.R: R script to obtained data from a subset of neuronal types from the entire transcriptomic database provided by the ABCA.

Source_Data_6_ Random_Forest_Umaps.ipynb: Python code to perform the Random Forest analysis on transcriptomic data from the ABCA.

Source_data_7_pathfinR_analysis_DE_Genes_Neuronal_Classes.r: R script to perform the pathfinder analysis and to generate the heatmaps from transcriptomics data of neuronal classes (ABCA).

Source_data_8_pathfinR_analysis_DE_genes_Neuronal_types.r: R script to perform the pathfinder analysis and to generate the heatmaps from transcriptomics data of neuronal types (ABCA).

We have also included a test set (see test_set folder) to explore how we call and identify differentially expressed synaptic genes. This test set include expression data among neuronal types from the CA3 class.

If you require more information on how to operate these do not hesitate contacting us: Àlex Bayés: abayesp@santpau.cat
