This repository contains custom-written code for the analyses reported in the
article: Stephani, T., Waterstraat, G., Haufe, S., Curio, G., Villringer, A., Nikulin, V.V. (2020). Temporal Signatures
of Criticality in Cortical Excitability as probed by early somtosensory responses. JNeurosci.

The Folder "Scripts" contains all analysis scripts which are structured as
follows:

1) Preprocessing:
*_preprocessing.m

2) Single-trial extraction (CCA) and scaling analysis of the early SEP (DFA):
*_CCA_DFA.m

3) Source reconstruction:
(a) *_head_models_brainstorm.m
(b) *_source_reconstruction_eLoreta.m

4) Control measures:
(a) *_DFA_periphery.m
(b) *_DFA_thalamus.m
(c) *_SNR_estimation.m
(d) *_SNR_simulations.m
(e) *_filter_simulation_DFA.m

5) Relationship between alpha activity and early SEP dynamics:
(a) *_prestimulus_alpha_amplitude_relationship_and_DFA.m
(b) *_DFA_continuous_alpha.m

6) Statistics:
(a) *_cluster_statistics_DFA.m
(b) *_model_comparisons_DFA.m
(c) *_LME_models_R.R

7) Figures displayed in manuscript:
*_figures.m

For most of the .m scripts, EEGLAB is required, as well as Brainstorm for the
source reconstruction.

The folder "Functions" contains Matlab functions which
are called from the scripts.





