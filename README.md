# Mindfulness-based therapy improves brain functional network reconfiguration efficiency

## Reference
>Yue, W.L., Ng, K.K., Koh, A.J.K., Perini, F., Doshi, K., Zhou, J.H., & Lim, J. (2023). Mindfulness-based therapy improves brain functional network reconfiguration efficiency. Translational Psychiatry.

## Usage
`FC_simlarity.m` is used to derive FC similarity between rest and task.

`lme_analyses.m` is used to obtain time effects and group x time interactions (from linear mixed models) for FC similarity and behavioral measures.

`lme_permutations.m` is used to create null distributions using permutations for time effects and group x time interactions.

`lm_analyses.m` is used to obtain associations (from linear models) between changes in FC similarity and behavioral measures after intervention.

`lm_permutations.m` is used to create null distributions using permutations for associations.


Outputs from earlier scripts are used as inputs in later scripts in this analysis pipeline.


Derivation of FC similarity measures is based on: 
>Schultz, D. H., & Cole, M. W. (2016). Higher intelligence is associated with less task-related brain network reconfiguration. Journal of neuroscience, 36(33), 8551-8561. https://doi.org/10.1523/JNEUROSCI.0358-16.2016
