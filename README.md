# Sex differences in autism

This repository is a code appendix for our in submission paper on sex differences in autism. The code is intended for examination only, if you would like to
reproduce the analysis please contact me (@cfhammill) and we can discuss how you can adapt this code to run on another computing system.

- ./setup contains code to organize CIVET results and create the imputed IQ data sets
- ./cortical contains the main cortical vertexwise analyses
- ./subcortical contains the equivalent analysis on both MAGeT structure volumes and AAL atlas structure volumes from CIVET
- ./laplace contains the repaired cortical thickness analysis switching from tlink thickness to tlaplace
- ./autismSexDifferences is an R package used by most files that had source code for an earlier version of the analysis
- ./figure has code for generating the figures for the paper
