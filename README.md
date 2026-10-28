# joint-qtl-mapping-pleiotropy-yeast
Github repo containing all code associated with the manuscript "Joint QTL mapping reveals patterns of pleiotropy across structured environments in budding yeast".

The Fitness inference folder contains all the scripts required to infer segregant fitness (jointly across replicate sets of flasks with the same environment/media).
Fitness inference is a multi-step process (eg: read filtering, parsing, error-correction, etc) and the large size of each individual input file and the large number
of such input files necessitates running the fitness inference scripts on a High Performance Computing (HPC) cluster in parallel. The batch scripts to do this are all
located in the 'batch' directory. The 'data' directory contains various txt and csv files that are required to run these fitness inference scripts (or are output these
scripts). 

After fitness inference, the next step in our pipeline is QTL inference. The scripts we used for QTL inference and the associated documentation can all be found in the
following github repo:
https://github.com/efenton32/qtl_tools

Following QTL inference, we annotate each QTL (and each SNP in the cross) features such as gene name, mutation type, gene essentiality, gradient-orderedness of effect sizes,
U-shapedness of effect sizes, etc. The scripts for this purpose can be found in the 'QTL annotation' folder. 
