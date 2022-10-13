## About

Nucleosome peaks from healthy controls were downloaded from: Snyder et al (2016): https://doi.org/10.1016/j.cell.2015.11.050

Link to data repository: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71378

The scripts in the runners folder are used to find the following:
- nucleosome_peaks_167.sh
  - The distance to the closest peak of all 167bp fragments.
- nucleosome_peaks_genome.sh
  - The distance to the closest peak of all fragments genome-wide.
- nucleosome_peaks_length.sh
  - The length of fragments within 500bp of a nucleosome peak.
- nucleosome_peaks_site.sh
  - The distance to the closest peak of all fragments within regions specified in a bed file.
