# atacCNV: Identifying copy number variations from single-cell ATAC-seq data
epiAneufinder is an algorithm used for calling Copy Number Variations (CNVs) from single-cell ATAC (scATAC) data. Single-cell open chromatin profiling via the single-cell Assay for Transposase-Accessible Chromatin using sequencing (scATAC-seq) assay has become a mainstream measurement of open chromatin in single-cells. epiAneufinder exploits the read count information from scATAC-seq data to extract genome-wide copy number variations (CNVs) for each individual cell. epiAneufinder allows the addition of single-cell CNV information to scATAC-seq data, without the need of additional experiments, unlocking a layer of genomic variation which is otherwise unexplored.

### Description
The algorithm works in three steps:

Data preparation (binning, GC correction, removal of blacklisted regions)
Genome segmentation based on maximum Anderson-Darling distance
Gain/loss assignments


### Authors

Akshaya Ramakrishnan (akshaya4r@gmail.com)

Aikaterini Symeonidi (aikaterini.symeonidi@helmholtz-muenchen.de)

Patrick Hanel (patrick.hanel@helmholtz-muenchen.de)

Michael Schubert

Maria Colomé-Tatché (maria.colome@helmholtz-muenchen.de)

