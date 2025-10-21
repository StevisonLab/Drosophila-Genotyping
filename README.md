# Drosophila Genotyping
Scripts and commands detailed here used for validation of DGRP (Drosophila Genetic Reference Panel) and DPSE (*Drosophila pseudoobscura*) lines. Both scripts are very similar to each other and contain extensive comments on usage of commands within.

Both scripts utilize specific VCF and reference FASTA genomes. For DGRP the VCF file can be downloaded [here]( https://zenodo.org/records/155396/files/dgrp2.vcf.gz) directly from DGRP. The reference FASTA file can be downloaded directly from UCSC Genome Browser [here](https://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/dm3.fa.gz). Both links will download the file immediately upon clicking. 

For DPSE the VCF file was produced in-house by Madison Watkins. The page for the reference FASTA file can be found on [NCBI](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_009870125.1/). This will not immediately download the file, though NCBI provides the following curl command to download a folder containing the file on that page:

`curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_009870125.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_009870125.1.zip" -H "Accept: application/zip"`
