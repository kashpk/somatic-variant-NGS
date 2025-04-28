# somatic-variant-NGS
This pipeline identifies somatic mutations from RNA-Seq data using matched normal samples and supports analysis with both T2T-CHM13 and GRCh38 reference genomes. It involves aligning RNA-Seq reads, preprocessing BAM files, detecting variants with MuTect2, and filtering for high-confidence somatic mutations. Input files include tumor and matched normal FASTQs, reference genomes, and gene annotations, while outputs consist of processed BAM files, somatic variant VCFs, and summary statistics. Proper indexing of reference files is essential, and caution is advised when interpreting RNA-Seq-derived variants due to inherent biases. 

For further details or collaboration, please contact Pragya Kashyap at kashyap.3@iitj.ac.in
