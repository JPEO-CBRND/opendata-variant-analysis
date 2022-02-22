# Variant Prioritization - Bioinformatics Analysis
<p>
The Jupyter Notebook in this repository performs Bioinformatics variant analysis on NCBI Samples that are related to SARS-CoV-2 Wastewater surveillance and identifies potential variants of interest in these samples. The notebook includes data analytics on search results and samples, variant results, and annotation results. Code from GISAID Analysis python utilities are utilized in this notebook for variant analysis.
</p>

Analysis Workflow:
![image](https://media.github.boozallencsn.com/user/10697/files/3c5e0253-1e6a-4e66-9877-0737a9a4eb56)

The following third-party tools are used for facilitating variant analysis and required to be installed locally:
* [ncov-ingest](https://github.com/nextstrain/ncov-ingest) (Optional): Nextstrain's `ncov-ingest` to get data from GISAID EpiCov Database. `ncov-ingest` requires a GISAID API endpoint for exporting data, which can be obtained through reaching out to GISAID support. 
* [minimap2](https://github.com/lh3/minimap2/blob/master/misc/README.md): Minimap2 aligns a sequence given a reference genome. Minimap2 includes a javascript utility called paftools.js which calls variants from assembly to reference genome. 
* [bgzip and tabix](http://www.htslib.org/download/): bgzip and tabix are used during VCF merging for compression of vcf and multivcf files and index file creation.
* [bcftools](https://samtools.github.io/bcftools/): Bcftools's merge software is used for VCF merging to create a multivcf for variant annotation analysis.
* [snpeff and snpsift](https://pcingola.github.io/SnpEff/download/): SnpEff and SnpSift are used for variant annotation prediction and data filtering respectively.

The steps for analysis in the Notebook are:
1. **Searching for SARS-CoV-2 Wastewater Samples**: 
2. **Variant Calling**: Run variant calling using Minimap2 on SARS-CoV-2 Wastewater FASTA sequences identified in NCBI searches
3. **VCF Merging**: Run merging on VCF files using bcftools
4. **Variant Annotation**: Run variant annotation and filtering using SnpEff and SnpSift on merged vcf
5. **Validation Against Clinical Variants**: Compare sewage variants in NCBI to clinical variants in NCBI

## Setup
It is recommended to run these python utilities in a virtual environment for package organization and avoiding versioning conflicts. Please make sure to replace local paths for installations and files when running the Notebook code.
1. **Install Python** (The Notebook Code was developed using Python3.8)
2. **Install PIP3**
3. **Install Jupyter and Setup Python Interpreter as Needed**
4. **Install Libraries listed in Notebook**

## Future Enhancements
* **Including Additional Databases**: This analysis focuses on NCBI Samples searches in the Entrez Nucleotide database. Including other databases such as the SRA database may include more SARS-CoV-2 wastewater-surveillance related samples. Also including samples that may be in other open-data sources such as GISAID may broaden identification of sewage-specific variants.
* **Additional Analysis Tools**: The variant analysis tools (Minimap2, Picard, Bcftools, SnpEff, and SnpSift) were based off of the analysis pipeline used by Rayko et al. in the references below. Additional tools may be included to facilitate/confirm analysis results
* **Comparing Sewage Annotation Results to Other Datasets**: Merging Sewage Annotation results with other datasets may derive more insights on SARS-CoV-2 variants and their affects on virulence, drug resistance, etc. For example, The National Center of Advanced Translational Sciences (NCATS) contains a repository of open-source data on Viral SARS-CoV-2 Mutations and relation to drug resistance. The variants in the sewage sample subsets could potentially be related back to this data set by amino acid mutation.


