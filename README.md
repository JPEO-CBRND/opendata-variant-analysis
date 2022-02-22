# GISAID Analysis Bioinformatics Pipeline
<p>
Assortment of standalone python utilities that run on the command line for variant analysis. Functionality includes: parsing GISAID data exports and FASTA file creation, genome alignment, variant calling, VCF (Variant Call Format) file merging, conversion to parquet, and variant annotation. These utilities were specifically used for GISAID EpiCoV samples, however some of these utilities can be applied to other variant analysis initiatives.
  
Please note that this repository is also similar to Bioinformatics analysis performed for the Variant Prioritization Pipeline.
</p>

Analysis Workflow:

<img width="836" alt="image" src="https://user-images.githubusercontent.com/99741809/155190687-329c4977-e840-4eb3-bfde-a64abec7b506.png">


Example Analysis Steps using GISAID EpiCov Samples:
1. **FASTA Creation from GISAID data exports**: Run `ncov_to_fasta.py` to obtain individual FASTA sequence files from `ncov-ingest` data exports
2. **Variant Calling**: Run `variant_calling.py` to perform alignment, variant calling based on the reference genome
3. **VCF Merging**: Run `variant_calling.py` to perform vcf merging using bcftools or picard's merge
4. **Variant Annotation**: Run `variant_annotation.py` to perform variant annotation on identified variants

## Setup
It is recommended to run these python utilities in a virtual environment for package organization and avoiding versioning conflicts. 
1. **Install Python** (These scripts were developed using Python3.8)
2. **Install PIP3**
3. **Setup Python Virtual Environment**
4. **Activate Python Virtual Environment**: `source {path to venv}/bin/activate`
5. **Install Python Dependencies (Includes libraries for ML Analysis)**: `pip3 install -r requirements.txt`
6. **Run Analysis Scripts in Virtual Environment**

The following third-party tools are used for facilitating variant analysis and required to be installed locally:
* [ncov-ingest](https://github.com/nextstrain/ncov-ingest) (Optional): Nextstrain's `ncov-ingest` to get data from GISAID EpiCov Database. `ncov-ingest` requires a GISAID API endpoint for exporting data, which can be obtained through reaching out to GISAID support. 
* [minimap2](https://github.com/lh3/minimap2/blob/master/misc/README.md): Minimap2 aligns a sequence given a reference genome. Minimap2 includes a javascript utility called paftools.js which calls variants from assembly to reference genome. 
* [bgzip and tabix](http://www.htslib.org/download/): bgzip and tabix are used during VCF merging for compression of vcf and multivcf files and index file creation.
* [picard](https://broadinstitute.github.io/picard/): Picard MergeVcf software is used for VCF merging to create a multivcf for variant annotation analysis.
* [bcftools](https://samtools.github.io/bcftools/): Bcftools's merge software is used for VCF merging to create a multivcf for variant annotation analysis.
* [snpeff and snpsift](https://pcingola.github.io/SnpEff/download/): SnpEff and SnpSift are used for variant annotation prediction and data filtering respectively.

## Usage 
  ### GISAID Data Ingest
  The `ncov_to_fasta.py` script assists parsing individual FASTA sequences from GISAID data downloaded using Nextstrain's `ncov-ingest` application. The script will require first running ncov-ingest's software for downloading GISAID exports. Please note that access to GISAID's API will be required and can be obtained by contacting GISAID support. To download GISAID data:
  
1. Download ncov-ingest: https://github.com/nextstrain/ncov-ingest
2. Make sure GISAID API endpoint is added as an environment variable in the following format (Example URI: https://www.epicov.org/epi3/3p/hcov-19/export/export.json.bz2:
   `GISAID_API_ENDPOINT={gisaid uri}`
3. Run `fetch-from-gisaid` and enter endpoint credentials (this downloads the data from the api endpoint locally in ndjson format)
   `pipenv run ncov-ingest/bin/fetch-from-gisaid > data/gisaid.ndjson`
4. Run `transform-gisaid` (this gives us the metadata.tsv and multifasta sequences.fasta files)
   `pipenv run ncov-ingest/bin/transform-gisaid data/gisaid.ndjson`

   The following output files will be generated and `sequences.fasta` and `metadata.tsv` can be used as input for `ncov_to_fasta.py` for FASTA creation:
   ```
    ubuntu@ip-10-194-21-241:~/gisaid/bi-pipeline/gisaid-pull/ncov-ingest/data/gisaid$ ls -lhtr
    total 128G
    -rw-rw-r-- 1 ubuntu ubuntu 1.3G Oct 26 21:17 metadata.tsv
    -rw-rw-r-- 1 ubuntu ubuntu 195M Oct 26 21:17 additional_info.tsv
    -rw-rw-r-- 1 ubuntu ubuntu 126G Oct 26 21:40 sequences.fasta
    -rw-rw-r-- 1 ubuntu ubuntu 1.2G Oct 26 21:40 metadata.tsv.raw
    ```
  
  
  #### ncov_to_fasta.py
Parse out individual fasta files from ncov-ingest's metadata.tsv and sequences.fasta output. Creates directories all_hcov19_fasta and ps_hcov19_fasta if user requires samples with known patient status to be analyzed separately. Run `python3 ncov_to_fasta.py -h` to see usage samples, below is an overview of the functionality.

- `ncovtofasta`: Parses out individual FASTA sequence files from input from ncov-ingest's metadata.tsv and sequences.fasta files. Will separate out FASTA sequences and metadata with known patient status if specified.

   Example:
   > python3 scripts/ncov_to_fasta.py ncovtofasta --metadata /home/ubuntu/gisaid/bi-pipeline/gisaid-pull/ncov-ingest/data/gisaid/metadata.tsv --multifasta /home/ubuntu/gisaid/bi-pipeline/gisaid-pull/ncov-ingest/data/gisaid/sequences.fasta --patientstatus y --outputdir /home/ubuntu/gisaid/bi-pipeline/output/gisaid-nov-19-2021-pull/fasta

   Expected output:
   ```
   {outputdir}/all_hcov19_fasta :
        output directory for all fasta sequence files parsed out from ncov-ingest's sequences.fasta file
   (Optional)
    {output_dir}/ps_hcov19_fasta :
        Output directory containing only fasta samples that have known patient status
    {output_dir}/ps_hcov19_fasta/ps_metadata.tsv :
        Output metadata for samples with known patient status
   ```

- `uploadtos3`: Uploads .fa files to AWS S3 bucket specified. Assumes that aws cli is configured locally.


### Variant Calling and VCF Merging
Alignment and variant calling are performed on input FASTA files using the `minimap2` software. We skip processing FASTA samples that are determined to have an ambiguous base percentage > 2% (percentage of 'N' nucleotides in the FASTA sequence), and the samples that are skipped are listed in `skipped_samples.txt`. The `minimap2` command that is run from Rayko et al.'s codebase is shown below. The first command performs long read alignment and gives us CIGAR notation and cs tags (sequence alignment strings). The second command sorts by the reference start coordinate. The third command runs a javascript utility called paftools.js that uses the cs tags for calling variants from the alignment (-L20000 is the minimum alignment length and -f is for the reference sequences).

`{minimap2} -c --cs {ref} {input .fa} 2>/dev/null | sort -k6,6 -k8,8n | {k8} {paftools} call -L20000 -f {ref} - > {output .vcf}`

After variant calling is performed, VCF merging can be performed to create a multivcf of all genomic sample .vcf files for downstream analysis. Both merging with `bcftools` and `picard` can optionally be used. Due to large sets of samples and local memory constraints, the algorithm divides the set of samples into smaller 200-file subsets of samples and first applies merging to those smaller sets. Then merging is re-applied to the subsets to derive a final, single multivcf file of all sample variants. **NOTE: Please make sure the number of open files on the local system is set to MORE THAN the # of VCFs/200 before running merge algorithms. The final merging of larger datasets requires accessing all subset files, therefore this system variable needs to be updated. On linux, this may be setting the open file limit using `ulimit -n {number of files}`**

The following `picard` command is used in the algorithm for merging. The parameters, `USE_JDK_INFLATER` and `USE_JDK_DEFLATER`, are included to avoid core dump error. It should be noted that picard's merge utility includes all occurrences of a variant across all samples.

`java -jar {picard jar} MergeVcfs -I {list of vcfs} -O {output directory} --USE_JDK_DEFLATER true --USE_JDK_INFLATER true`

The following `bcftools` command is used in the algorithm for merging. The parameter `-m none` is included to avoid multiallelic variant reporting. Bcftools will only include a single occurrence of a variant across all samples.

`{bcftools_path} merge -m none --file-list {list of .vcf.gz files} -Oz -o {output multivcf file}`

The following `bgzip` and `tabix` commands are used for file compression and index file creation in the algorithm:

`{bgzip path} -c {.vcf} {out_dir}/bcftools_output_comp_vcfs/{.vcf}.gz`
`{tabix_path} -p vcf {.vcf.gz file}`

#### variant_calling.py
  The script variant_call.py includes functionality for invoking Minimap2 for alignment and variant calling on sample FASTA files, merging variant call format (VCF) files with Bcftools or Picard, and convert multivcf format to parquet format. The output of the `picardmergevcf` and `bcfmergevcf` files will be a file named `final_multivcf.vcf.gz`, which can be used in further variant annotation analysis. Below is an overview of functionality.

- `variantcall`:

  Example:
  >python3 /home/ubuntu/gisaid/bi-pipeline/scripts/variant_calling.py variantcall --ref /home/ubuntu/gisaid/bi-pipeline/reference/NC_045512.2.fa --fastadir /home/ubuntu/gisaid/bi-pipeline/output/gisaid-nov-19-2021-pull/fasta/all_hcov19_fasta --outdir /variant-data/gisaid-nov-19-2021-pull/all_hcov19_vcf --minimap2 /home/ubuntu/gisaid/tools/minimap2-2.22_x64-linux/minimap2 --k8 /home/ubuntu/gisaid/tools/minimap2-2.22_x64-linux/k8 --paftools /home/ubuntu/gisaid/tools/minimap2-2.22_x64-linux/paftools.js
  
  Expected Output:
  ```
  {out_dir}/*.vcf :
    All vcf files for input FASTA samples
  {out_dir}/skipped_samples.txt :
    Tab-delimitted file listing out samples that were skipped due to ambiguous base percentage > 2%
  ```

- `picardmergevcf`:
   
   Example:
   >python3 variant_calling.py picardmergevcf --picard /Users/roshnaagarwal/picard/build/libs/picard.jar --vcfdir /Users/roshnaagarwal/Documents/JPEO/GISAID-workflow/vcf_output --outdir /Users/roshnaagarwal/Documents/JPEO/GISAID-workflow/vcf_output

   Expected Output
   ```
   {out_dir}/picard_output_vcf_lists :
        Intermediary directory containing list(s) of vcfs (all vcfs, subsets, merged)
   {out_dir}/picard_output_merged_vcfs :
        Directory containing .vcf.gz compressed files and .tbi index files of all merged vcfs and merged subsets
   {out_dir}/picard_output_merged_vcfs/final_multivcf.vcf.gz :
        Final compressed merged multivcf file
   ```
   
- `bcfmergevcf`:

   Example:
   > python3 /home/ubuntu/gisaid/bi-pipeline/scripts/variant_calling.py bcfmergevcf --bcf /home/ubuntu/gisaid/tools/bcftools-1.13/bcftools --bgzip /home/ubuntu/gisaid/tools/htslib-1.13/bgzip --tabix /home/ubuntu/gisaid/tools/htslib-1.13/tabix --vcfdir /home/ubuntu/gisaid/bi-pipeline/output/gisaid-oct-26-2021-pull/vcf/all_hcov_19_vcf --outdir /home/ubuntu/gisaid/bi-pipeline/output/gisaid-oct-26-2021-pull/vcf/all_hcov_19_vcf

   Expected Output:
   ```
    {out_dir}/bcftools_output_comp_vcfs :
        Directory containing compressed .vcf files and index files
    {out_dir}/bcftools_output_vcf_lists :
        Directory containing list(s) of vcfs (all vcfs, subsets, merged)
    {out_dir}/bcftools_output_merged_vcfs :
        Directory containing .vcf.gz compressed files and .tbi index files of all merged vcfs and merged subsets
    {out_dir}/bcftools_output_merged_vcfs/final_multivcf.vcf.gz :
        Final compressed merged multivcf file
   ```
  - `vcftoparquet`: Convert vcf file to parquet format
  - `removevcfdup`: Remove duplicate variants in picard output
  - `uploadtos3`: Upload vcf files to s3

### Variant Annotation
Variant Annotation is performed using SnpEff and filtering annotation data into a tab-delimited format is performed using SnpSift.
The following SnpEff command is used within the algorithm to perform variant annotation on the compressed multivcf obtained from running vcf merging.

`java -Xmx8g -jar {SnpEff jar} NC_045512.2 {compressed multivcf} > {snpeff output vcf}`

The following SnpSift command is run within the algorithm to filter SnpEff annotation results and extract the following fields into a tab-delimited format: CHROM, POS, REF, ALT, QUAL, FILTER, QNAME, ANN[*].IMPACT (this will have values of HIGH, MODERATE, LOW, MODIFIER), "EFF[*].FUNCLASS", "ANN[*].EFFECT" (Effect in Sequence ontology terms e.g. 'missense_variant', 'synonymous_variant', 'stop_gained', etc.), "ANN[*].GENE" (Gene name, e.g. 'PSD3'), "ANN[*].CODON" (three nucleotide codon), "ANN[*].AA" (Amino acid), "ANN[*].AA_POS" (the position of amino acid in genome).

`cat {snpeff output vcf file} | {vcfEffOnePerLine.pl script} | java -jar {SnpSift jar} extractFields - CHROM POS REF ALT QUAL FILTER QNAME \"EFF[*].IMPACT\" \"EFF[*].FUNCLASS\" \"EFF[*].EFFECT\" \"EFF[*].GENE\" \"EFF[*].CODON\" \"EFF[*].AA\" \"ANN[*].AA_POS\"  > {snpsift output file}`

 - `snpeff`:
 
    Example:
    >python3 /home/ubuntu/gisaid/bi-pipeline/scripts/variant_annotation.py snpeff --jar /home/ubuntu/gisaid/tools/snpEff/snpEff.jar --vcf /home/ubuntu/gisaid/bi-pipeline/output/gisaid-nov-19-2021-pull/vcf/ps_hcov19_vcf/bcftools_output_merged_vcfs/final_multivcf.vcf.gz --outdir /home/ubuntu/gisaid/bi-pipeline/output/gisaid-nov-19-2021-pull
    
    Expected Output:
    ```
    {out_dir}/snpeff-annotation-{timestamp}.vcf :
        Output SnpEff annotation .vcf file
    {out_dir}/snpEff_summary.html :
        Default SnpEff output summary statistics file, default output from SnpEff.
    {out_dir}/snpEff_genes.txt :
        Default SnpEff output tab-separated file having counts of number of variants affecting each transcript and gene
    ```
    
 - `snpsift`:
 
    Example:
    >python3 /home/ubuntu/gisaid/bi-pipeline/scripts/variant_annotation.py snpsift --jar /home/ubuntu/gisaid/tools/snpEff/SnpSift.jar  --vcf /home/ubuntu/gisaid/bi-pipeline/output/gisaid-nov-19-2021-pull/snpeff-annotation-11242021220504.vcf --effscript /home/ubuntu/gisaid/tools/snpEff/scripts/vcfEffOnePerLine.pl --outdir /home/ubuntu/gisaid/bi-pipeline/output/gisaid-nov-19-2021-pull
    
    Expected Output:
    ```
    {out_dir}/snpsift_output_{timestamp}.txt
    ```
 
## Future Enhancements

Current timing on the Bioinformatics workflow steps are listed below. This is based on running the scripts on a single, AWS EC2 c5.9xlarge instance.

<img width="722" alt="image" src="https://user-images.githubusercontent.com/99741809/155190793-ab0e36bc-b76e-4669-80b7-ec1a0f874e8a.png">

* **Pyspark**: Integrating PySpark may introduce scalability and distribution across a larger compute cluster for larger sample datasets.
* **Additional Analysis Tools**: The variant analysis tools (Minimap2, Picard, Bcftools, SnpEff, and SnpSift) were based off of the analysis pipeline used by Rayko et al. in the references below. Additional tools for variant analysis can be included for more comprehensive verification of variant results. 


## References
1. Rayko, Mikhail, and Aleksey Komissarov. "Quality Control Of Low-Frequency Variants In SARS-Cov-2 Genomes". 2020. Cold Spring Harbor Laboratory, doi:10.1101/2020.04.26.062422. Accessed 29 Nov 2021.
