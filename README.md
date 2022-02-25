<h2>Air Force COVID-19 Wastewater Surveillance:
Variant Prioritization Pipeline Summary</h2>

**BLUF:** A Variant Prioritization Pipeline was developed to help identify and track publicly deposited data and available literature from user-supplied topics of interest in order to correlate them with SARS-CoV-2 viral genomic variants. The pipeline will help users collate open-source literature and data on the evolving COVID-19 pandemic in a faster, easier, and automated fashion for various downstream purposes (e.g., SARS-CoV-2 assay design and development).  

**Background:** Since the onset of the COVID-19 pandemic, SARS-CoV-2 genomic variants are typically characterized by their impacts on patients, including disease transmissibility, disease severity, immune system evasion, and vaccine efficacy (reviewed in Tao et al., 2021). As SARS-CoV-2 variants emerge more rapidly than epidemiological studies can be conducted, correlates of protection against SARS-CoV-2 variants have become of great interest. A comprehensive understanding of these laboratory proxy measurements derives from epidemiological research paired with animal models and in vitro studies, including computational, biochemical, and protein structure assays. These scientific findings are often disseminated in the form of pre-print and peer-reviewed publications.

**Results:** We present a Variant Prioritization Pipeline to help track and potentially flag emerging SARS-CoV-2 viral genomic variants. The pipeline is composed of two workflows and currently focuses on wastewater surveillance. The first workflow queries NCBI Genbank, the major domestic source for virus nucleotide sequences, for SARS-CoV-2 sequences collected from wastewater or environmental sources. The workflow uses Minimap2, Bcftools, and SnpEff/SnpSift to call, merge, and annotate genomic variants.    The second workflow accepts search terms from a user to automatically retrieve and annotate relevant publications. The workflow uses NCBI tools, Pubtator, and custom Python scripting to display search terms mapped to SARS-CoV-2 viral variants. 

 ![image](https://user-images.githubusercontent.com/99741809/155772174-4a8f4c7c-1a54-4f56-9a6f-cbc17f1d2716.png)
 <img width="470" alt="image" src="https://user-images.githubusercontent.com/99741809/155772264-2f74b138-3f9a-4862-b9a0-17cde73e0a4d.png">

 
**References**

Tao, K., Tzou, P.L., Nouhin, J. et al. The biological and clinical significance of emerging SARS-CoV-2 variants. Nat Rev Genet 22, 757–773 (2021). https://doi.org/10.1038/s41576-021-00408-x
