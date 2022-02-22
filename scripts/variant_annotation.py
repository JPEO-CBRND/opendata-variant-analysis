#!/usr/bin/python3

"""Apply variant annotation to multivcf or vcf file

This script is written for Variant Analysis efforts for the JPEO Wastewater Surveillance Project.
This script uses SnpEff for variant annotation and Snpsift for filtering resulting annotation files.
This script uses Python3 and requires that `argparse` and `subprocess` be installed within the Python
environment you are running this script in. This script also requires local installations of SnpEff.

Run `python3 variant_annotation.py -h` for usage information

This file contains the following functions:
    * fasta_to_vcf - Uses Minimap2 to run alignment and variant calling on individual fasta files and outputs vcf files
    * main - the main function of the script
"""

import sys, os
import argparse
import time
import subprocess
from Bio import SeqIO

def snpeff_annotation(jar, vcf, out_dir):
    """Runs Snpeff for variant annotation given a vcf file

    SnpEff produces a vcf file as output with annotations, snpeff_summary.html file, and a snpeff_genes.txt file.

    Parameters
    ----------
    jar : str
        snpeff jar path
    vcf : str
        input vcf/multivcf for annotation (.vcf, .vcf.gz format)
    out_dir :
        output directory path for annotation files

    Output
    ------
    {out_dir}/snpeff-annotation-{timestamp}.vcf :
        Output SnpEff annotation .vcf file
    {out_dir}/snpEff_summary.html :
        Default SnpEff output summary statistics file, default output from SnpEff.
    {out_dir}/snpEff_genes.txt :
        Default SnpEff output tab-separated file having counts of number of variants affecting each transcript and gene
    """

    try:
        timestamp = time.strftime("%m%d%Y%H%M%S")
        start = time.time()

        snpeff_cmd = "java -Xmx8g -jar {} NC_045512.2 {} > {}".format(jar, vcf, os.path.join(out_dir, "snpeff-annotation-{}.vcf".format(timestamp)))
        os.system(snpeff_cmd)
        end = time.time()
        print("Variant annotation with snpeff took: {} minutes".format((end - start)/60))
    except Exception as e:
        print(e)
        sys.exit(1)

def snpsift_annotation(jar, vcf, oneperline_script, out_dir):
    """Runs Snpsift extractFields utility for filtering against Snpeff annotation files

    SnpSift's extractFields utility allows us to create a tab-delimited format of vcf annotation output from SnpEff.
    We can use the vcfEffOnePerLine.pl to list one effect per line

    Parameters
    ----------
    snpsnift_jar : str
        snpeff jar path
    vcf : str
        input vcf/multivcf for annotation (.vcf, .vcf.gz format)
    oneperline_script :
        path to vcfEffOnePerLine.pl script from snpeff
    out_dir :
        output directory path for annotation files

    Output
    ------
    {out_dir}/snpsift_output_{timestamp}.txt
    """

    try:
        timestamp = time.strftime("%m%d%Y%H%M%S")
        start = time.time()

        snpsift_cmd = "cat {} | {} | java -jar {} extractFields - CHROM POS REF ALT QUAL FILTER QNAME \"EFF[*].IMPACT\" \"EFF[*].FUNCLASS\" \"EFF[*].EFFECT\" \"EFF[*].GENE\" \"EFF[*].CODON\" \"EFF[*].AA\" \"ANN[*].AA_POS\"  > {}"\
            .format(vcf, oneperline_script, jar, os.path.join(out_dir, "snpsift_output_{}.txt".format(timestamp)))
        os.system(snpsift_cmd)
        end = time.time()
        print("Filtering with SnpSift took: {} minutes".format((end - start)/60))
    except Exception as e:
        print(e)
        sys.exit(1)

def main():
    args_parser = argparse.ArgumentParser(prog='python3 variant_annotation.py',
                                                  description='Run variant annotation using snpeff or snpsift ')
    sub_parser = args_parser.add_subparsers(dest='command')

    # Add the arguments
    va_parser = sub_parser.add_parser('snpeff', help='run variant annotation with snpeff')
    required_args = va_parser.add_argument_group('required arguments')
    required_args.add_argument('--jar', required=True, help='Path to snpeff jar')
    required_args.add_argument('--vcf', required=True, help='Path to vcf file (.vcf or .vcf.gz)')
    required_args.add_argument('--outdir', required=True, help='Path to output directory for annotation files')

    va_parser = sub_parser.add_parser('snpsift', help='run filtering on variant annotation results with snpsift')
    required_args = va_parser.add_argument_group('required arguments')
    required_args.add_argument('--jar', required=True, help='Path to snpeff jar')
    required_args.add_argument('--vcf', required=True, help='Path to vcf file (.vcf or .vcf.gz)')
    required_args.add_argument('--effscript', required=True, help='Path to SnpEff\'s vcfEffOnePerLine.pl')
    required_args.add_argument('--outdir', required=True, help='Path to output directory for annotation files')


    # Parse args and invoke workflows
    args = args_parser.parse_args()
    if args.command == 'snpeff':
        if args.jar and args.vcf and args.outdir:
            snpeff_annotation(jar=args.jar, vcf=args.vcf, out_dir=args.outdir)
    elif args.command == 'snpsift':
        if args.jar and args.vcf and args.effscript and args.outdir:
            snpsift_annotation(jar=args.jar, vcf=args.vcf, oneperline_script=args.effscript, out_dir=args.outdir)


if __name__ == "__main__":
    main()