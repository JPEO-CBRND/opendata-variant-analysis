#!/usr/bin/python3

"""
This script is written for GISAID data analysis efforts for the JPEO Wastewater Surveillance Project.
Expected input is Nextstrain's ncov-ingest output files, metadata.tsv and sequences.fasta and converts to individual
FASTA files and isolates samples with known patient status if desired. This script will also upload output fasta files to s3,
given the s3 bucket and assumes that the local environment already has aws cli configured.

This script uses Python3 and requires that `pandas`,`argparse`,`shutil`, and `copyfile` be installed within the Python
environment you are running this script in.

Run `python3 ncov_to_fasta.py -h` for usage information.

"""

import sys, os
import argparse
import time
import pandas as pd
import shutil
from shutil import copyfile

def parse_ncov_ingest(metadata, multifasta, outputdir, pat_status):
    """Parses all fasta sequences from Nextstrain's ncov-ingest and separates out samples with patient status if desired for analysis.
    1. Parse out individual fasta sequences from sequences.fasta output
    2. If samples with known patient status is desired, invokes isolate_patient_fastas() function to isolate samples with known patient status

    Parameters
    ----------
    metadata : str
        local metadata.tsv file path from ncov-ingest
    multifasta : str
        local sequences.fasta file path from ncov-ingest
    outputdir :
        output directory path for individual fasta files
    pat_status :
        Either `p` for patient status only sequences or `a` for all sequences

    Output
    ------
    {outputdir}/all_hcov19_fasta :
        output directory for all fasta sequence files parsed out from ncov-ingest's sequences.fasta file
    """

    # Check input and output paths
    if not os.path.isfile(metadata):
        print("The file path specified does not exist: {}".format(metadata))
        sys.exit()
    if not os.path.isfile(multifasta):
        print("The file path specified does not exist: {}".format(multifasta))
        sys.exit()
    if not os.path.isdir(outputdir):
        print("The directory path specified does not exist: {}".format(outputdir))
        sys.exit()

    try:
        timestamp = time.strftime("%m%d%Y%H%M%S")

        # Create new directory for all fasta files
        output_path = os.path.join(outputdir, "all_hcov19_fasta")

        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        os.mkdir(output_path)

        # Open multifasta and parse out sequences
        # Parse out all fasta sequences
        with open(multifasta) as f_input:
            for line in f_input:
                if line.startswith('>'):
                    acc = line.split(">")[1].strip()
                    print("Creating .fa for: {}".format(acc))
                    f = open("{}/{}.fa".format(output_path, acc), "w")
                    f.write(line)
                elif line.startswith(('A', 'T', 'G', 'C')):
                    f.write(line)
                    f.close()

        # If patient status samples desired, separate these out
        if (pat_status == 'y'):
            isolate_patient_fastas(outputdir, metadata, output_path)

    except Exception as e:
        print(e)
        sys.exit(1)

def isolate_patient_fastas(output_dir, metadata, fasta_dir):
    '''
    If fasta output with patient status is required, create a separate directory of patient status only samples and copy
    FASTA files to that directory for consumption. Use ncov-ingest's metadata.tsv file as reference.

    Parameters
    ----------
    output_dir : str
        Output directory for new directory to save fasta samples with known patient status
    metadata : str
        Path to ncov-ingest's metadata.tsv file
    fasta_dir : str
        FASTA directory with parsed FASTA files

    Output
    ------
    {output_dir}/ps_hcov19_fasta :
        Output directory containing only fasta samples that have known patient status
    {output_dir}/ps_hcov19_fasta/ps_metadata.tsv :
        Output metadata for samples with known patient status
    '''

    try:
        # Create directory for FASTAs with known patient status
        ps_output_path = os.path.join(output_dir, "ps_hcov19_fasta")
        if os.path.exists(ps_output_path):
            shutil.rmtree(ps_output_path)
        os.mkdir(ps_output_path)

        # Load metadata into dataframe
        # Filter to known patient status
        # Save patient status metadata output to new file
        metadata_file_path = os.path.join(ps_output_path, "ps_metadata.tsv")
        metadata_df = pd.read_csv(metadata, sep='\t')
        metadata_with_ps = metadata_df[(metadata_df['covv_patient_status'] != "unknown") & (metadata_df['covv_patient_status'].notnull())]
        metadata_with_ps.to_csv(metadata_file_path, sep="\t", index=False)
        print("Patient Status Metadata saved at: {}".format(metadata_file_path))

        # Copy FASTAs to ps_hcov19_fasta
        for index,row in metadata_with_ps.iterrows():
            src = os.path.join(fasta_dir, "{}.fa".format(row['gisaid_epi_isl']))
            dst = os.path.join(ps_output_path, "{}.fa".format(row['gisaid_epi_isl']))

            if os.path.isfile(src):
                copyfile(src, dst)

        print("Total samples with known patient status: {}".format(len(metadata_with_ps)))

    except Exception as e:
        print(e)
        sys.exit(1)

def upload_to_s3(fastadir, bucket):
    """ Uploads fasta files to s3

    Parameters
    ----------
    fastadir : str
        directory path for fasta files
    bucket : str
        s3 bucket to upload fasta files to
    """
    try:
        if not os.path.isdir(fastadir):
            print("The directory path specified does not exist: {}".format(fastadir))
            sys.exit()

        print("Uploading FASTA from {} to {}".format(fastadir, bucket))

        for f in os.listdir(fastadir):
            if f.endswith(".fa"):
                print("Uploading {}".format(os.path.join(fastadir,f)))
                os.system("aws s3 cp {} {}".format(os.path.join(fastadir, f), bucket))

    except Exception as e:
        print(e)
        sys.exit(1)

def main():
    args_parser = argparse.ArgumentParser(prog='ncov_to_fasta',
                                          description='Parse out individual fasta files from ncov-ingest metadata.tsv and .sequences output. Isolate samples with known patient status if desired.')

    sub_parser = args_parser.add_subparsers(dest='command')

    # ncovtofasta arguments
    ncov_args_parser = sub_parser.add_parser('ncovtofasta', help='Parse out individual fasta files from ncov-ingest metadata.tsv and .sequences output')
    required_args = ncov_args_parser.add_argument_group('ncov to fasta arguments')
    required_args.add_argument('--metadata', required=True, help='Input metadata.tsv path from ncov_ingest')
    required_args.add_argument('--multifasta', required=True, help='Input sequences.fasta path from ncov_ingest')
    required_args.add_argument('--patientstatus', required=True, help='Identify samples with known patient status? Either `y` for yes or `n` no. ')
    required_args.add_argument('--outputdir', required=True, help='The path to the output directory for fasta files')

    # uploadtos3 arguments
    s3_args_parser = sub_parser.add_parser('uploadtos3', help='Upload FASTA output files to s3 bucket')
    s3_args = s3_args_parser.add_argument_group('s3 upload arguments')
    s3_args.add_argument('--fastadir', required=True, help='Directory with fasta files')
    s3_args.add_argument('--s3bucket', required=True, help='s3 Bucket to upload output fasta files to')

    # Invoke utilities
    args = args_parser.parse_args()
    if args.command == 'ncovtofasta':
        if args.metadata and args.multifasta and args.patientstatus and args.outputdir:
            parse_ncov_ingest(metadata=args.metadata, multifasta=args.multifasta, outputdir=args.outputdir, pat_status=args.patientstatus)
    elif args.command == 'uploadtos3':
        if args.fastadir and args.s3bucket:
            upload_to_s3(fastadir=args.fastadir, bucket=args.s3bucket)

if __name__ == "__main__":
    main()
