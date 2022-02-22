#!/usr/bin/python3

"""Set of python functions to run alignment, variant calling, vcf merges, vcf to parquet conversion

This script is written for Variant Analysis efforts for the JPEO Wastewater Surveillance Project. It utilizes Minimap2 for
alignment and variant calling (paftools.js in Minimap2). Merging can also be performed on output vcfs using Picard or Bcftools.
Vcf or multivcf files can also be converted to parquet format.

This script uses Python3 and requires that `argparse`, `subprocess`, `shutil`, `glob`, `Bio`, `SeqIO`, `pyarrow`,
 `fastparquet`, `allel`, `sys`, `hashlib`, and `os` be installed within the Python environment you are running this script in.
 This script also requires local software installations of Minimap2, picard, bcftools, bgzip, and tabix.

Run `python3 variant_calling.py -h` for usage information

This file contains the following functions:
    * fasta_to_vcf - Uses Minimap2 to run alignment and variant calling on individual fasta files and outputs vcf files
    * merge_vcfs_picard - Merges vcfs using picard's MergeVcfs tool
    * merge_vcfs_bcf - Merges vcfs using bcftools merge tool
    * vcf_to_parquet - Converts vcf/multivcf to parquet format
    * main - The main function of the script
"""

import sys, os
import argparse
import time
import subprocess
import glob
import shutil
import pandas as pd
import pyarrow
import fastparquet
import allel
from Bio import SeqIO
import hashlib

def fasta_to_vcf(ref, fasta_dir, out_dir, minimap2, k8, paftools):
    """
    Runs Minimap2 alignment and variant calling on input fasta files and outputs resulting .vcf files.

    Pipeline executed:
    The minimap2 command that is run from Rayko et al.'s codebase is shown below. The first command performs long read alignment and gives us CIGAR notation and cs tags (sequence alignment strings).
    The second command sorts by the reference start coordinate.
    The third command runs a javascript utility called paftools.js that uses the cs tags for calling variants from the alignment (-L20000 is the minimum alignment length and -f is for the reference sequences)

    `{minimap2} -c --cs {ref} {input .fa} 2>/dev/null | sort -k6,6 -k8,8n | {k8} {paftools} call -L20000 -f {ref} - > {output .vcf}`

    If a FASTA sample has greater than 2% ambiguous 'N' bases, the sample is skipped and not processed.

    Parameters
    ----------
    ref : str
        reference genome (.fa) path
    fasta_dir : str
        input directory containing fasta (.fa) files
    out_dir : str
        output directory path for .vcf files
    minimap2 : str
        Minimap2 script path
    k8 : str
        K8 javascript shell path (this is included in Minimap2 installation)
    paftools : str
        paftools path (this is a utility provided by Minimap2)

    Output
    ------
    {out_dir}/*.vcf :
        All vcf files for input FASTA samples
    {out_dir}/skipped_samples.txt :
        Tab-delimitted file listing out samples that were skipped due to ambiguous base percentage > 2%
    """

    try:
        timestamp = time.strftime("%m%d%Y%H%M%S")
        start = time.time()

        skipped_samples = open(os.path.join(out_dir, "skipped_samples.txt"), 'w')
        skipped_samples.write("SAMPLE\tSEQ_LENGTH\tN_COUNT\tN_PERCENTAGE\n")

        for f in os.listdir(fasta_dir):
            if f.endswith(".fa"):

                # Skip samples with ambiguous 'N' bases percentage > 2%
                record = SeqIO.read(os.path.join(fasta_dir, f), "fasta")
                n_count = record.seq.count("N")
                seq_len = len(record.seq)

                # Error check
                if seq_len < 1:
                    continue

                n_percentage = n_count / seq_len
                if n_percentage > 0.02:
                    skipped_samples.write("{}\t{}\t{}\t{:.2f}\n".format(f, seq_len, n_count, n_percentage))
                    continue

                # Process sample
                print('Now working on: {}'.format(f))
                minimap2_cmd = subprocess.Popen(
                    ['{}'.format(minimap2), '-c', '--cs', '{}'.format(ref), '{}'.format(os.path.join(fasta_dir, f)), ]
                    , stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, )
                minimap2_cmd.wait()

                sort_cmd = subprocess.Popen(['sort', '-k6,6', '-k8,8n'], stdin=minimap2_cmd.stdout,
                                            stdout=subprocess.PIPE, )
                sort_cmd.wait()

                f_out = open(os.path.join(out_dir, '{}.vcf'.format(os.path.splitext(os.path.basename(f))[0])), "w")
                paftools_cmd = subprocess.Popen(
                    ['{}'.format(k8), '{}'.format(paftools), 'call', '-L20000', '-f', '{}'.format(ref),
                     '-', '>'], stdin=sort_cmd.stdout, stdout=f_out)
                paftools_cmd.wait()
                f_out.close()
            else:
                continue

        skipped_samples.close()
        end = time.time()
        print("Variant calling took: {} minutes".format((end - start) / 60))
    except Exception as e:
        print(e)
        sys.exit(1)


def merge_vcfs_picard(picard, vcf_dir, out_dir):
    """
    Use picard to merge individual .vcf files into a single multivcf file. Merging is first applied to 200-file vcf subsets
    and then merging is re-applied to these subsets to derive a single multifasta file. This optimization was included to
    remediate local memory issues.

    NOTE: Please make sure the number of open files on the local system is set to MORE THAN the # of VCFs/200 before running merge algorithms.
    The final merging of larger datasets requires accessing all subset files, therefore this system variable needs to be updated.
    On linux, this may be setting the open file limit using `ulimit -n {number of files}`

    Command executed: `java -jar {path to picard installation}/picard.jar MergeVcfs -I {vcf list} -O {path to output .gz multivcf}`

    Parameters
    ----------
    picard : str
        path to picard installation for picard.jar
    bgzip : str
        path to bgzip installation
    vcf_dir : str
        input directory path for .vcf files
    out_dir :
        output directory path for multifasta .gz file

    Output
    ------
    {out_dir}/picard_output_vcf_lists :
        Intermediary directory containing list(s) of vcfs (all vcfs, subsets, merged)
    {out_dir}/picard_output_merged_vcfs :
        Directory containing .vcf.gz compressed files and .tbi index files of all merged vcfs and merged subsets
    {out_dir}/picard_output_merged_vcfs/final_multivcf.vcf.gz :
        Final compressed merged multivcf file
    """

    timestamp = time.strftime("%m%d%Y%H%M%S")
    start = time.time()
    try:

        # Pre-check for ulimit to avoid final merging issues
        # The shell's ulimit should be greater than the value of #_VCF_files/200 (200 is the size of merging subsets)
        vcfs_num = len(os.listdir(vcf_dir))
        ulimit_val = subprocess.check_output("ulimit -n", shell=True)
        ulimit_thres = vcfs_num / 200
        if (int(ulimit_val.strip()) < ulimit_thres):
            print("ERROR: System number of open files (ulimit) needs to be increased to more than {} else merging may fail.".format(ulimit_thres))
            sys.exit(1)

        # 1. Create list(s) of .vcf file to process
        print("1. Creating output directory for vcf lists: picard_output_vcf_lists")
        list_path = os.path.join(out_dir, "picard_output_vcf_lists")
        if os.path.exists(list_path):
            shutil.rmtree(list_path)
        os.mkdir(list_path)

        # Split files into lists of vcf subsets if number of vcfs > 200
        print("2. Splitting files into lists of 200 subsets")
        if len(glob.glob1(vcf_dir, "*.vcf")) <= 200:
            print("Working on list of .vcfs to merge: {}/subset_vcfs".format(list_path))
            os.system("find {} -name '*.vcf' > {}/subset_vcfs".format(vcf_dir, list_path))
        else:
            print("Splitting .vcfs into 200-file subsets to merge: {}/subset_vcfs*".format(list_path))
            print(len(glob.glob1(vcf_dir, "*.vcf")))
            os.system("find {} -name '*.vcf' | split -l 200 - {}/subset_vcfs".format(vcf_dir, list_path))
        sys.exit(0)

        # 2. Apply merging to .vcf files or multi-merging if subsets are created
        print("3. Create output directory for merged vcfs: picard_output_merged_vcfs")
        merge_path = os.path.join(out_dir, "picard_output_merged_vcfs")
        if os.path.exists(merge_path):
            shutil.rmtree(merge_path)
        os.mkdir(merge_path)
        subset_vcfs_list = glob.glob("{}/subset_vcfs*".format(list_path))
        print(list_path)

        if len(subset_vcfs_list) == 1:
            # 2a. Single vcf list was created
            print("Merging single vcf list...")
            multivcf_path = os.path.join(merge_path, "final_multivcf.vcf.gz")
            os.system(
                "java -jar {} MergeVcfs -I {}/subset_vcfs -O {} --USE_JDK_DEFLATER true --USE_JDK_INFLATER true".format(
                    picard, list_path, multivcf_path))
        elif len(subset_vcfs_list) > 1:
            # 2b. Subset lists of vcfs was created
            # Merge each subset into multivcfs
            print("Subsets detected... merging subsets of vcfs")
            for f in os.listdir(list_path):
                multivcf_path = os.path.join(merge_path, "multivcf_{}.vcf.gz".format(os.path.splitext(f)[0]))
                if f.startswith("subset_vcfs"):
                    os.system(
                        "java -jar {} MergeVcfs -I {} -O {} --USE_JDK_DEFLATER true --USE_JDK_INFLATER true".format(
                            picard, os.path.join(list_path, f), multivcf_path))

            # Create new list of merged subsets
            # Picard is creating merged subsets of 0 bytes, skip these subsets
            print("Creating new list of merged subsets...")
            vl_path = os.path.join(list_path, "merged_subsets.list")
            vcf_list = open(vl_path, "w")
            for f in os.listdir(merge_path):
                if f.startswith("multivcf_") and f.endswith(".vcf.gz") and os.path.getsize(
                        os.path.join(merge_path, f)) > 0:
                    vcf_list.write(os.path.join(merge_path, f) + "\n")
            vcf_list.close()
            print("VCF list file added: {}".format(vl_path))

            # Merge subsets into single vcf
            print("Merging subsets into single vcf...")
            final_multivcf = os.path.join(merge_path, "final_multivcf.vcf.gz")
            os.system(
                "java -jar {} MergeVcfs -I {} -O {} --USE_JDK_DEFLATER true --USE_JDK_INFLATER true".format(picard,
                                                                                                            vl_path,
                                                                                                            final_multivcf))
        else:
            print("Error: no vcf file lists were produced. Exiting.")
            sys.exit(1)
    except Exception as e:
        print(e)

    end = time.time()
    print("Merge with Picard MergeVcfs took: {} minutes".format((end - start) / 60))


def merge_vcfs_bcf(bcf, bgzip, tabix, vcf_dir, out_dir):
    """
    Use bcftools to merge individual .vcf files into a single multivcf file. Merging is first applied to 200-file vcf subsets
    and then merging is re-applied to these subsets to derive a single multifasta file. Bcftools requires bgzip compression of
    .vcf files and index file creation before running merge operations. Bcftools merge does not include multiallelic variants with the '-m none' options.

    NOTE: Please make sure the number of open files on the local system is set to MORE THAN the # of VCFs/200 before running merge algorithms.
    The final merging of larger datasets requires accessing all subset files, therefore this system variable needs to be updated.
    On linux, this may be setting the open file limit using `ulimit -n {number of files}`

    Commands executed:
     Compression : `{bgzip path} -c {.vcf} {out_dir}/bcftools_output_comp_vcfs/{.vcf}.gz`
     Index Creation: `{tabix_path} -p vcf {.vcf.gz file}`
     Merging : `{bcftools_path} merge -m none --file-list {list of .vcf.gz files} -Oz -o {output multivcf file}`

    Parameters
    ----------
    bcf : str
        path to bcftools installation
    bgzip : str
        path to bgzip installation
    tabix : str
        path to tabix installation for index files
    vcf_dir : str
        input directory path for .vcf files
    out_dir :
        output directory path for multifasta .gz file

    Output
    ------
    {out_dir}/bcftools_output_comp_vcfs :
        Directory containing compressed .vcf files and index files
    {out_dir}/bcftools_output_vcf_lists :
        Directory containing list(s) of vcfs (all vcfs, subsets, merged)
    {out_dir}/bcftools_output_merged_vcfs :
        Directory containing .vcf.gz compressed files and .tbi index files of all merged vcfs and merged subsets
    {out_dir}/bcftools_output_merged_vcfs/final_multivcf.vcf.gz :
        Final compressed merged multivcf file
    """

    timestamp = time.strftime("%m%d%Y%H%M%S")
    start = time.time()
    try:

        # Pre-check for ulimit to avoid final merging issues
        # The shell's ulimit should be greater than the value of #_VCF_files/200 (200 is the size of merging subsets)
        vcfs_num = len(os.listdir(vcf_dir))
        ulimit_val = subprocess.check_output("ulimit -n", shell=True)
        ulimit_thres = vcfs_num / 200
        if (int(ulimit_val.strip()) < ulimit_thres):
            print("ERROR: System number of open files (ulimit) needs to be increased to more than {} else merging may fail.".format(ulimit_thres))
            sys.exit(1)

        # 1. Apply bgzip compression to vcf files and create index files using tabix
        comp_path = os.path.join(out_dir, "bcftools_output_comp_vcfs")
        if os.path.exists(comp_path):
            shutil.rmtree(comp_path)
        os.mkdir(comp_path)
        print("1. Created directory for bgzip and index files: {}".format(comp_path))

        # Compress vcf files and create index file
        print("2. Applying bgzip compression and created index files...")
        for f in os.listdir(vcf_dir):
            if f.endswith(".vcf"):
                print("-- working on {}".format(f))
                os.system("{} -c {} > {}".format(bgzip, os.path.join(vcf_dir, f), os.path.join(comp_path, f + ".gz")))
                os.system("{} -p vcf {}".format(tabix, os.path.join(comp_path, f + ".gz")))

        # 2. Create list(s) of .vcf file to process
        list_path = os.path.join(out_dir, "bcftools_output_vcf_lists")
        if os.path.exists(list_path):
            shutil.rmtree(list_path)
        os.mkdir(list_path)
        print("3. Created directory for file lists: {}".format(list_path))

        # Split .vcf.gz files into lists of compressed vcf subsets if number of vcfs > 200
        print("4. Splitting files")
        if len(glob.glob1(comp_path, "*.vcf.gz")) <= 200:
            print("Working on list of .vcf.gz files to merge: {}/subset_vcfs".format(list_path))
            os.system("find {} -name '*.vcf.gz' > {}/subset_vcfs".format(comp_path, list_path))
        else:
            print("Splitting .vcf.gz files into 200-file subsets to merge: {}/subset_vcfs*".format(list_path))
            print(vcf_dir)
            os.system("find {} -name '*.vcf.gz' | split -l 200 - {}/subset_vcfs".format(comp_path, list_path))

        # 3. Apply merging to .vcf files or multi-merging if subsets are created
        print("5. Working on merging with bcftools...")
        merge_path = os.path.join(out_dir, "bcftools_output_merged_vcfs")
        if os.path.exists(merge_path):
            shutil.rmtree(merge_path)
        os.mkdir(merge_path)
        subset_vcfs_list = glob.glob("{}/subset_vcfs*".format(list_path))
        print(list_path)

        if len(subset_vcfs_list) == 1:
            # 3a. Single vcf list was created
            print("Merging single vcf list...")
            multivcf_path = os.path.join(merge_path, "final_multivcf.vcf.gz")
            os.system("{} merge -m none --file-list {}/subset_vcfs -Oz -o {}".format(bcf, list_path, multivcf_path))
        elif len(subset_vcfs_list) > 1:
            # 3b. Subset lists of vcfs was created
            # Merge each subset into multivcfs and create new index file
            print("Subsets detected... merging subsets of .vcf.gz files")
            for f in os.listdir(list_path):
                multivcf_path = os.path.join(merge_path, "multivcf_{}.vcf.gz".format(os.path.splitext(f)[0]))
                if f.startswith("subset_vcfs"):
                    os.system("{} merge --force-samples -m none --file-list {} -Oz -o {}".format(bcf,
                                                                                                 os.path.join(list_path,
                                                                                                              f),
                                                                                                 multivcf_path))
                    os.system("{} -p vcf {}".format(tabix, multivcf_path))

            # Create new list of merged subsets
            print("Creating new list of merged subsets...")
            vl_path = os.path.join(list_path, "merged_subsets.list")
            vcf_list = open(vl_path, "w")
            for f in os.listdir(merge_path):
                if f.startswith("multivcf_") and f.endswith(".vcf.gz") and os.path.getsize(
                        os.path.join(merge_path, f)) > 0:
                    vcf_list.write(os.path.join(merge_path, f) + "\n")
            vcf_list.close()
            print("VCF list file added: {}".format(vl_path))

            # Merge subsets into single vcf and create index file
            print("Merging subsets into single vcf...")
            final_multivcf = os.path.join(merge_path, "final_multivcf.vcf.gz")
            os.system("{} merge --force-samples -m none --file-list {} -Oz -o {}".format(bcf, vl_path, final_multivcf))
            os.system("{} -p vcf {}".format(tabix, final_multivcf))
        else:
            print("Error: no vcf file lists were produced. Exiting.")
            sys.exit(1)
    except Exception as e:
        print(e)
        sys.exit(1)

    end = time.time()
    print("Merge with bcftools merge took: {} minutes".format((end - start) / 60))


def vcf_to_parquet(vcf_path, out_dir):
    """
    Use scikit-allel and pandas to convert a multivcf or single vcf file to parquet format

    Parameters
    ----------
    vcf_dir : str
        input directory path for vcf file (.vcf or .vcf.gz)
    out_dir :
        output directory path for multivcf parquet file

    Output
    ------
    {out_dir}/multivcf_parquet.parquet :
        parquet compressed format of multivcf or single vcf file
    """

    try :
        start = time.time()

        print('Using scikit-allel for vcf to dataframe conversion...')
        df = allel.vcf_to_dataframe(vcf_path)
        df.to_parquet(os.path.join(out_dir, 'multivcf_parquet.parquet'))

        df_test = pd.read_parquet(os.path.join(out_dir, 'multivcf_parquet.parquet'), engine='fastparquet')
        print("Verifying total number of variants in parquet file: {}".format(len(df_test.index)))

        end = time.time()
        print("Conversion to parquet took: {} minutes".format((end - start) / 60))
    except Exception as e:
        print(e)

def remove_dup_from_vcf(input_vcf, output_vcf_path, bgzip, tabix):
    """
    Helper function to remove duplicate variants in vcf files generated by picard.
    Ref: https://www.codevscolor.com/python-remove-duplicate-lines-text-file

    Parameters
    ----------
    input_vcf :
    output_vcf_path :

    Output
    ------
    {output_vcf_path}/final_vcf_cleaned.vcf :
        VCF with duplicate variants removed:
    """

    try:
        # Create a set of hash'ed lines in the file
        completed_lines_hash = set()
        output_vcf_name = os.path.join(output_vcf_path, "final_vcf_cleaned.vcf")
        output_file = open(output_vcf_name, "w")

        # Remove duplicate variants using hashvalue
        for line in open(input_vcf, "r"):
            hashValue = hashlib.md5(line.rstrip().encode('utf-8')).hexdigest()
            if hashValue not in completed_lines_hash:
                output_file.write(line)
                completed_lines_hash.add(hashValue)
        output_file.close()

        # Compress final vcf and create index file
        print("{} {}".format(bgzip, output_vcf_name))
        print("{} -p vcf {}".format(tabix, output_vcf_name))
        os.system("{} {}".format(bgzip, output_vcf_name))
        os.system("{} -p vcf {}.gz".format(tabix, output_vcf_name))

    except Exception as e:
        print(e)

def upload_to_s3(vcfdir, bucket):
    """ Uploads vcf files to s3

    Parameters
    ----------
    vcfdir : str
        directory path for fasta files
    bucket : str
        s3 bucket to upload fasta files to
    """
    try:
        if not os.path.isdir(vcfdir):
            print("The directory path specified does not exist: {}".format(vcfdir))
            sys.exit()

        print("Uploading VCFs from {} to {}".format(vcfdir, bucket))

        for f in os.listdir(vcfdir):
            if f.endswith(".vcf"):
                print("Uploading {}".format(f))
                os.system("aws s3 cp {} {}".format(os.path.join(vcfdir, f), bucket))

    except Exception as e:
        print(e)


def main():
    args_parser = argparse.ArgumentParser(prog='python3 variant_calling.py',
                                          description='Run minimap2 for alignment and variant calling on individual fasta files. Also allows users to merge vcfs using picard ')
    sub_parser = args_parser.add_subparsers(dest='command')

    # Add the arguments
    vc_parser = sub_parser.add_parser('variantcall', help='run alignment and variant calling')
    required_args = vc_parser.add_argument_group('variant calling arguments')
    required_args.add_argument('--ref', required=True, help='Path to reference genome')
    required_args.add_argument('--fastadir', required=True, help='Directory path to fasta file(s)')
    required_args.add_argument('--outdir', required=True, help='Path to output directory for vcfs')
    required_args.add_argument('--minimap2', required=True, help='Path to minimap2 script')
    required_args.add_argument('--k8', required=True, help='Path to K8')
    required_args.add_argument('--paftools', required=True, help='Path to paftools.js')

    mergevcf_picard_parser = sub_parser.add_parser('picardmergevcf', help='merge individual vcf files with picard')
    mergevcf_picard_args = mergevcf_picard_parser.add_argument_group('merge vcf arguments')
    mergevcf_picard_args.add_argument('--picard', required=True, help='Path to picard jar')
    mergevcf_picard_args.add_argument('--vcfdir', required=True, help='Directory path to vcf files')
    mergevcf_picard_args.add_argument('--outdir', required=True,
                                      help='Directory path where output file should be saved')

    mergevcf_bcf_parser = sub_parser.add_parser('bcfmergevcf', help='merge individual vcf files with bcf')
    mergevcf_bcf_args = mergevcf_bcf_parser.add_argument_group('merge vcf arguments')
    mergevcf_bcf_args.add_argument('--bcf', required=True, help='Path to bcf installation')
    mergevcf_bcf_args.add_argument('--bgzip', required=True, help='Path to bgzip script')
    mergevcf_bcf_args.add_argument('--tabix', required=True, help='Path to tabix script')
    mergevcf_bcf_args.add_argument('--vcfdir', required=True, help='Directory path to vcf files')
    mergevcf_bcf_args.add_argument('--outdir', required=True, help='Directory path where output file should be saved')

    vcftoparquet_parser = sub_parser.add_parser('vcftoparquet', help='convert vcf to parquet')
    vcftoparquet_bcf_args = vcftoparquet_parser.add_argument_group('convert vcf file to parquet')
    vcftoparquet_bcf_args.add_argument('--vcf', required=True, help='Path to vcf file')
    vcftoparquet_bcf_args.add_argument('--outdir', required=True, help='Path to output directory')

    vcfdup_parser = sub_parser.add_parser('removevcfdup', help='remove duplicate variants in vcf file')
    vcfdup_args = vcfdup_parser.add_argument_group('remove duplicate variants in vcf file')
    vcfdup_args.add_argument('--vcf', required=True, help='Path to vcf file')
    vcfdup_args.add_argument('--outdir', required=True, help='Path to output directory')
    vcfdup_args.add_argument('--bgzip', required=True, help='Path to bgzip script')
    vcfdup_args.add_argument('--tabix', required=True, help='Path to tabix script')

    upload_s3_parser = sub_parser.add_parser('uploadtos3', help='Upload vcf files to s3')
    uploads3_args = upload_s3_parser.add_argument_group('Upload vcf files to s3')
    uploads3_args.add_argument('--vcfdir', required=True, help='Path to vcf files')
    uploads3_args.add_argument('--s3bucket', required=True, help='s3 bucket path')

    # Parse script arguments and invoke functionality
    args = args_parser.parse_args()
    if args.command == 'variantcall':
        if args.ref and args.fastadir and args.outdir and args.minimap2 and args.k8 and args.paftools:
            fasta_to_vcf(ref=args.ref, fasta_dir=args.fastadir, out_dir=args.outdir, minimap2=args.minimap2, k8=args.k8,
                         paftools=args.paftools)
    elif args.command == 'picardmergevcf':
        if args.picard and args.vcfdir and args.outdir:
            merge_vcfs_picard(picard=args.picard, vcf_dir=args.vcfdir, out_dir=args.outdir)
    elif args.command == 'bcfmergevcf':
        if args.bcf and args.bgzip and args.tabix and args.vcfdir and args.outdir:
            merge_vcfs_bcf(args.bcf, args.bgzip, args.tabix, args.vcfdir, args.outdir)
    elif args.command == 'vcftoparquet':
        if args.vcf and args.outdir:
            vcf_to_parquet(args.vcf, args.outdir)
    elif args.command == 'removevcfdup':
        if args.vcf and args.outdir and args.bgzip and args.tabix:
            remove_dup_from_vcf(input_vcf=args.vcf, output_vcf_path=args.outdir, bgzip=args.bgzip, tabix=args.tabix)
    elif args.command == 'uploadtos3':
        if args.vcfdir and args.s3bucket:
            upload_to_s3(vcfdir=args.vcfdir, bucket=args.s3bucket)


if __name__ == "__main__":
    main()
