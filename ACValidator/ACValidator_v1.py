
# =================================================================================================================
# Assembly based Circular RNA Validator (ACValidator) workflow
#
# Input: Circular RNA coordinate to be validated
#        Input SAM file to validate the junction in
#
# Output: SAM file with extracted reads and its corresponding BAM, FASTQ
#      Pseudo-reference file with the extracted circRNA sequence extracted from the reference
#      Trinity FASTA and alignment files
#      .txt files containing the overlapping contig sequence, if an overlap exists
#
# Invocation by the wrapper script ACValidator_launcher.sh
#
# Author:
# Shobana Sekar, Liang lab
# Translational Genomics Research Institute, Phoenix, Arizona
#
# Date:
# April 2018
# =================================================================================================================

from __future__ import print_function
import argparse
import pysam
import logging
import os
from subprocess import check_call
import textwrap as tw
import time

log = logging.getLogger(__name__)

# Make sure to change these paths accordingly
# REFPATH should contain the individual chromosome reference FASTA files (eg: 1.fa, 2.fa etc.)
SAMTOOLSPATH = "/packages/samtools/1.4.1/bin/samtools"
BWAPATH = "/packages/bwa/0.7.12/bwa"
BAMTOFASTQPATH = "/packages/BEDTools/2.26.0/bin/bamToFastq"
REFPATH = "/home/ssekar/splitFastas"
TRINITYOUTDIR = "trinity_out"
TRINITYPATH = "/packages/trinityrnaseq/2.4.0/Trinity"


def launch_sh(cmd):
    """Launch a shell command but log the cmd line first"""
    log.info(cmd)
    return check_call(cmd, shell=True)


def convertSamToBam(input_sam):
    """Function to convert SAM to BAM"""
    global output_bam

    output_bam = input_sam.split(".")[0]
    bamfile = os.path.join(output_bam + ".bam")
    sorted_bamfile = os.path.join(output_bam + ".sorted.bam")
    bam_index = os.path.join(output_bam + ".sorted.bam.bai")

    log.info('{}: {}'.format("Output prefix is",output_bam))

    if not os.path.exists(bamfile):
        log.info("Running samtools..")
        cmd = "{} view -b -o {}.bam {}".format(SAMTOOLSPATH, output_bam, input_sam)
        launch_sh(cmd)
    else:
        log.info("Bam exists.. moving on")

    if not os.path.exists(sorted_bamfile):
        log.info("Samtools sorting..")
        cmd = "{} sort -o {} {}.bam".format(SAMTOOLSPATH, sorted_bamfile, output_bam)
        launch_sh(cmd)
    else:
        log.info("Sorted bam exists.. moving on")

    if not os.path.exists(bam_index):
        log.info("Indexing..")
        cmd = "{} index {}.sorted.bam".format(SAMTOOLSPATH, output_bam)
        launch_sh(cmd)
    else:
        log.info("Index exists..")

    time.sleep(2)


def createPseudoRef(circ_coord):
    """Function to create the pseudo-reference containing the circRNA
    sequence"""
    log.info("Creating pseudoref..")

    interval1 = '{}:{}-{}'.format(chrom, endwdw, end)
    interval2 = '{}:{}-{}'.format(chrom, start+1, startwdw+1)

    cmd = "{} faidx {}/{}.fa {} > Scrambled_{}.fa".format(SAMTOOLSPATH, REFPATH, str(chrom), interval1, circ_coordinates)
    launch_sh(cmd)

    cmd = "{} faidx {}/{}.fa {} >> Scrambled_{}.fa".format(SAMTOOLSPATH, REFPATH, str(chrom), interval2, circ_coordinates)
    launch_sh(cmd)

    log.info('Calling parse_fasta')
    parse_fasta('Scrambled_' + circ_coordinates + '.fa', circ_coordinates)

    cmd = "{} index -a bwtsw Pseudoref_{}.fa".format(BWAPATH, circ_coordinates)
    launch_sh(cmd)

    time.sleep(2)


def parse_fasta(infasta, coordinates):
    """Function to re-format the pseudo-reference FASTA file, so that each line in
    FASTA contains the same number of letters (required for IGV checks)"""
    global fasta_string
    fasta_string = ""

    out_f1 = open('Pseudoref_' + circ_coordinates + '.fa', 'w')

    with open(infasta, 'r') as infile:
        line = infile.readline()

        while line != "":
            if line.startswith('>'):
                line = infile.readline()
                continue

            fasta_string = fasta_string + line.strip()
            line = infile.readline()
    out_f1.write('>' + str(coordinates) + '\n')
    out_f1.write(tw.fill(fasta_string, 70))
    out_f1.close()

    time.sleep(2)


def extractRegions(circ_coord):
    """Function to extract reads from a defined window on either side of the
    circRNA junction"""
    global coord_split
    global chrom
    global start
    global end
    global startwdw
    global endwdw

    coord_split = circ_coord.split(":")
    chrom = coord_split[0]
    start = int(coord_split[1].split("-")[0])
    end = int(coord_split[1].split("-")[1])
    
    startwdw = start + int(w)
    endwdw = end - int(w)

    log.info(str(start) + ":" + str(end))
    log.info(str(startwdw) + ":" + str(endwdw))

    if not os.path.exists(circ_coord):
        os.makedirs(circ_coord)

    os.chdir(circ_coord)

    interval1 = '{}:{}-{}'.format(chrom, endwdw, end)
    interval2 = '{}:{}-{}'.format(chrom, start, startwdw)

    cmd = "{} view -h ../{}.sorted.bam \"{}\" > Region_reads_{}.sam".format(SAMTOOLSPATH, output_bam, interval1, circ_coord)
    launch_sh(cmd)

    cmd = "{} view ../{}.sorted.bam \"{}\" >> Region_reads_{}.sam".format(SAMTOOLSPATH, output_bam, interval2, circ_coord)
    launch_sh(cmd)

    cmd = "{} view -b -o Region_reads_{}.bam Region_reads_{}.sam".format(SAMTOOLSPATH, circ_coord, circ_coord)
    launch_sh(cmd)

    cmd = "{} sort -n -o Region_reads_{}.qsorted.bam Region_reads_{}.bam".format(SAMTOOLSPATH, circ_coord, circ_coord)
    launch_sh(cmd)

    cmd = "{} -i Region_reads_{}.qsorted.bam -fq Region_reads_{}.fq".format(BAMTOFASTQPATH, circ_coord, circ_coord)
    launch_sh(cmd)

    cmd = "sed -i '1~4 s/$/\/1/g' Region_reads_{}.fq".format(circ_coord)
    launch_sh(cmd)

    time.sleep(2)


def runTrinity():
    """Function to run trinity on the extracted reads"""
    cmd = "{} --seqType fq --single Region_reads_{}.fq --output {}_{} --CPU 20 --max_memory 100G --full_cleanup".\
           format(TRINITYPATH, circ_coordinates, TRINITYOUTDIR, circ_coordinates)
    launch_sh(cmd)

    time.sleep(2)


def runBwaAlignment():
    """Function to align the trinity assembled FASTA against the
    pseudo-reference"""
    cmd = "{} mem -T 19 Pseudoref_{}.fa trinity_out_{}.Trinity.fasta > Trinity_bwa_{}.sam".\
           format(BWAPATH, circ_coordinates, circ_coordinates, circ_coordinates)
    launch_sh(cmd)

    convertSamToBam("Trinity_bwa_{}.sam".format(circ_coordinates))

    time.sleep(2)


def checkOverlap(infile, fs='\t'):
    """Function to check whether the assembled contigs cross over the created
    pseudo-reference junction. Four different stringency criteria defined, can
    be modified per user requirements"""
    out_f2 = open('Check_overlap_out_highStringency_' + circ_coordinates + '.txt', 'w')
    out_f3 = open('Check_overlap_out_medStringency_' + circ_coordinates + '.txt', 'w')
    out_f4 = open('Check_overlap_out_lowStringency_' + circ_coordinates + '.txt', 'w')
    out_f5 = open('Check_overlap_out_vlowStringency_' + circ_coordinates + '.txt', 'w')

    for read in infile:
        readname = read.query_name
        sequence = read.query_sequence
        log.info("Sequence is: {}".format(sequence))
        log.info("Fasta string is: {}".format(fasta_string))

        junction_seq_hs = fasta_string[(int(w)-30):(int(w)+30)]
        if sequence.find(junction_seq_hs) != -1:
            out_f2.write(fs.join((str(readname), sequence, "Found overlap", "\n")))

        junction_seq_ms = fasta_string[(int(w)-20):(int(w)+20)]
        if sequence.find(junction_seq_ms) != -1:
            out_f3.write(fs.join((str(readname), sequence, "Found overlap", "\n")))

        junction_seq_ls = fasta_string[(int(w)-10):(int(w)+10)]
        if sequence.find(junction_seq_ls) != -1:
            out_f4.write(fs.join((str(readname), sequence, "Found overlap", "\n")))

        junction_seq_vls = fasta_string[(int(w)-5):(int(w)+5)]
        if sequence.find(junction_seq_vls) != -1:
            out_f5.write(fs.join((str(readname), sequence, "Found overlap", "\n")))


def arg_parser():
    """Builds the argparse object"""
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--infile', required=True, help='Input Sam file')

    parser.add_argument('-c', '--coordinate', required=True, help='Input coordinate file')

    parser.add_argument('-w', '--window', required=True, help='Window size')

    parser.add_argument('--log-filename', default=None, help='Filename to save logs')

    parser.add_argument('--log-filemode', default='a', help='File mode for log file')

    parser.add_argument('--log-format',
                        default='[%(asctime)s] %(levelname)s [%(name)s.%(funcNa'
                                'me)s:%(lineno)d] %(message)s',
                        help='Formatting template string')

    parser.add_argument('--log-level', default='INFO', help='Minimum level for log emission')

    return parser


def main(args=None):
    """This is the main function from where every other function is called."""
    parser = arg_parser()
    args = parser.parse_args()

    logging.basicConfig(
        filename=args.log_filename,
        filemode=args.log_filemode,
        format=args.log_format,
        level=getattr(logging, args.log_level)
    )

    global circ_coordinates
    global w

    # We assume the SAM file is in the current directory, and create a custom
    # out directory for the input SAM file
    in_sample = args.infile
    circ_coordinates = args.coordinate
    w = args.window
    in_sam = in_sample + ".sam"
    dir_name = in_sample + "_validation_tests_" + str(w)

    log.info("Infile is: {}".format(in_sam))
    log.info("Coordinates: {}".format(circ_coordinates))

    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    os.chdir(dir_name)

    bamf = os.path.join(in_sample + ".bam")
    sorted_bamf = os.path.join(in_sample + ".sorted.bam")
    bamf_index = os.path.join(in_sample + ".sorted.bam.bai")

    if not os.path.exists(in_sam):
        log.info("mv ../" + in_sam + " $PWD")
        check_call("mv ../" + in_sam + " $PWD", shell=True)

    if not os.path.exists(bamf):
        log.info("mv ../" + bamf + " $PWD")
        check_call("mv ../" + bamf + " $PWD", shell=True)

    if not os.path.exists(sorted_bamf):
        log.info("mv ../" + sorted_bamf + " $PWD")
        check_call("mv ../" + sorted_bamf + " $PWD", shell=True)

    if not os.path.exists(bamf_index):
        log.info("mv ../" + bamf_index + " $PWD")
        check_call("mv ../" + bamf_index + " $PWD", shell=True)

    convertSamToBam(in_sam)
    extractRegions(circ_coordinates)
    createPseudoRef(circ_coordinates)
    runTrinity()
    runBwaAlignment()
    input_bam = pysam.AlignmentFile('Trinity_bwa_' + circ_coordinates + '.sam')
    checkOverlap(input_bam)


if __name__ == '__main__':
    main()
