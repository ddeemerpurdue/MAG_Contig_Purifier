'''
Program designed to grab all recycled contigs and see how much nucleotide
content they overlap with in original binners

Arguments needed:
Master ANI file
File with recycled contigs, bins, etc.
Example usage:
$ python recycled_contigs_overlap.py -a <ani.txt> -t 90 -r 0.75
-y <recycled.txt> -o <output_file.txt>
'''

import os
import argparse


""" Arguments """
parser = argparse.ArgumentParser(description="Parser")
parser.add_argument("-a", "--ANI", help="ANI master file", required=True)
parser.add_argument("-t", "--Threshold", help="ANI threshold to consider \
matches", default=90)
parser.add_argument("-r", "--Ratio", help="Ratio of orthos to matches",
                    default=0.75)
parser.add_argument("-y", "--Recycled", help="Suspect recycled contigs file",
                    required=True)
parser.add_argument("-o", "--Output", help="Output file to write to",
                    required=True)
argument = parser.parse_args()


def read_ANI_line(rline):
    line = rline.split('\t')
    quer = line[0]
    ref = line[1]
    c_rec = line[2]
    c_ori = line[3]
    ani = float(line[4])
    rat = (float(line[5]) / float(line[6]))
    qbin = line[8].strip()
    return line, quer, ref, c_rec, c_ori, ani, rat, qbin


def read_Rec_line(rline):
    line = rline.split('\t')
    quer = line[0]
    ref = line[1]
    c_rec = line[2]
    ref_bin = line[3]
    quer_bin = line[4].strip()
    return line, quer, ref, c_rec, ref_bin, quer_bin


def return_overlapping_contigs(anifile, threshold, ratio):
    """
    Function takes in the ANI file and returns a
    dictionary of potential overlapping contigs
    """
    ani_d = {'R1':{}, 'R2':{}, 'R3':{}, 'R4':{},
             'W1':{}, 'W2':{}, 'W3':{}, 'W4':{}}
    with open(anifile) as af:
        line = af.readline()
        while line:
            line, quer, ref, c_rec, c_ori, ani, rat, qbin = read_ANI_line(line)
            if (quer == ref) and (c_rec != c_ori):
                if (ani > float(threshold)) and (rat >= ratio):
                    ani_d[quer][c_rec] = [c_ori, ani, rat, qbin]
            else:
                pass
            line = af.readline()
    return ani_d


def read_recycled_contigs(rec):
    rec_d = {'R1':{}, 'R2':{}, 'R3':{}, 'R4':{},
             'W1':{}, 'W2':{}, 'W3':{}, 'W4':{}}
    with open(rec) as r:
        line = r.readline()  # Skip header
        line = r.readline()
        while line:
            line, quer, ref, c_rec, ref_bin, quer_bin = read_Rec_line(line)
            try:
                # Remain open to possibility a contig could be directed
                # towards to multiple bins
                if quer_bin in rec_d[quer][c_rec]:
                    pass
                else:
                    rec_d[quer][c_rec].append(quer_bin)
            except KeyError:
                rec_d[quer][c_rec] = [quer_bin]
            line = r.readline()
    return rec_d


def compare_rec_overlapping(anifile, threshold, ratio, rec):
    master = []
    overlaps = return_overlapping_contigs(anifile, threshold, ratio)
    recycled = read_recycled_contigs(rec)
    for quer in overlaps.keys():
        for r_con in overlaps[quer].keys():
            if (r_con in recycled[quer] and
                overlaps[quer][r_con][3] == recycled[quer][r_con][0]):
                c_ori = overlaps[quer][r_con][0]
                ani = overlaps[quer][r_con][1]
                ratio = overlaps[quer][r_con][2]
                cbin = overlaps[quer][r_con][3]
                rbin = recycled[quer][r_con][0]
                line = f"{quer}\t{r_con}\t{c_ori}\t{ani}\t{ratio}\t{cbin}\t{rbin}\n"
                master.append(line)
    return master


def write_suspect_contigs(anifile, threshold, ratio, rec, output):
    outlist = compare_rec_overlapping(anifile, threshold, ratio, rec)
    with open(output, 'w') as o:
        header = f"Sample\tRecycled_Contig\tNearest_Match\tANI\tRatio\tNearest_Bin\tRecycled_Bin\n"
        o.write(header)
        for line in outlist:
            o.write(line)
    return 'Complete'


if __name__ == "__main__":
    write_suspect_contigs(argument.ANI, float(argument.Threshold),
                          float(argument.Ratio), argument.Recycled,
                          argument.Output)
