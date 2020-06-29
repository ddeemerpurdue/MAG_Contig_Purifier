"""
Program to subset the dataframe based on specific
attributes that we are looking for
"""
import pandas as pd

myfile = 'Skeleton/Nobin_Bin.txt'
binfile = 'Skeleton/Bin_Bin.txt'
testbin = 'Skeleton/BB1000.txt'

def get_cois(data):
    cnames = ['Query', 'Reference', 'Qnode', 'Rnode', 'ANI',
          'Orthos', 'Total', 'Qbin', 'Rbin']
    df = pd.read_csv(data, sep='\t',
                 names=cnames)
    df = df[df.ANI >= 90.0]
    nobin_nobin = df[(df.Qbin == 'NoBin') & (df.Rbin == 'NoBin')]
    nobin_bin = df[(df.Qbin == 'NoBin') & (df.Rbin != 'NoBin')]
    bin_bin = df[(df.Qbin != 'NoBin') & (df.Rbin != 'NoBin')]
    return nobin_nobin, nobin_bin, bin_bin


def write_cois(data, out1, out2, out3):
    nb_nb, nb_b, b_b = get_cois(data)
    nb_nb.to_csv(out1, sep = '\t', index=False)
    nb_b.to_csv(out2, sep = '\t', index=False)
    b_b.to_csv(out3, sep = '\t', index=False)
    print('Done!')


def create_nobin_dictionary(file):
    df = pd.read_csv(file, sep='\t')
    large_dictionary = {}
    for quer in df.Query.unique():
        reference_dict = {}
        for refer in df.Reference.unique():
            if quer == refer:
                pass
            else:
                a = df[(df.Query == quer) & (df.Reference == refer)]
                contig_dict = {}
                for i in range(a.shape[0]):
                    qcontig = a.iloc[i, 2]
                    rbin = a.iloc[i, 8]
                    if qcontig in contig_dict:
                        contig_dict[qcontig].append(rbin)
                    else:
                        contig_dict[qcontig] = [rbin]
                    # large_dictionary[quer] = {refer}
                reference_dict[refer] = contig_dict
        large_dictionary[quer] = reference_dict
    return large_dictionary


# a = create_nobin_dictionary(myfile)
# print(a)


def create_bin_dictionary(file):
    df = pd.read_csv(file, sep='\t')
    large_dictionary = {}
    for quer in df.Query.unique():
        print(quer)
        reference_dict = {}
        for refer in df.Reference.unique():
            if quer == refer:
                pass
            else:
                contig_dict = {}
                for bins in df.Qbin.unique():
                    a = df[(df.Query == quer) &
                           (df.Reference == refer) &
                           (df.Qbin == bins)]
                    for i in range(a.shape[0]):
                        qbin = a.iloc[i, 7]
                        rbin = a.iloc[i, 8]
                        # Make a dictionary where key == qbin and values
                        # are associated rbin
                        if qbin in contig_dict:
                            contig_dict[qbin].append(rbin)
                        else:
                            contig_dict[qbin] = [rbin]
                # Count occurences of each Rbin per Qbin
        # At the end, contig_dict would contain all unique Qbins
        # for Quer/Refer
                # Go through each Qbin key
                for key in contig_dict.keys():
                    cnt_list = []
                    # Loop through the set of values (unique)
                    for l in set(contig_dict[key]):
                        # If it's first value...
                        if cnt_list == []:
                            cnt_list.append(l)
                            cnt_list.append(contig_dict[key].count(l))
                            # Ex. ['Bin013', 13]
                        else:
                            # If multiple Rbins
                            if contig_dict[key].count(l) > cnt_list[1]:
                                cnt_list.append(l)
                                cnt_list.append(contig_dict[key].count(l))
                            # Ex. ['Bin013', 13, 'Bin011', 4]
                            else:
                                pass
                    # Now have key go from each occurence to counts
                    contig_dict[key] = cnt_list
                    # Ex. {'qBin002': ['rBin13', 13, 'rBin011', 4]}
                    # Proceed for all qBins
                reference_dict[refer] = contig_dict
        large_dictionary[quer] = reference_dict
    return large_dictionary


# b = create_bin_dictionary(testbin)
# print(b)


def map_unbinned_contigs(nobin_dict, bin_dict):
    lines = []
    nonbins = create_nobin_dictionary(nobin_dict)
    binners = create_bin_dictionary(bin_dict)
    # 1 - Loop through all queries
    for query in nonbins.keys():
        # 2 - Loop through all references
        for reference in nonbins[query].keys():
            # 3 Go through all unbinned contigs
            for contig in nonbins[query][reference].keys():
            # 4 Capture their nearest bin match
                bins_to_match = nonbins[query][reference][contig]
            # Ex. {R1: {R2: {Node_10_blah: ['bin015', 'bin001']}}}
            # 5 Since sometimes multiple matches, loop through each match
                for bin_to_match in bins_to_match:
            # Ex. 'bin015' - Saying contig from R1 matched bin 15 from R2
            # 6 Search 'binners' to find match
            # Where do other contigs from R1 match R2 bins?
                    # Go through all query bins from bin-bin dict
                    for quer_bin in binners[query][reference].keys():
            # 7 If there's a match between Rbin and Qbin/Rbin
                        if (bin_to_match in binners[query][reference][quer_bin] and
                                binners[query][reference][quer_bin][1] > 10):
            # This would match binners {R1: {R2: {Qbin:['Rbin']}}} -- 'Rbin'
                            bin_match = quer_bin
                        else:
                            bin_match = "No_Recycle"
                    lines.append(f"{query}\t{reference}\t{contig}\t{bin_to_match}\t{bin_match}\n")
    return lines


def write_recycled_bins(nobin, bindic, output):
    a = map_unbinned_contigs(nobin, bindic)
    with open(output, 'w') as o:
        for line in a:
            o.write(line)
    return None


# print(f"Query: {quer}\tReference: {refer}")

# write_recycled_bins(myfile, binfile, 'ex.txt')

with open('ex.txt') as f:
    lengths = []
    line = f.readline()  # Skip header
    line = f.readline()
    while line:
        if line.split('\t')[4].strip() == 'No_Recycle':
            print(line.split('\t')[4])
            line = f.readline()
        else:
            length = line.split('\t')[2].split('_')[3]
            length = int(length)
            lengths.append(length)
            line = f.readline()
mean = sum(lengths)/len(lengths)
print(mean)
print(len(lengths))

import matplotlib.pyplot as plt

plt.hist(x=lengths, bins = 100)
plt.show()