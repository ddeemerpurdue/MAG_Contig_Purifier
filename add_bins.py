'''
-BLAH
'''

def get_bindict(binfile):
    # Start with query bins
    bindict = {}
    with open(binfile) as qf:
        line = qf.readline()
        while line:
            node = line.split('\t')[0].split('_')[1]
            node = str(int(node))
            bin_name = line.split('\t')[1].strip()
            bindict[node] = bin_name
            line = qf.readline()
    return bindict


def get_all_bindicts(binlist):
    master = {}
    for bins in binlist:
        name = bins.split('/')[1].split('-')[0]
        print(name)
        master[name] = get_bindict(bins)
    return master


def append_bins_to_ani(binlist, anifile):
    writelist = []
    bins = get_all_bindicts(binlist)
    # Parse through ANI file
    with open(anifile) as af:
        constant_line = af.readline()
        while constant_line:
            line = constant_line.split('\t')
            quer_num = line[0].split('/')[2].split('_')[1]
            ref_num = line[1].split('/')[2].split('_')[1]
            query = line[0].split('/')[1]
            reference = line[1].split('/')[1]
            quer_node = line[0].split('/')[2].rsplit('.', 1)[0]
            reference_node = line[1].split('/')[2].rsplit('.', 1)[0]
            ani = line[2]
            orths = line[3]
            total = line[4].strip()
            # print(f"Query={query}, Ref={reference}, Nums={quer_num}\t{ref_num},\nNodes:\t{quer_node}\t{reference_node}\nANI:\t{ani}")
            # See if there's a bin for the nodes
            if quer_num in bins[query]:
                quer_bin = bins[query][quer_num]
            else:
                quer_bin = 'NoBin'
            # Repeat for reference
            if ref_num in bins[reference]:
                ref_bin = bins[reference][ref_num]
            else:
                ref_bin = 'NoBin'
            # Write all information to a file
            new_info = '\t'.join([query, reference, quer_node, reference_node,
                                  ani, orths, total, quer_bin, ref_bin])
            writelist.append(new_info)
            constant_line = af.readline()
    return writelist




def write_new_ani(binlist, anifile, output):
    writelist = append_bins_to_ani(binlist, anifile)
    with open(output, 'w') as o:
        for line in writelist:
            o.write(line + '\n')
    print('Done')


mylist = ['Bin-Files/W1-bins.txt', 'Bin-Files/W2-bins.txt',
          'Bin-Files/W3-bins.txt', 'Bin-Files/W4-bins.txt',
          'Bin-Files/R1-bins.txt', 'Bin-Files/R2-bins.txt',
          'Bin-Files/R3-bins.txt', 'Bin-Files/R4-bins.txt']
# file = 'W1/W1.R1.ANI.5000f.txt'
# qb = 'Bin-Files/W1-bins.txt'
# rb = 'Bin-Files/R1-bins.txt'


write_new_ani(mylist, 'Skeleton/All.master.5000.txt',
              'Skeleton/All.5000.txt')
