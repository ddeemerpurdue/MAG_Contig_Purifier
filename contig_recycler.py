'''
Program to recycle contigs

Example usage:
$ python contig_recycler.py -b <blastfile.csv> -a <assembly.fasta>
  -t 95 -o <newfile.fasta>
OR
$ python contig_recycler.py -b <blastfile.csv> -a <assembly.fasta>
  -t 95 -o <newfile.fasta> --Add -p <path/to/bin/dir>
'''


import subprocess
import argparse
import os
print(os.getcwd())

''' The following is how to add arguments to use the program directly
in the BASH shell '''

''' ARGUMENTS '''
parser = argparse.ArgumentParser(description="Parser")
parser.add_argument("-b", "--Blastfile", help="Blast file (.csv)",
                    required=True)
parser.add_argument("-a", "--Assembly", help="Original assembly file",
                    required=True)
parser.add_argument("-t", "--Threshold", help="Percent identity threshold.\
                    Default = 95%", required=False, default=95)
parser.add_argument("-o", "--Output", help="Output file to write non-binned, \
                    high matching entries to", required=True)
parser.add_argument("-x", "--Add", help="Auto-create new binfiles",
                    required=False, action='store_true')
parser.add_argument("-p", "--Binpath", help="Path to binfiles",
                    required=False)

argument = parser.parse_args()
''' ARGUMENTS '''


def get_nobin_list(blastfile, thresh):
    ''' Function to capture contig names that
    are unbinned, yet align with an identity above
    a user-defined threshold '''
    savecontigs = []
    with open(blastfile) as blast:
        line = blast.readline()

        while line:
            line = line.split(',')
            query = line[0]
            subject = line[2]
            sub_bin = line[3]
            identity = float(line[4])
            # Add to list query contig
            if 'N/A' in sub_bin:
                if identity >= float(thresh):
                    savecontigs.append(subject + '_' + current_bin)
                else:
                    print('Contig does not pass identity threshold')
            else:
                current_bin = sub_bin
            line = blast.readline()
    return savecontigs


def get_nobin_seqs(file, thresh, conts):
    ''' A function that takes the contig names from the above
    function and finds a match in the assembly file. Next, it
    writes - in .fasta format - the sequences to a file.
    Note: since there is discontinuity between contig names in
    blast file and W1.contigs.fasta, there is a step that resolved
    that issue by matching contig numbers '''
    writedict = {}
    count = 0
    # nobin_list is a list of contig names needed to be appended
    nobin_list = get_nobin_list(file, thresh)
    # Open the contigs.fasta file provided
    with open(conts) as cont:
        line = cont.readline()
        while line:
            # Find header file
            if line.startswith('>'):
                # Split it into a list object everywhere a '_' occurs
                # e.g. - 'NODE_10_Length_20' turns into:
                # ['NODE', '10', 'Length', '20']
                line_identifier = line.split('_')[1]  # Grab the node number
                line_identifier = float(line_identifier)
                # Change type above from string to integer
                # Look in nobin_list to see if there's a header that matches
                for nb in nobin_list:
                    # Need to split up the identifier and grab the number
                    # e.g. - 'c_000012' becomes ['c', '000012']
                    nb_match = float(str(nb.split('_')[1]))  # Change type
                    # If there is a match between nobin contig and assembly
                    if line_identifier == nb_match:
                # Keep track of how many contigs match, as this will be used
                # at the end to assert all of our contigs were added
                        count += 1
                        sequence = ''
                        seq = cont.readline().strip()  # strip rms newline char
                # Below is one (of many ways) to capture the sequence.
                # Idea is once you come across a defline ('>'), you read
                # the new line, and while the line does not start with
                # '>', keep reading lines and appending new sequence lines   
                        while not seq.startswith('>'):
                            sequence = sequence + seq
                            seq = cont.readline().strip()
                        else:
                # Once a new entry ('>') starts, write the previous sequence
                # to the dictionary with key == name and value == sequence
                            writedict[nb] = sequence
                # Erase the sequence variable so its good for next round
                            sequence = ''
                    else:
                        pass
            else:
                pass
            line = cont.readline()
    if count != len(nobin_list):
        print('Not all non-binned sequences were found in contig file!')
        # Above, just want to be notified if - for whatever reason - the
        # nonbinned file wasn't in the original assembly. This would be
        # cause for concern and somewhere an error would exist.
    else:
        print('All good.')
    return writedict


# After these two functions, we have 1. grabbed all names of contigs
# without a bin but with a % identity > a threshold. Then took that
# list and compare to the original assembly in order to get the
# sequence and store all .fasta entries into a single key. The next
# function writes those entries to a .fasta file of your choosing.


def write_nobin_seqs(file, thresh, conts, outputfile):
    count = 0
    writedict = get_nobin_seqs(file, thresh, conts)
    print(f'All {len(writedict)} matches include:\n')
    for key in writedict.keys():
        print(key)
    # Open a new file to write data to
    with open(outputfile, 'w') as out:
    # Here, we will want to write key (identifier) on one line and
    # then the following line will be the sequence. This will be
    # repeated for all fasta entries we previously grabbed
        # Step 1: loop through all keys (identifiers)
        for key in writedict.keys():
            # Step 2: write key on first line, followed by newline char
            out.write(str('_'.join(key.split('_')[0:-1])) + '\n')
            # Step 3: write the sequence, followed by newline char
            out.write(writedict[key] + '\n')
            count += 1
    return 'New file created.'


# The next function is only used if the '--add' flag is used, followed
# by the directory containing all of the bin files. This will create a
# new copy of each binfile, containg the original entries and the new
# entries based off blast results.


def write_new_binfiles(file, thresh, conts, outputfile, bin_directory):
    cnt = 0
    # Obtain a dictionary of nonbinners identity/sequence
    nonbinners = get_nobin_seqs(file, thresh, conts)
    # Go through each identifier...
    for key in nonbinners.keys():
        # This will be, for example, bin001 (below)
        ident = key.split('_')[-1]
        # Compare it to files in the provided bin path
        for file in os.listdir(argument.Binpath):
            # If the identity (bin001, etc.) is in the filename
            if ident in file.lower():
                fullfile = os.path.join(argument.Binpath, file)
                print(f'Found match for {key}')
                ident2 = '_'.join(key.split('_')[0:-1])
                # Now check to see sequence isn't already in there
                with open(fullfile) as f:
                    line = f.readline()
                    while line:
                        if ident2 in line:
                            print('Sequence already in file!!!')
                            cnt += 1
                            break
                        else:
                            line = f.readline()
                if cnt < 1:
                    with open(fullfile, 'a') as fasta:
                        fasta.write(str('_'.join(key.split('_')[0:-1])) + '\n')
                        fasta.write(str(nonbinners[key]) + '\n')


if __name__ == '__main__':
    write_nobin_seqs(argument.Blastfile, argument.Threshold,
                     argument.Assembly, argument.Output)
    if argument.Add and (argument.Binpath is not None):
        write_new_binfiles(argument.Blastfile, argument.Threshold,
                           argument.Assembly, argument.Output,
                           argument.Binpath)
