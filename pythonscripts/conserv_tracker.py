'''
conserv_tracker.py tracks the number of aromatic (F,W,Y) and glycine residues within the orthologue/human pairwise alignment ('info_table_pairwise_tdpmammal.tx')
and reports them in the conserv_tracker_tdpmammal.txt
'''

AroDict = {
'F': ["Phe", "Phenylalanine"],
'W': ["Trp", "Tryptophan"],
'Y': ["Tyr", "Tyrosine"],
}

def extract_str(input_str):
    if input_str is None or input_str == '':
        return 0

    out_string = ''
    for ele in input_str:
        if ele.isdigit():
            continue
        else:
            out_string += ele
    return(out_string)

def extract_nbr(input_str):
    if input_str is None or input_str == '':
        return 0

    out_number = ''
    for ele in input_str:
        if ele.isdigit():
            out_number += ele
    return int(out_number)

E = open("protein_seq_tdp.txt", "r")
F = open("info_table_pairwise_tdpmammal.txt", "r")
W = open("conserv_tracker_tdpmammal.txt", "w")

W.write("NULL")
W.write("\t")
W.write("GLY")
W.write("\t")
W.write("ARO")
W.write('\t')
W.write('TOTAL')
W.write("\n")

'''
The following reports the number of substitutions/deletions at residues that have glycine or aromatics in the human sequence
and places them in the holder array special_holder....the residues that dont have a gly or aro are going to be marked with
'''

residue_holder = []

for line in E:
    line = line.strip()
    line = line.split('\t')
    residue_holder.append(line)

print(residue_holder)

#reading pairwise table and adding up subs for aros, glys, and total
residue_holder_hold = []
letters = []
numbers = []
sub_number_holder = 0
number_total = 0
counter = 0

for line in F:
    if (counter == 0):
        counter =+ 1
        continue
    counter += 1
    number_total = 0
    print(counter)
    residue_holder_hold = str(residue_holder[counter-2])
#taking away the weird residual ['blah'] when taking the item out of residue_holder
    residue_holder_hold = residue_holder_hold.strip("['|']")
    line = line.strip()
    line = line.split('\t')
    letters = extract_str(line)
    print(letters)
    print(line)
    if (letters[0] == 'G'):
        for n in range(len(line)):
            if (n == 0):
                continue
            sub_number_holder = line[n].strip("'")
            if (sub_number_holder.isdigit()):
                print(sub_number_holder)
                sub_number_holder = int(sub_number_holder)
                number_total += sub_number_holder
        W.write(residue_holder_hold)
        W.write('\t')
        W.write(str(number_total))
        W.write('\t')
        W.write('0')
        W.write('\t')
        W.write(str(number_total))
        W.write('\n')
    elif (letters[0]in AroDict):
        for n in range(len(line)):
            if (n == 0):
                continue
            sub_number_holder = line[n].strip("'")
            if (sub_number_holder.isdigit()):
                print(sub_number_holder)
                sub_number_holder = int(sub_number_holder)
                number_total += sub_number_holder
        W.write(residue_holder_hold)
        W.write('\t')
        W.write('0')
        W.write('\t')
        W.write(str(number_total))
        W.write('\t')
        W.write(str(number_total))
        W.write('\n')
    elif (len(line) > 2):
        for n in range(len(line)):
            if (n == 0):
                continue
            sub_number_holder = line[n].strip("'")
            if (sub_number_holder.isdigit()):
                print(sub_number_holder)
                sub_number_holder = int(sub_number_holder)
                number_total += sub_number_holder
        W.write(residue_holder_hold)
        W.write('\t')
        W.write('0')
        W.write('\t')
        W.write('0')
        W.write('\t')
        W.write(str(number_total))
        W.write('\n')
    else:
        W.write(residue_holder_hold)
        W.write('\t')
        W.write('0')
        W.write('\t')
        W.write('0')
        W.write('\t')
        W.write('0')
        W.write('\n')


W.close()
E.close()
F.close()
