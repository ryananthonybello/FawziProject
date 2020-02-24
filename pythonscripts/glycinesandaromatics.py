aminoConversion = {
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


F = open("info_table_pairwise_vertFUS.txt", "r")
X = open("gly_and_aroFUSvert.txt", "w")

#GLYCINE TRACKER
residue_holder = []
special_holder = []
number = 0

for line in F:
    number += 1
    if (number == 1):
        continue
    line = line.strip()
    line = line.split('\t')
    residue_holder = extract_str(line)
    if (residue_holder[0] == "G" and len(residue_holder) > 1):
        special_holder.append(line)
    else:
        continue

holder1 = []
holder2 = []
x = 0
X.write("GLYCINES")
X.write("\n")
X.write("Residue")
X.write("\t")
X.write("Variants")
X.write("\n")
for n in range(len(special_holder)):
    holder1 = str(special_holder[n][0]) + str(special_holder[n][1])
    X.write(holder1)
    x = 2
    while (len(special_holder[n]) > x):
        X.write("\t")
        holder2 = str(special_holder[n][x]) + str(special_holder[n][x+1])
        X.write(holder2)
        x += 2
    X.write("\n")

#AROMATIC TRACKER
residue_holder1 = []
special_holder1 = []
number = 0

F.close()
F = open("info_table_pairwise_vertFUS.txt", "r")

for line in F:
    number += 1
    if (number == 1):
        continue
    line = line.strip()
    line = line.split('\t')
    residue_holder1 = extract_str(line)
    if (residue_holder1[0] in aminoConversion.keys() and len(residue_holder1) > 1):
        special_holder1.append(line)
    else:
        continue

holder3 = []
holder4 = []
x = 0
X.write("AROMATICS")
X.write("\n")
X.write("Residue")
X.write("\t")
X.write("Variants")
X.write("\n")
for n in range(len(special_holder1)):
    holder3 = str(special_holder1[n][0]) + str(special_holder1[n][1])
    X.write(holder3)
    x = 2
    while (len(special_holder1[n]) > x):
        X.write("\t")
        holder4 = str(special_holder1[n][x]) + str(special_holder1[n][x+1])
        X.write(holder4)
        x += 2
    X.write("\n")


F.close()
X.close()
