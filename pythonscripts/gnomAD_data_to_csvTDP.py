import csv

human_protein = "MSEYIRVTEDENDEPIEIPSEDDGTVLLSTVTAQFPGACGLRYRNPVSQCMRGVRLVEGILHAPDAGWGNLVYVVNYPKDNKRKMDETDASSAVKVKRAVQKTSDLIVLGLPWKTTEQDLKEYFSTFGEVLMVQVKKDLKTGHSKGFGFVRFTEYETQVKVMSQRHMIDGRWCDCKLPNSKQSQDEPLRSRKVFVGRCTEDMTEDELREFFSQYGDVMDVFIPKPFRAFAFVTFADDQIAQSLCGEDLIIKGISVHISNAEPKHNSNRQLERSGRFGGNPGGFGNQGGFGNSRGGGAGLGNNQGSNMGGGMNFGAFSINPAMMAAAQAALQSSWGMMGMLASQQNQSGPSGNNQNQGNMQREPNQAFGSGNNSYSGSNSGAAIGWGSASNAGSGSGFNGGFGSSMDSKSSGWGM"

aminoConversion = {
'A': ["Ala","Alanine"],
'C': ["Cys", "Cysteine"],
'D': ["Asp","Aspartic Acid"],
'E': ["Glu", "Glutamic Acid"],
'F': ["Phe", "Phenylalanine"],
'G': ["Gly", "Glycine"],
'H': ["His", "Histidine"],
'I': ["Ile", "Isoleucine"],
'K': ["Lys", "Lycine"],
'L': ["Leu", "Leucine"],
'M': ["Met", "Methionine"],
'N': ["Asn","Asparagine"],
'P': ["Pro", "Proline"],
'Q': ["Gln", "Glutamine"],
'R': ["Arg","Arginine"],
'S': ["Ser", "Serine"],
'T': ["Thr", "Threonine"],
'V': ["Val", "Valine"],
'W': ["Trp", "Tryptophan"],
'Y': ["Tyr", "Tyrosine"],
'-': ["Gap", "Gap"]
}

def extract_nbr(input_str):
    if input_str is None or input_str == '':
        return 0

    out_number = ''
    for ele in input_str:
        if ele.isdigit():
            out_number += ele
    return int(out_number)

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


#O = open("ExACFUS.csv", "r")
F = open("protein_seq_tdp.txt","r")
K = open("EMPTYCSV.txt", "w")
W = open("trash.txt", "w")

K.write("NULL")
K.write('\t')

#filling out first row of CSV document
first_line = list(aminoConversion.keys())

for n in range(len(first_line)-1):
    K.write(first_line[n])
    K.write('\t')

K.write('DEL')
K.write('\n')
'''
#hold ExAC csv data in csv_lines
csv_lines = []
counter = 0
'''
protein_consequence_hold = []
allele_frequency = []
occurence_number = []
singleton_check = []
with open("TDPgnomADv3CSV.csv", newline='') as csvfile:
    csv_lines = csv.reader(csvfile)
    count = 0
    for row in csv_lines:
        if(count == 0):
            count +=1
            continue
        print(row[9])
        protein_consequence_hold.append(row[9])
        '''allele_frequency = float(row[16])'''
        allele_frequency = row[15]
        occurence_number.append(allele_frequency)
        '''
        singleton_check = int(row[11])
        if (singleton_check == 1):
            occurence_number.append('1')
            continue
        if (float(allele_frequency) < 1/10000):
            occurence_number.append('2')
            continue
        if (float(allele_frequency) < 1/1000):
            occurence_number.append('3')
            continue
        if (allele_frequency < .01):
            occurence_number.append('4')
            continue
        if (allele_frequency < .05):
            occurence_number.append('5')
            continue
        if (allele_frequency < .5):
            occurence_number.append('6')
            continue
        if (allele_frequency >= .5):
            occurence_number.append('7')
            continue
        '''

for n in range(len(protein_consequence_hold)):
    W.write(protein_consequence_hold[n])
    W.write('\t')
    W.write(occurence_number[n])
    W.write('\n')


'''
aa_residue = []
aa_residue_num = []
exac_poly = []
exac_poly_hold = []
exac_poly_num = []
exac_poly_num_track = 0
aminoConversion_check = aminoConversion.items()
print(aminoConversion_check)
for n in range(len(human_protein)):
    exac_poly_num = extract_nbr(protein_consequence_hold[exac_poly_num_track])
    exac_poly_hold = extract_str(protein_consequence_hold[exac_poly_num_track])
    exac_poly = exac_poly_hold[5:8]
    print(exac_poly)
    print(exac_poly_num)
    aa_residue_num = n
    aa_residue = str(aa_residue_num) + human_protein[n]
    K.write(aa_residue)
    K.write('\t')
    if(aa_residue_num == exac_poly_num):
        print(aa_residue_num)
        for n in range(len(first_line)-1):
            if (exac_poly in aminoConversion.get(first_line[n])):
                print('in')
                K.write(occurence_number[n])
                K.write('\t')
            else:
                K.write('0')
                K.write('\t')
        exac_poly_num_track += 1
    else:
        for n in range(len(first_line)-1):
            K.write('0')
            K.write('\t')
    K.write('\n')
'''

aa_residue = []
aa_residue_num = []
aminoConversion_check = aminoConversion.items()
print(aminoConversion_check)
for n in range(len(human_protein)):
    aa_residue_num = n+1
    aa_residue = str(aa_residue_num) + human_protein[n]
    K.write(aa_residue)
    K.write('\t')
    for n in range(len(first_line)):
        K.write('0')
        K.write('\t')
    K.write('\n')

K.close()
K = open("EMPTYCSV.txt","r")
X = open("TDPgnomADv3.txt","w")
aa_residue_num = 0
exac_poly = []
exac_poly_hold = []
exac_poly_num = []
exac_poly_num_track = 0
aminoConversion_check = aminoConversion.items()
protein_length = range(len(human_protein))
aa_residue_num = 0
line_holder = []
print(len(occurence_number))
print(len(protein_consequence_hold))
print(protein_consequence_hold)
for line in K:
    line = line.strip()
    line_holder = line.split('\t')
    if (aa_residue_num == 0):
        for n in range(len(line_holder)):
            X.write(line_holder[n])
            if(n < len(line_holder)):
                X.write('\t')
        X.write('\n')
        aa_residue_num +=1
        continue
    if (exac_poly_num_track <= len(protein_consequence_hold)-1):
        exac_poly_num = extract_nbr(protein_consequence_hold[exac_poly_num_track])
        exac_poly_hold = extract_str(protein_consequence_hold[exac_poly_num_track])
        exac_poly = exac_poly_hold[5:8]
    while(aa_residue_num == exac_poly_num and exac_poly_num_track <= len(protein_consequence_hold)-1):
        for n in range(len(first_line)-1):
            if(exac_poly_num_track <= len(occurence_number)-1):
                if (exac_poly_num_track <= len(protein_consequence_hold)-1):
                    exac_poly_num = extract_nbr(protein_consequence_hold[exac_poly_num_track])
                    exac_poly_hold = extract_str(protein_consequence_hold[exac_poly_num_track])
                    exac_poly = exac_poly_hold[5:8]
                    print(protein_consequence_hold[exac_poly_num_track])
                if (exac_poly in aminoConversion.get(first_line[n]) and aa_residue_num == exac_poly_num):
                    line_holder[n+1] = occurence_number[exac_poly_num_track]
                    exac_poly_num_track += 1
                print(exac_poly_num_track)
    for n in range(len(line_holder)):
        X.write(line_holder[n])
        if(n < len(line_holder)):
            X.write('\t')
    X.write('\n')
    aa_residue_num += 1



for line in K:
    line = line.strip()
    line = line.split('\t')
    print(line)


'''
#protein_consequence will hold the protein consequence column in a list
protein_consequence = []

#build protein_consequence list
for n in range(len(csv_lines)):
    protein_consequence.append(csv_lines[n][5])
'''


K.close()
F.close()
W.close()
X.close()
'''
consequence = []
#holds the number taken from protein consequnce column in csv
poly_loc = []
# holds ints of poly_loc to check against residue list
poly_loc_num = []
#going through each protein_consquence entry
for index in range(len(protein_consequence)):
    consequence = protein_consequence[index]
    poly_loc = "".join([n for n in consequence if n.isdigit()])
    poly_loc = int(poly_loc)
    poly_loc_num[in]

    #for n in range(len(protein_consequence[index]))




for index in range(len(human_protein)):
    counter = str(index+1)
    F.write(counter+human_protein[index])
    F.write('\n')

F.close()
'''
