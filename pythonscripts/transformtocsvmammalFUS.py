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

human_protein = "MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNSSSGGGGGGGGGGNYGQDQSSMSSGGGSGGGYGNQDQSGGGGSGGYGQQDRGGRGRGGSGGGGGGGGGGYNRSSGGYEPRGRGGGRGGRGGMGGSDRGGFNKFGGPRDQGSRHDSEQDNSDNNTIFVQGLGENVTIESVADYFKQIGIIKTNKKTGQPMINLYTDRETGKLKGEATVSFDDPPSAKAAIDWFDGKEFSGNPIKVSFATRRADFNRGGGNGRGGRGRGGPMGRGGYGGGGSGGGGRGGFPSGGGGGGGQQRAGDWKCPNPTCENMNFSWRNECNQCKAPKPDGPGGGPGGSHMGGNYGDDRRGGRGGYDRGGYRGRGGDRGGFRGGRGGGDRGGFGPGKMDSRGEHRQDRRERPY"

F = open("info_table_pairwise_vertFUS.txt", 'r')
K = open("csvforGNUvertFUS.txt", 'w')

#import the data from info_table_pairwise without first line
F_lines = []
#will hold data for polymorphisms
F_polymorphism = []


for line in F:
    line = line.strip()
    line = line.split('\t')
    line_check = 0
    if (line[1] in aminoConversion):
        line_check += 1
    if (line_check == 0):
        continue
    F_lines.append(line)

print(F_lines)

K.write('\t')
for key in aminoConversion:
    K.write(key)
    K.write ('\t')
K.write('\n')


for index in range(len(F_lines)):
    #line_progress = 0
    aminonumber = F_lines[index][0]
    int_aminonumber = int(aminonumber)
    K.write(aminonumber+human_protein[int_aminonumber-1])
    K.write("\t")
    #line_progress +=1
    line_holder = F_lines[index]
    # initial_line_progress = line_progress
    # if (F_lines[index][line_progress] in aminoConversion):
    for key in aminoConversion:
        match = 0
        for index2 in range(len(line_holder)):
            if (len(line_holder) == 2):
                break
            if (index2 < 2):
                continue
            if (line_holder[index2] == key):
                K.write(line_holder[index2+1])
                K.write('\t')
                match +=1
        if (match == 0):
            K.write('0')
            K.write('\t')
        '''
        line_progress = initial_line_progress
        while (F_lines[index][line_progress] != key):
            line_progress +=2
        if (F_lines[index][line_progress] == key):
                K.write(F_lines[index][line_progress+1])
                #print(F_lines[index][line_progress+1])
                K.write('\t')
                continue
        else:
            K.write('0')
            K.write('\t')
        '''
    K.write('\n')



    '''
    if (isinstance(F_lines[line_progress])
        aminonumber = F_lines[line_progress][index]
    K.write(F_lines[line_progress][index]+human_protein[aminonumber])
    if (F_lines[1][index] in aminoConversion)
    else
        K.write('\t')
    '''
F.close()
