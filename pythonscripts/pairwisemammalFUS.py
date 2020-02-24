from Bio import Align
from Bio import AlignIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

matrix= matlist.blosum62

#########################################
# TODO: UPDATE FOR DATASET EXCLUDED FROM THE MSA BECUASE OF CUTOFF,
# DoNotInclude = ["Ovis aries","Fukomys damarensis", "Camelus ferus","Homo sapiens","Condylura cristata","Tupaia chinensis"]

# TODO: UPDATE WITH SEQUENCE
# source: isoform1 from uniprot...claimed to be the canonical sequence
human_protein = "MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNSSSGGGGGGGGGGNYGQDQSSMSSGGGSGGGYGNQDQSGGGGSGGYGQQDRGGRGRGGSGGGGGGGGGGYNRSSGGYEPRGRGGGRGGRGGMGGSDRGGFNKFGGPRDQGSRHDSEQDNSDNNTIFVQGLGENVTIESVADYFKQIGIIKTNKKTGQPMINLYTDRETGKLKGEATVSFDDPPSAKAAIDWFDGKEFSGNPIKVSFATRRADFNRGGGNGRGGRGRGGPMGRGGYGGGGSGGGGRGGFPSGGGGGGGQQRAGDWKCPNPTCENMNFSWRNECNQCKAPKPDGPGGGPGGSHMGGNYGDDRRGGRGGYDRGGYRGRGGDRGGFRGGRGGGDRGGFGPGKMDSRGEHRQDRRERPY"

aminoConversion = {
'-': ["Gap", "Gap"],
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
}

# OPENS THE APPROPRIATE FILES
# TODO: change F to the fasta format of the new data set
# currently using old fasta without integrated species names
F = open("FUSVertFastaFINAL.txt", 'r')
# TODO: chage J to the tabf format of the new data set
J = open("FUSVertINFO.txt", 'r')
'''
# Amber's old intialization....infers that the fasta/tab delimited from orthodb is usable as downloaded
F_lines = []

#  PARSES THE FASTA FILE WITH THE ANIMAL NAMES
for line in F:
    line = line.strip()
    if (line[0]=='>'):
        continue
    # makes a list of the animal names and the associated sequences
    F_lines.append(line)

#  WILL HOLD THE SEQUENCES FOR  THAT ANIMAL
line_content = []

# GOES THROUGH THE OTHER FILE
for line in J:
    line = line.strip()
    line_content.append(line.split('\t'))
'''

F_lines = []
#holds the int_prot_id from the metadata lines in the fasta file to compare to those in the tab delimited file later
F_codeholder = []
#holds all the int_prot_id from F_codeholder
F_code = []

#  PARSES THE FASTA FILE WITH THE ANIMAL NAMES
for line in F:
    line = line.strip()
    if (line[0]=='>'):
    #saving int_prot_id from file to compare to those in tab delimited to eliminate names from sequences that were not used
        F_codeholder = line[line.find(">")+1:line.find(" {")]
        F_code.append(F_codeholder)
        continue
    # makes a list of the animal names and the associated sequences
    F_lines.append(line)
    # print F_code

#  WILL HOLD THE SEQUENCES FOR  THAT ANIMAL
line_content = []
# use this to progressively go through F_code as you check ids
F_index = 0

# GOES THROUGH THE OTHER FILE
for line in J:
    line = line.strip()
    #checks the int_prot_id from F_code against the one in the current line; if there is not a match, then we move onto the next line to check
    line = line.split('\t')
    #print line
    #print F_code[F_index]
    if (F_code[F_index]!= line[5]):
        continue
    line_content.append(line)
    #print line_content
    #upon successfull check and append, the counter for F_code goes up, moving onto the next int_prot_id
    F_index += 1

mammalSEQ = {}
for gene in range(len(line_content)):
    if (gene == 0):
        continue
    if(not line_content[gene][4] in mammalSEQ.keys()):
        mammalSEQ[line_content[gene][4]] = F_lines[gene-1]
    else:
        name = line_content[gene][4]
        name += "1"
        mammalSEQ[name] = F_lines[gene-1]

F.close()
J.close()

################################################################################
################################################################################
################################################################################
# Creates a new FASTA format file
P = open("upadted_sequences_vertFUS.txt", 'w')


# will hold the list of mutations at each location
mutationLIST = []
for amino in human_protein:
    mutationLIST.append({})

counter = 0;
# for each sequence in the mammals sequences
for key in mammalSEQ.keys():
    has_unk = False
    # double check if the sequence has errors
    for amino in mammalSEQ[key]:
        # if so skip them
        if (amino == 'X'):
            has_unk = True
    '''
    # only update the counter for ones that stay in
    if (has_unk or key in DoNotInclude):
        continue
    '''
    P.write('>')
    P.write(key)
    P.write('\n')
    P.write(mammalSEQ[key])
    P.write('\n')
    counter += 1
    #TODO: update with desired alignment metrics
    # get the pairwise alignment
    alignments2 = pairwise2.align.globalds(human_protein,mammalSEQ[key],matrix,-10,-0.5)
    human_holder = alignments2[0][0]
    animal_holder = alignments2[0][1]
    # testing
#    if(counter ==1):
#        print(human_holder)
#        print(animal_holder)
    # new counter for the human sequence
    counter2 = 0
    # for each amino acid in the alignment
    for aminoacid in range(len(human_holder)):
        # if there is a gap in the human side
        if(human_holder[aminoacid] == '-'):
            continue
        # update the number of amino acids without the gaps
        counter2 +=1
        # if mutation found
        if(human_holder[aminoacid] != animal_holder[aminoacid]):
            if(animal_holder[aminoacid] in mutationLIST[counter2-1].keys()):
                mutationLIST[counter2-1][animal_holder[aminoacid]]+=1
            else:
                mutationLIST[counter2-1][animal_holder[aminoacid]]=1

print(counter)
P.close()
'''
Q = open("mammalalignment.txt", 'w')
Q.write(alignments2)
Q.close()
'''
# FOR TESTING print(mutationLIST)

##################### PARSES THE PTMS FILE #####################################

################################################################################
'''
# TODO: update to your file format
# READ THE FILE OF PTMS
# IFNEW Change file name
Q = open("PTMs.csv", 'r')
holder = []
PTMS = []
# FOR EACH LINE IN THE FILE
for line in Q:
    #STRIP THE LINE OF WHITE SPACE
    line = line.strip()
    #SPLIT THE LINE ON A COMMA
    holder = line.split(",")
    # ONLY HOLD ON TO THE RELEVENT STUFF
    holder=holder[:3]
    # IF ITS EMPTY DON'T INCLUDE IT
    if(holder[0]==""):
        continue
    # ADD THAT TO OUR LIST OF PTM LISTS
    PTMS.append(holder)
# FOR SOME REASON THIS SETS TO AN UNK CHAR AND THIS HARDCODE RESETS IT
PTMS[0][0]=7
Q.close()
'''
##################### PARSES THE MUTASTIONS FILE ###############################

################################################################################
'''
# reads a file of mutations would need to be updated to your file format
M = open("mutationslisted.csv",'r')
holder = []
MUTATIONS = []
#FOR EACH LINE IN THE FILE
for line in M:
    #STRIP THE LINE OF WHITE SPACE
    line = line.strip()
    #SPLIT THE LINE ON A COMMA
    holder = line.split(",")
    # IF ITS EMPTY DON'T INCLUDE IT
    if(holder[0]==""):
        continue
    MUTATIONS.append(holder)
MUTATIONS = MUTATIONS[1:]
M.close()
'''
################################################################################
##################### CREATES A TABLE OF INFO ##################################
################################################################################

updated = 0
amino = ''
O = open("info_table_pairwise_vertFUS.txt",'w')
O.write("loc")
O.write('\t')
O.write("AA")
O.write('\t')
#O.write("PTM?")
#O.write('\t')
O.write("total")
O.write('\t')
#O.write("MUT?")
#O.write('\t')
O.write('\t')
O.write("Vars")
O.write('\n')
for aminoacid in range(len(mutationLIST)):
    num_variants=0
    amino = human_protein[aminoacid]
    # write the original number
    O.write(str(aminoacid+1))
    O.write('\t')
    # write the amino acid at that location

    O.write(human_protein[aminoacid])
    O.write('\t')

    '''
    # write if ptms were FOUND
    done = False
    # For each PTM in the list of PTMS
    for mod in range(len(PTMS)):
        # check if it applies to that specific amino acid
        if((int(PTMS[mod][0])-1)==aminoacid and PTMS[mod][1]==amino and not done):
            # if so update the bool holder and write to that file
            O.write('T')
            done = True
    if (not done):
        O.write('F')
    O.write('\t')
    '''
    '''
    for key in mutationLIST[aminoacid]:
        num_variants += int(mutationLIST[aminoacid][key])
    O.write(str(num_variants))
    O.write('\t')
    '''
    '''
    found = False
    for value in MUTATIONS:
        if(int(value[0])-1==aminoacid and not found):
            O.write(value[1]+"-"+value[2])
            found = True
    O.write('\t')
    '''
    for key in mutationLIST[aminoacid]:
        if (key in aminoConversion):
            O.write(key)
            O.write('\t')
    #    elif(key =="-"):
    #        O.write("del")
    #        O.write('\t')
        else:
            continue
        O.write(str(mutationLIST[aminoacid][key]))
        O.write('\t')
    O.write('\n')

O.close()
