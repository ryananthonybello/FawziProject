# OPENS THE APPROPRIATE FILES
# IF NEW FILES: change F to the fasta format of your data set
F = open("FUSVertFastaFINAL.txt", 'r')
# IF NEW FILES: chage J to the tabf format of your data set
J = open("FUSVertINFO.txt", 'r')
# IF NEW FILES: change name to new file name
K = open("orthofusvert.txt", 'w')

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

# MAKES A FASTA FILE FORMAT WITH JUST THE ANIMAL NAMES AND  THE SEQUENCES
for gene in range(len(line_content)):
    if (gene == 0):
        continue
    K.write('>')
    K.write(line_content[gene][4])
    K.write('>  ')
    K.write(line_content[gene][7])
    K.write('\n')
    K.write(F_lines[gene-1])
    K.write('\n')


F.close()
J.close()
K.close()
