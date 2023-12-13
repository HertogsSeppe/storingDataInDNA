from config import codons


def bits_to_bases(bits, nr_codons=3):
    if nr_codons == 3:
        wl = 16
    if nr_codons == 4:
        wl = 22

    id_secuence = ""
    id = 0

    for i in range(3):
        index = (id // (47**i)) % 47
        id_secuence = codons[index] + id_secuence

    padding = (wl - len(bits) % wl) % wl

    bits = bits + padding * "0"

    id_secuence = id_secuence + codons[padding]

    strand = id_secuence

    for i in range(len(bits) // wl):
        word = int(bits[wl * i : wl * i + wl], 2)

        res = ""

        for i in range(nr_codons):
            index = (word // (47**i)) % 47

            res = codons[index] + res
        strand = strand + res

    return strand

#""" 
#Stack overflow, Brunão: https://stackoverflow.com/questions/75611160/
#converting-any-file-to-binary-1-and-0-and-then-back-to-original-file-without-cor
#"""

#functions below here are for decoding process. need to be reorganized and maybe put somewhere else than util?

def DNA_matching(DNA_strings):
    #assumes DNA_strings is a list containing each sequenced nucleotide string. Can be changed if necessary.
    
    newly_sorted = []
    strings_copy = DNA_strings.copy()

    for i in range(0,len(DNA_strings,1)):

        DNA_seq = strings_copy(0)
        complement_seq = complement_str(DNA_seq)

        if complement_seq in strings_copy:
            
            newly_sorted = newly_sorted + [DNA_seq, complement_seq]

            strings_copy.remove(DNA_seq)
            strings_copy.remove(complement_seq)

    return '_'


def complement_str(line):
    #creates the complement for a DNA string.

    complement = ""
    for i in range(0,len(line),1):
        if line[i] == "A":
            complement = complement + "T"

        elif line[i] == "T":
            complement = complement + "A"

        elif line[i] == "G":
            complement = complement + "C"

        elif line[i] == "C":
            complement = complement + "G"

        else:
            #since illumina sequencing machines can also give an ambiguous result this filters those results out.
                #e.g. returns "R" if it could be A or G, or returns "Y" if it assumes C or T
            #for our purposes we assume no ambiguity (i think), so this could be left out.
            complement = complement + "X"

    return complement


def restructure_B_to_A(base47_list, index_len = 3, redun_B_len = 3):
    #assumes base_list is a list of base 47 column lists that have already been correctly ordered by index nr. if not, has to be changed.
    #assumes all base string have the same length. If not, has to be changed.
    #   - redun_A_len is the horizontal (row) redundancy.
    #   - redun_B_len is the vertical (column) redundancy.

    base47_Btrim_list = [0]*len(base47_list)

    for i in range(0,len(base47_list),1):

        base47_str = base47_list[i]

        base47_Btrim_list[i] = base47_str[index_len:(len(base47_str)-redun_B_len)]

    col_len = len(base47_Btrim_list[0])
    row_len = len(base47_Btrim_list)

    b47_rows_list = [0]*col_len

    for j in range(0, col_len, 1):

        row = [0]*row_len

        for k in range(0, row_len, 1):

            col = base47_Btrim_list[k]
            row[k] = col[j]
        
        b47_rows_list[j] = row

    return b47_rows_list###


def restructure_A_to_final(base47_list, redun_A_len = 3):
    #assumes base47_list is a list of base 47 row lists.

    col_len = len(base47_list)
    row_len = len(base47_list[0]-redun_A_len)

    restructured_list = []

    for j in range(0, col_len, 1):
        base47_trimmed_row = base47_list[k][redun_A_len:]

        restructured_list.extend(base47_trimmed_row)

    return restructured_list


def bases_to_47(bases, codon_len=3):
    #assumes "bases" is a string of bases, so only processes 1 DNA string at a time. can be changed if necessary.

    bases = bases[12:] #staat nu voor een index van 3 nrs en met 1 extra bit voor zero-indication. Moet nog variabel / aanpasbaar.
    index = 0
    base_47_numbers = []

    for i in range(0, int(len(bases)/codon_len), 3):

        codon = bases[index: index + codon_len]

        if codon in codons:
            base_47_numbers.append(codons.index(codon))
        else:
            base_47_numbers.append('XXX') #moet nog kijken wat precies te zetten, of om gewoon weg te laten en te tellen als deletie
    
    return base_47_numbers


def b47_to_binary(b47s, codon_len=3):
    #assumes b47s is a list of base 47 numbers. if not, can be changed.

    if codon_len == 3:
        nr_bits = 16
    if codon_len == 4:
        nr_bits = 22

    index = 0
    base_10_numbers = []

    for i in range(0, int(len(b47s)/codon_len), 3):

        base_10_val = 0
        
        for j  in range(0, codon_len, 1):

            b47_nr = b47s[index]

            base_10_val = base_10_val + b47_nr*47**(codon_len-j-1)

            index = index + 1
        
        base_10_numbers.append(base_10_val)

    bit_string = ""

    for i in range(0, len(base_10_numbers), 1):
        
        decimal = base_10_numbers[i]

        bit_len = ""

        for i in range(nr_bits):
            bit = str((decimal // (2**i)) % 2)

            bit_len = bit + bit_len

        bit_string = bit_string + bit_len

    return bit_string


def file_to_binary_string(file_path):
    with open(file_path, "rb") as file:
        binary_code = file.read()
        binary_string = "".join(format(byte, "08b") for byte in binary_code)
    return binary_string


def binary_string_to_file(binary_string, file_path):
    with open(file_path, "wb") as file:
        bytes_list = [
            int(binary_string[i : i + 8], 2) for i in range(0, len(binary_string), 8)
        ]
        bytes_arr = bytearray(bytes_list)
        file.write(bytes_arr)


def DNA_strand_to_file(strand, file_path):
    with open(file_path, "w") as file:
        file.write(strand)
