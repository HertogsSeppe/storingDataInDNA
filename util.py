from config import codons
import galois


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

# -------------------------------------------------------------------------------------------------------------
# FUNCTIONS BELOW HERE ARE FOR DECODING PROCESS. MAYBE NEED TO BE REORGANIZED OR PUT SOMEWHERE ELSE THAN UTIL?
# -------------------------------------------------------------------------------------------------------------

def DNA_matching(DNA_strings):
    # Assumes DNA_strings is a list containing each sequenced nucleotide string. Can be changed if necessary.
    
    newly_sorted = []
    not_matched = []
    strings_copy = DNA_strings.copy()

    # Keep matching until either all complements are found or all (remaining) strands are determined unmatched.
    while len(strings_copy) > 0:

        # Selects DNA string to find the complement to.
        DNA_seq = strings_copy.pop(0)
        
        # Creates the exact complement to the DNA sequence.
        complement_seq = complement_str(DNA_seq)

        # Finds the (closest) complement in the DNA list.
        # --> Here does not yet specify a length for or a max error for the initial matching. Might need to be added.
        found_complement = find_complementary(complement_seq, strings_copy)

        # If a complement was found puts it with the DNA strand into a list and removes the complement from the still to be matched list.
        if found_complement != False:
            newly_sorted.append([DNA_seq, found_complement])
            strings_copy.remove(found_complement)
        
        # If no complement was found in the DNA list the DNA strand is added to a list for furter inspection.
        else:
            not_matched.append(DNA_seq)

    while len(not_matched) > 0:

        




    return newly_sorted


def levenshtein_distance(a,b):
    # Performs levenshtein distance calculation between two strings a and b.
    # Good explanation: https://blog.paperspace.com/measuring-text-similarity-using-levenshtein-distance/
    
    LS_matrix = np.zeros((len(b)+1,len(a)+1))
    
    for j in range(0, len(a)+1, 1):
        
        LS_matrix[0,j] = j
            
    for k in range(0, len(b)+1, 1):

        LS_matrix[k,0] = k
    
    for l in range(0, len(b), 1):

        for m in range(0, len(a), 1):
            
            vals = [int(LS_matrix[l,m]),int(LS_matrix[l+1,m]),int(LS_matrix[l,m+1])]
            
            if a[m] == b[l]:
                
                LS_matrix[l+1,m+1] = vals[0]
                
            else:
                LS_matrix[l+1,m+1] = min(vals) + 1

    # Returns the Levenshtein distance between the two strings. 
    return int(LS_matrix[len(b),len(a)])


def find_complementary(DNA_strand, DNA_list, comp_len=10, max_comp_error=5):
    # Takes as inputs:
    # - DNA_strand to which to find the complement.
    # - DNA_list with all other DNA strands.
    # - comp_len the length of strand for first comparison.
    # - max_comp_error the max error allowed in the first comparison.

    comp_vals = [0]*len(DNA_list)

    # Performs the levenshtein distance calculation between the target DNA and the DNA list
    # where it only compares the first X nucleotides to save on amount of computation.
    for i in range(0,len(DNA_list),1):

        comp_strand = DNA_list[i]

        comp_vals[i] = levenshtein_distance(DNA_strand[:comp_len],comp_strand[:comp_len])


    min_diff = min(comp_vals)

    # Checks if the found minimal error is smaller than the given maximal allowed error.
    # If not, assumes no clear complement present at current time.
    if min_diff <= max_comp_error:

        # Since at least one strand with less error than the given max allowed error is found
        # it assumes the complement strand present in the DNA_list so can be found.
        found = True

        # Checks if there is only one found complement. If not, performs levenshtein distance
        # calculation for all found complements with the minimal error to determine which of those has minimal error.
        if comp_vals.count(min_diff) > 1:

            # Find indexes of all found complements with minimal error.
            index_list = [i for i, val in enumerate(comp_vals) if val == min_diff]
            #https://www.geeksforgeeks.org/python-ways-to-find-indices-of-value-in-list/

            comp_vals_long = [0]*len(index_list)

            # Perform levenshtein for all minimal error found complements
            for i in range(0,len(index_list),1):
                comp_vals_long[i] = levenshtein_distance(DNA_strand,DNA_list[index_list[i]])

            # Find index of the strand with minimal error and select complement
            comp_index = index_list[comp_vals_long.index(min(comp_vals_long))]
            complement = DNA_list[comp_index]
        
        # If there is only one found complement with minimal error this is the complement.
        else:
            complement = DNA_list[comp_vals.index(min_diff)]

    # If none of the DNA strands have an error smaller than the given max allowed error it is
    # assumed that the complementary strand is either very damaged or lost.
    else:
        # Assumes no complement can be found (at present time) due to either to much damage or the complement being lost.
        complement = False

    # Returns the closest match from the DNA list to the generated exact DNA string.
    return complement


def complement_str(line):
    # Creates the complement for a DNA string.

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
            # Since Illumina sequencing machines can also give an ambiguous result this filters those results out.
                # (e.g. returns "R" if it could be A or G, or returns "Y" if it assumes C or T)
            # For our purposes we assume no ambiguity (i think), so this could be left out.
            complement = complement + "X"

    # Returns the generated complement to the given DNA string.
    return complement


def restructure_B_to_A(base47_list, index_len = 3, redun_B_len = 3):
    # Assumes base47_list is a list of base 47 column lists that have already been correctly ordered by index number
    # and error checked. If not, has to be changed.
    #
    # Assumes all base string have the same length. If not, has to be changed.
    #   - redun_B_len is the vertical (column) redundancy.

    base47_Btrim_list = [0]*len(base47_list)

    # Strips the DNA strings of redundancy numbers and index numbers.
    for i in range(0,len(base47_list),1):

        base47_str = base47_list[i]

        base47_Btrim_list[i] = base47_str[index_len:(len(base47_str)-redun_B_len)]

    col_len = len(base47_Btrim_list[0])
    row_len = len(base47_Btrim_list)

    b47_rows_list = [0]*col_len

    # Cycles through columns selecing appropriate row numbers
    for j in range(0, col_len, 1):

        row = [0]*row_len

        for k in range(0, row_len, 1):

            col = base47_Btrim_list[k]
            row[k] = col[j]
        
        b47_rows_list[j] = row

    # Returns a list the rows of the data.
    return b47_rows_list


def restructure_A_to_final(base47_list, redun_A_len = 3):
    # Assumes base47_list is a list of base 47 row lists that have already been error checked.

    col_len = len(base47_list)
    row_len = len(base47_list[0]-redun_A_len)

    restructured_list = []

    for j in range(0, col_len, 1):
        base47_trimmed_row = base47_list[j][redun_A_len:]

        restructured_list.extend(base47_trimmed_row)

    # Returns a single list containing all the base 47 numbers.
    return restructured_list


def bases_to_47(bases, codon_len=3):
    # Converts the strings of bases into lists of base 47 numbers.
    # Assumes bases is a string of bases, so only processes 1 DNA string at a time. Can be changed if necessary.

    bases = bases[12:] #staat nu voor een index van 3 nrs en met 1 extra bit voor zero-indication. Moet nog variabel / aanpasbaar.
    index = 0
    base_47_numbers = []

    for i in range(0, int(len(bases)/codon_len), 3):

        codon = bases[index: index + codon_len]

        # Checks if the selected codon is valid.
        if codon in codons:
            base_47_numbers.append(codons.index(codon))

        # If the selected codon is not valid:
        else:
            base_47_numbers.append('XXX') #moet nog kijken wat precies te zetten, of om gewoon weg te laten en te tellen als deletie
    
    # Returns a list of base 47 numbers.
    return base_47_numbers


def b47_to_binary(b47s, codon_len=3):
    # Converts a list of base 47 numbers into a string of bits.
    # Assumes b47s is a list of base 47 numbers. if not, can be changed.

    if codon_len == 3:
        nr_bits = 16
    if codon_len == 4:
        nr_bits = 22

    index = 0
    base_10_numbers = []

    # Converts the 3-long base 47 number into decimal.
    for i in range(0, int(len(b47s)/codon_len), 3):

        base_10_val = 0
        
        for j  in range(0, codon_len, 1):

            b47_nr = b47s[index]

            base_10_val = base_10_val + b47_nr*47**(codon_len-j-1)

            index = index + 1
        
        base_10_numbers.append(base_10_val)

    bit_string = ""

    # Converts the decimal number into binary.
    for i in range(0, len(base_10_numbers), 1):
        
        decimal = base_10_numbers[i]

        bit_len = ""

        for i in range(nr_bits):
            bit = str((decimal // (2**i)) % 2)

            bit_len = bit + bit_len

        bit_string = bit_string + bit_len

    # Returns a string with all the bits from the data block.
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
