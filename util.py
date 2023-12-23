""" 
Stack overflow, Brunão: https://stackoverflow.com/questions/75611160/
converting-any-file-to-binary-1-and-0-and-then-back-to-original-file-without-cor
"""


def file_to_binary_string(file_path):
    # Obtains binary data from a file.

    with open(file_path, "rb") as file:
        binary_code = file.read()
        binary_string = "".join(format(byte, "08b") for byte in binary_code)
    return binary_string


def binary_string_to_file(binary_string, file_path):
    # Writes a string of binary data to a file.

    with open(file_path, "wb") as file:
        bytes_list = [
            int(binary_string[i : i + 8], 2) for i in range(0, len(binary_string), 8)
        ]
        bytes_arr = bytearray(bytes_list)
        file.write(bytes_arr)


def DNA_strand_to_file(strand, file_path):
    # Writes a string of nucleotides to a file.

    with open(file_path, "w") as file:
        file.write(strand)


# Transpose of a matrix


def flip_matrix(matrix):
    new_matrix = []
    for i in range(len(matrix[0])):
        new_matrix.append([int(row[i]) for row in matrix])
    return new_matrix


def GF_to_ints(matrix):
    ints_matrix = [
        [int(matrix[i][j]) for j in range(len(matrix[0]))] for i in range(len(matrix))
    ]
    return ints_matrix


def b47_to_binary(b47s, codon_len=3):
    # assumes b47s is a list of base 47 numbers. if not, can be changed.

    if codon_len == 3:
        nr_bits = 16
    if codon_len == 4:
        nr_bits = 22

    index = 0
    base_10_numbers = []

    for i in range(0, int(len(b47s) / codon_len), 3):
        base_10_val = 0

        for j in range(0, codon_len, 1):
            b47_nr = b47s[index]

            base_10_val = base_10_val + b47_nr * 47 ** (codon_len - j - 1)

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


def base47_to_bin(data, codon_len=4):
    if codon_len == 3:
        nr_bits = 16
    if codon_len == 4:
        nr_bits = 22

    padding = data[0]

    data = data[1:]

    bits = ""

    itters = len(data) - len(data) % codon_len

    for i in range(0, itters, codon_len):
        value = 0

        for j in range(codon_len):
            value += data[i + j] * (47 ** (codon_len - 1 - j))

        bit = str("{0:022b}".format(value))
        bits += bit

    if padding == 0:
        return bits
    bits = bits[:-padding]
    return bits


def pair_columns(columns):

    paired_cols = []
    for i in range(0, len(columns) - 1, 2):
        new_col = []
        for j in range(len(columns[0])):
            new_col.append(columns[i][j] + 47 * columns[i + 1][j])
        paired_cols.append(new_col)
    if len(columns) % 2 == 1:
        paired_cols.append(columns[-1])
    return paired_cols


def separate_columns(columns):
    sep_cols = []
    for col in columns:
        new_col1 = []
        new_col2 = []
        # Separate columns
        for j in col:
            new_col1.append(j % 47)
            new_col2.append(j // 47)
        sep_cols.append(new_col1)
        sep_cols.append(new_col2)
    return sep_cols


def complementary_strand(strand):
    compl_strand = (
        strand.replace("G", "c").replace("C", "g").replace("A", "t").replace("T", "a")
    )
    return compl_strand.upper()


def reversed_strand(strand):
    return strand[::-1]

def levenshtein_distance(a,b, error_check=False):
    # Performs levenshtein distance calculation between two strings a and b.
    # If error_check is set to true, make sure that "a" is the string of correct length and "b" is the string of wrong legth.
    # Good explanation: https://blog.paperspace.com/measuring-text-similarity-using-levenshtein-distance/
    
    LS_matrix = np.zeros((len(b)+1,len(a)+1))
    len_a = len(a)
    len_b = len(b)
    
    for j in range(0, len_a+1, 1):
        LS_matrix[0,j] = j
            
    for k in range(0, len_b+1, 1):
        LS_matrix[k,0] = k
    
    for l in range(0, len_b, 1):
        for m in range(0, len_a, 1):
            vals = [LS_matrix[l,m],LS_matrix[l+1,m],LS_matrix[l,m+1]]
            
            if a[m] == b[l]:
                LS_matrix[l+1,m+1] = vals[0]
            else:
                LS_matrix[l+1,m+1] = min(vals) + 1

    # Performs an error correction for insertions and deletions if one string is of correct length and one string is of wrong length.
    if error_check == True:
        errors = []
        index_row = 1
        index_col = 1
        
        # Checks which error occurs at which position. At the moment also indicates "subs" if no error occures
        # at a certain position. For the purposes of our use it is not important or necessary to fix this.
        for i in range(1,len(b)+1,1):
            
            index_a = len_a - index_col
            index_b = len_b - index_row
            error_check_vals = [LS_matrix[index_b+1, index_a], LS_matrix[index_b, index_a+1], LS_matrix[index_b, index_a]]
            
            if error_check_vals[2] == min(error_check_vals):
                errors.append((index_b,index_a,"subs"))
                index_col = index_col+1
                index_row = index_row+1
                
            elif error_check_vals[1] == min(error_check_vals):
                errors.append((index_b,index_a,"del"))
                index_row = index_row+1
            
            elif error_check_vals[0] == min(error_check_vals):
                errors.append((index_b,index_a,"ins"))
                index_col = index_col+1
        
        # Cycles through the errors to check which errors are insertions and which are deletions and counters them.
        for j in range(0,len(errors),1):
            error_nr = errors[j]
            
            if error_nr[2] == "del":
                split_b = [*b]
                del_index = error_nr[0]
                del split_b[del_index]
                b = ''.join(split_b)

            elif error_nr[2] == "ins":
                ins_index = error_nr[0]
                part_1 = b[:ins_index+1]
                part_2 = b[ins_index+1:]
                b = part_1 + "A" + part_2
        # If error_check = True, returns the corrected string b.
        return b
    
    # If error_check = False, Returns the Levenshtein distance between the two strings. 
    return int(LS_matrix[len(b),len(a)])
