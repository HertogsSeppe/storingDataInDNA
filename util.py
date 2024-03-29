from config import codon_dict, codons
import hashlib

""" 
Stack overflow, Brunão: https://stackoverflow.com/questions/75611160/
converting-any-file-to-binary-1-and-0-and-then-back-to-original-file-without-cor
"""


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


def base47_to_bin(data, codon_len=4, m=23):
    if codon_len == 3:
        nr_bits = 16
    if codon_len == 4:
        nr_bits = 22

    col_padding = data[0]
    padding = data[1]

    if sum(data[-m:]) == 0:
        data = data[:-m]

    data = data[2:-col_padding]

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


def DNA_to_base47(DNA_data):
    bases47 = []
    for i in range(0, len(DNA_data), 3):
        codon = DNA_data[i : i + 3]
        if codon not in codons:
            bases47.append(-1)
            # print("False codon")
            continue
        bases47.append(codon_dict[codon])

    return bases47


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


def complementairy_strand(strand):
    compl_strand = (
        strand.replace("G", "c").replace("C", "g").replace("A", "t").replace("T", "a")
    )
    return compl_strand.upper()[::-1]


def compare_files(file1, file2):
    # Compare file1 and file2 using SHA256 checksum, if they are the same return True else return False
    with open(file1, "rb") as fi:
        data1 = fi.read()
    hash1 = hashlib.sha256(data1).hexdigest()
    with open(file2, "rb") as fi:
        data2 = fi.read()
    hash2 = hashlib.sha256(data2).hexdigest()
    return hash1 == hash2
