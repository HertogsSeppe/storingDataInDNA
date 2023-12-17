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
