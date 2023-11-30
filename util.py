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


""" 
Stack overflow, Brunão: https://stackoverflow.com/questions/75611160/
converting-any-file-to-binary-1-and-0-and-then-back-to-original-file-without-cor
"""

def bases_to_binary(bases, nr_codons=3):
    
    bases = bases[12:]
    index = 0
    base_10_numbers = []

    for i in range(0, int(len(bases)/nr_codons), 3):

        base_10_val = 0
        
        for j  in range(0, nr_codons, 1):

            codon = bases[index: index + nr_codons]

            base_10_val = base_10_val + codons.index(codon)*47**(nr_codons-j-1)

            index = index + nr_codons
        
        base_10_numbers.append(base_10_val)
    
    if nr_codons == 3:
        nr_bits = 16
    if nr_codons == 4:
        nr_bits = 22

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
