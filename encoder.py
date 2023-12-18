import galois

from util import (
    file_to_binary_string,
    DNA_strand_to_file,
    flip_matrix,
    pair_columns,
    separate_columns,
)
from config import codons


class Encoder:
    def __init__(self):
        print("Encoder")
        self.k = 40
        self.m = 37
        self.red_frac_A = 0.1
        self.redB = 6

        self.GF = galois.GF(47)
        self.GF2 = galois.GF(47**2)

        self.rs_col = galois.ReedSolomon(46, 46 - self.redB, field=self.GF)
        self.rs_row = galois.ReedSolomon(
            47**2 - 1, int((47**2 - 1) * self.red_frac_A), field=self.GF2
        )

    def encode(self, inputPath, outputPath):
        # Read in the binary string
        binaryString = file_to_binary_string(inputPath)

        # Convert the binary string to a list of values base 47
        base47_list = self.bits_to_base47(binaryString, 4)

        # Apply reed solomon error correction and cut up in lists of lenth k + redA
        encrypted_data = self.apply_reed_solomon(base47_list)

        # Convert the base 47 lists to dna strands
        dnaStrand = self.base47_to_DNA(encrypted_data)

        DNA_strand_to_file(dnaStrand, outputPath)
        # binary_string_to_file(base47_to_bin(base47_list), outputPath)

        # Bits per base, before reed solomon
        print(
            len(binaryString),
            len(base47_list),
            len(binaryString) / (len(base47_list) * 3),
        )
        # Bits per base, after reed solomon
        print(len(binaryString), len(dnaStrand), len(binaryString) / (len(dnaStrand)))

    def bits_to_base47(self, bits, nr_codons=3):
        if nr_codons == 3:
            wl = 16
        elif nr_codons == 4:
            wl = 22
        else:
            print("ERROR: Invalid sequence")
            return ""

        data_list = []

        if wl == 22:
            # Add zero paddig to the end of the binary string
            zeroPad = (wl - len(bits) % wl) % wl
            bits = bits + zeroPad * "0"

            # Set the first element of the list to the added padding
            data_list = [zeroPad]

        # Go over the binary string in steps of "wl" and convert to a list of base 47 numbers
        for i in range(len(bits) // wl):
            word = int(bits[wl * i : wl * i + wl], 2)
            data_list += self.word_to_base47(word, nr_codons)

        return data_list

    def apply_reed_solomon(self, raw_data):
        # Split the strands into segments of length m
        columns = self.split_strands(raw_data)

        # Pair every two columns for encoding with exension field
        paired_cols = pair_columns(columns)

        # Reed Solomon encoding along rows
        self.update_row_rs_field(len(paired_cols))
        rows = flip_matrix(paired_cols)
        encoded_rows = self.rs_row.encode(rows)

        # Separate the paired columns
        paired_encoded_cols = flip_matrix(encoded_rows)
        separated_cols = separate_columns(paired_encoded_cols)

        # Add index
        for i in range(len(separated_cols)):
            separated_cols[i] = self.generate_index_base47(i) + separated_cols[i]

        # Reed solomon in the column direction
        result = self.rs_col.encode(separated_cols)

        return result

    def split_strands(self, raw_data):
        strands = []
        remainder = len(raw_data) % self.m

        for i in range(len(raw_data) // self.m):
            strand = raw_data[i * self.m : (i + 1) * self.m]
            strands.append(strand)
        if remainder == 0:
            return strands

        # If there is a remainder, add padded list at the end
        strand = raw_data[-remainder:] + (self.m - remainder) * [0]
        strands.append(strand)

        return strands

    def generate_index_base47(self, id=0):
        id_sequence = [1, id // 47, id % 47]
        return id_sequence

    def word_to_base47(self, word, nr_codons):
        strand = []

        for i in range(nr_codons - 1, -1, -1):
            index = (word // (47**i)) % 47

            strand.append(index)

        return strand

    def base47_to_DNA(self, data):
        nucleotides = ""
        for strand in data:
            for index in strand:
                nucleotides += codons[int(index)]
            nucleotides += "\n"
        return nucleotides

    def update_row_rs_field(self, nr_cols):
        red = int(round(nr_cols * self.red_frac_A))

        self.rs_row = galois.ReedSolomon(
            47**2 - 1, (47**2 - 1) - red, field=self.GF2
        )
