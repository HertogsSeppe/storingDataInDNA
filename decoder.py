import galois
from config import codons
from util import (
    flip_matrix,
    b47_to_binary,
    binary_string_to_file,
    base47_to_bin,
    GF_to_ints,
)
import random


class Decoder:
    def __init__(self):
        print("Decoder")

        k = 40
        m = 37
        self.red_frac_A = 0.1
        self.redB = 6

        self.GF = galois.GF(47)
        self.GF2 = galois.GF(47**2)

        self.rs_col = galois.ReedSolomon(46, 46 - self.redB, field=self.GF)
        self.rs_row = galois.ReedSolomon(
            47**2 - 1, int((47**2 - 1) * self.red_frac_A), field=self.GF2
        )

    def decode(self, inputPath, outputPath):
        print("Decoding...")

        # Read in the DNA strand
        with open(inputPath, "r") as file:
            DNA_data = file.read().splitlines()
        file.close()

        # Convert the DNA strand to base 47
        DNA_strands = []
        for strand in DNA_data:
            base47_line = self.DNA_to_base47(strand)
            DNA_strands.append(base47_line)

        # Reed Solomon decoding along colums and index
        decoded_strands = []
        indexes = []
        for strand in DNA_strands:
            decoded_strand = self.rs_col.decode(strand)
            decoded_strands.append(decoded_strand[3:])
            indexes.append(self.get_index_from_base47(decoded_strand[:3]))

        sorted_cols = (
            decoded_strands  # self.sort_columns_on_index(decoded_strands, indexes)
        )

        sorted_cols = GF_to_ints(sorted_cols)

        paired_cols = []
        for i in range(0, len(sorted_cols) - 1, 2):
            # TODO: make work for odd nr cols
            new_col = []
            for j in range(len(sorted_cols[0])):
                new_col.append(sorted_cols[i][j] + 47 * sorted_cols[i + 1][j])
            paired_cols.append(new_col)
        if len(sorted_cols) % 2 == 1:
            paired_cols.append(sorted_cols[-1])

        # Reed Solomon decoding along rows
        self.update_row_rs_field(len(paired_cols))
        encoded_rows = flip_matrix(paired_cols)
        result_rows = self.rs_row.decode(encoded_rows)
        decoded_cols = flip_matrix(result_rows)

        res_cols = []
        for i, col in enumerate(decoded_cols):
            new_col1 = []
            new_col2 = []

            for j in col:
                new_col1.append(j % 47)
                new_col2.append(j // 47)
            res_cols.append(new_col1)
            res_cols.append(new_col2)

        total = []
        for col in res_cols:
            total += col

        bits = base47_to_bin(total)

        binary_string_to_file(bits, outputPath)
        # Convert list of values base 47 to binary string

        print("Decoding...")

    def DNA_to_base47(self, DNA_data):
        bases47 = []
        for i in range(0, len(DNA_data), 3):
            codon = DNA_data[i : i + 3]
            if codon not in codons:
                bases47.append(0)
                print("False codon")
                continue
            # This can be optimised with a dict
            bases47.append(codons.index(codon))

        return bases47

    def RS_col_decoder(self, strand_base47):
        result_col = []
        decoded_col = self.rs_col.decode(strand_base47)
        result_col.append(decoded_col)
        return result_col

    def RS_row_decoder(self, decoded_col):
        result_row = []
        decoded_row = self.rs_row.decode(decoded_col)
        result_row.append(decoded_row)
        return result_row

    def get_index_from_base47(self, base47_list):
        id = base47_list[0] * 47**2 + base47_list[1] * 47 + base47_list[2]
        return int(id)

    def sort_columns_on_index(self, data, indexes):
        return [x for _, x in sorted(zip(indexes, data))]

    def update_row_rs_field(self, nr_cols):
        red = 6  # int(round(nr_cols * self.red_frac_A))

        self.rs_row = galois.ReedSolomon(
            47**2 - 1, (47**2 - 1) - red, field=self.GF2
        )
