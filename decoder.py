import galois
from config import codons
from util import (
    flip_matrix,
    binary_string_to_file,
    base47_to_bin,
    GF_to_ints,
    complementairy_strand,
    reversed_strand,
    pair_columns,
    separate_columns,
)


class Decoder:
    def __init__(self):
        print("Decoder")

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

    def decode(self, inputPath, outputPath):
        print("Decoding...")

        # Read in the DNA strand
        DNA_strands = self.read_strands(inputPath)

        # Reed Solomon decoding along colums and index
        decoded_cols, indexes = self.RS_col_decoder(DNA_strands)

        # Sort the columns by index
        sorted_cols = self.sort_columns_on_index(decoded_cols, indexes)

        # Convert the galois field back to an array of integers
        sorted_cols = GF_to_ints(sorted_cols)

        # Reed Solomon decoding along rows
        res_cols = self.RS_row_decoder(sorted_cols)

        # Concatenate all the DNA-strands
        total = []
        for col in res_cols:
            total += col

        # Convert list of values base 47 to binary string
        bits = base47_to_bin(total)
        binary_string_to_file(bits, outputPath)

        print("Done!")

    def read_strands(self, inputPath):
        with open(inputPath, "r") as file:
            DNA_data = file.read().splitlines()
        file.close()

        DNA_strands = []
        for strand in DNA_data:
            if len(strand) != 46 * 3:
                continue
            compl_strand = complementairy_strand(strand)

            # Convert the DNA strand to base 47
            DNA_strands.append(self.DNA_to_base47(strand))
            DNA_strands.append(self.DNA_to_base47(compl_strand))
            DNA_strands.append(self.DNA_to_base47(reversed_strand(strand)))
            DNA_strands.append(self.DNA_to_base47(reversed_strand(compl_strand)))

        return DNA_strands

    def DNA_to_base47(self, DNA_data):
        bases47 = []
        for i in range(0, len(DNA_data), 3):
            codon = DNA_data[i : i + 3]
            if codon not in codons:
                bases47.append(0)
                # print("False codon")
                continue
            # This can be optimised with a dict
            bases47.append(codons.index(codon))

        return bases47

    def RS_col_decoder(self, strand_base47):
        decoded_cols = []
        index_list = []
        for col in strand_base47:
            if col[0] != 1:
                continue
            decoded_col = self.rs_col.decode(col)
            decoded_cols.append(decoded_col[3:])
            index_list.append(self.get_index_from_base47(decoded_col[:3]))
        return decoded_cols, index_list

    def RS_row_decoder(self, columns):
        paired_cols = pair_columns(columns)

        self.update_row_rs_field(len(paired_cols))

        encoded_rows = flip_matrix(paired_cols)
        result_rows = self.rs_row.decode(encoded_rows)
        decoded_cols = flip_matrix(result_rows)

        res_cols = separate_columns(decoded_cols)

        return res_cols

    def get_index_from_base47(self, base47_list):
        id = base47_list[1] * 47 + base47_list[2]
        return int(id)

    def sort_columns_on_index(self, data, indexes):
        return [x for _, x in sorted(zip(indexes, data))]

    def update_row_rs_field(self, nr_cols):
        red = int(round(nr_cols * self.red_frac_A))

        self.rs_row = galois.ReedSolomon(
            47**2 - 1, (47**2 - 1) - red, field=self.GF2
        )
