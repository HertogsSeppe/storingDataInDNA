import galois
from config import codons, codon_dict
from util import (
    flip_matrix,
    binary_string_to_file,
    base47_to_bin,
    GF_to_ints,
    complementairy_strand,
    pair_columns,
    separate_columns,
    DNA_to_base47,
)


class Decoder:
    def __init__(self):
        print("Decoder")

        self.k = 26
        self.m = 23
        self.red_frac_A = 0.1
        self.redB = 20

        self.GF = galois.GF(47)
        self.GF2 = galois.GF(47**2)

        self.rs_col = galois.ReedSolomon(46, 46 - self.redB, field=self.GF)
        self.rs_row = galois.ReedSolomon(47**2 - 1, 47**2 - 1 - 20, field=self.GF2)

    def decode(self, inputPath, outputPath):
        print("Decoding...")

        # Read in the DNA strand and convert to strand with values base 47
        DNA_strands = self.read_strands(inputPath)

        # Reed Solomon decoding along colums and index
        decoded_cols, indexes, column_errors = self.RS_col_decoder(DNA_strands)
        print("Col errors:")
        print(column_errors)

        # Sort the columns by index
        sorted_cols = self.sort_columns_on_index(decoded_cols, indexes)

        # Convert the galois field back to an array of integers
        sorted_cols = GF_to_ints(sorted_cols)

        # Reed Solomon decoding along rows
        res_cols, row_errors = self.RS_row_decoder(sorted_cols)
        print("Row errors:")
        print(row_errors)

        succes = -1 not in row_errors

        # Concatenate all the DNA-strands
        total = []
        for col in res_cols:
            total += col

        # Convert list of values base 47 to binary string
        bits = base47_to_bin(total, codon_len=4, m=self.m)
        binary_string_to_file(bits, outputPath)

        print("Done!")
        return column_errors, row_errors, succes

    def read_strands(self, inputPath):
        with open(inputPath, "r") as file:
            DNA_data = file.read().splitlines()
        file.close()

        DNA_strands = []
        for strand in DNA_data:
            compl_strand = complementairy_strand(strand)
            DNA_strands.append(strand)
            DNA_strands.append(compl_strand)

        return DNA_strands

    def DNA_to_base47(self, DNA_data):
        bases47 = []
        for i in range(0, len(DNA_data), 3):
            codon = DNA_data[i : i + 3]
            if codon not in codons:
                bases47.append(0)
                # print("False codon")
                continue
            bases47.append(codon_dict[codon])

        return bases47

    def RS_col_decoder(self, DNA_strands):
        decoded_cols = []
        index_list = []
        errors = []
        for strand in DNA_strands:
            if len(strand) == 46 * 3 - 1:
                decoded_col, error = self.del_error_correction(strand)
            elif len(strand) == 46 * 3 + 1:
                decoded_col, error = self.add_error_correction(strand)
            elif len(strand) == 46 * 3:
                col = DNA_to_base47(strand)
                decoded_col, error = self.rs_col.decode(col, errors=True)
            else:
                continue

            if decoded_col[0] != 1:
                continue

            errors.append(error)
            decoded_cols.append(decoded_col[3:])
            index_list.append(self.get_index_from_base47(decoded_col[:3]))

        return decoded_cols, index_list, errors

    def RS_row_decoder(self, columns):
        paired_cols = pair_columns(columns)

        self.update_row_rs_field(len(paired_cols))

        encoded_rows = flip_matrix(paired_cols)
        result_rows, errors = self.rs_row.decode(encoded_rows, errors=True)
        decoded_cols = flip_matrix(result_rows)

        res_cols = separate_columns(decoded_cols)

        return res_cols, errors

    def get_index_from_base47(self, base47_list):
        id = int(base47_list[1]) * 47 + int(base47_list[2])
        return id

    def sort_columns_on_index(self, data, indexes):
        sorted_data = []

        for i in range(max(indexes)):
            data_index = indexes.index(i) if i in indexes else -1

            if data_index == -1:
                sorted_data.append([0] * self.m)
                print("Missing")
                continue

            sorted_data.append(data[data_index])
        return sorted_data

    def update_row_rs_field(self, nr_cols):
        red = int(round(nr_cols * self.red_frac_A))

        red = 40

        self.rs_row = galois.ReedSolomon(
            47**2 - 1, (47**2 - 1) - red, field=self.GF2
        )

    def del_error_correction(self, strand):
        for i in range(40):
            added_strand = strand[: 3 * i] + "A" + strand[3 * i :]
            base_47 = DNA_to_base47(added_strand)
            col, err = self.rs_col.decode(base_47, errors=True)

            if err != -1:
                return col, err
        return col, err

    def add_error_correction(self, strand):
        for i in range(40):
            subtr_strand = strand[: 3 * i] + strand[3 * i + 1 :]
            base_47 = DNA_to_base47(subtr_strand)
            col, err = self.rs_col.decode(base_47, errors=True)

            if err != -1:
                return col, err
        return col, err
