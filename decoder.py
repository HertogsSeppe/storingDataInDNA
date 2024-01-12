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
)


class Decoder:
    def __init__(self):
        print("Decoder")
        # Redundancy in the row direction
        self.redA = 30
        self.red_frac_A = 0
        self.k = (47**2 - 1) - self.redA

        # Redundancy in the column direction and index length
        self.redB = 20
        self.index_len = 4
        self.m = 46 - self.index_len - self.redB

        # Initialize galois fields
        self.GF = galois.GF(47)
        self.GF2 = galois.GF(47**2)

        # Initialize Reed Solomon encoders with the right redundancy
        self.rs_col = galois.ReedSolomon(46, 46 - self.redB, field=self.GF)
        self.rs_row = galois.ReedSolomon(
            47**2 - 1, (47**2 - 1) - self.redA, field=self.GF2
        )

    def decode(self, inputPath, outputPath):
        # Read in the DNA strand and convert to strand with values base 47
        DNA_strands = self.read_strands(inputPath)

        # Reed Solomon decoding along colums and index
        decoded_cols, indexes, column_errors = self.RS_col_decoder(DNA_strands)
        # print("Col errors:")
        # print(column_errors)

        # Sort the columns by index
        sorted_cols = self.sort_columns_on_index(decoded_cols, indexes)

        # Convert the galois field back to an array of integers
        sorted_cols = GF_to_ints(sorted_cols)

        res_cols = []
        row_errors = []
        succes = []

        # Decode each block in the row direction
        for i in range(len(sorted_cols) // (2 * (47**2 - 1)) + 1):
            # Reed Solomon decoding along rows
            block = sorted_cols[i * 2 * (47**2 - 1) : (i + 1) * 2 * (47**2 - 1)]
            res_block_cols, row_block_errors = self.RS_row_decoder(block)
            # print("Row errors:")
            # print(row_errors)

            # Succes if all the rows have succesfully been decoded
            succes.append(-1 not in row_block_errors)
            # Add the decoded columns to the result columns
            res_cols += res_block_cols
            # Store the errors
            row_errors.append(row_block_errors)

        # Concatenate all the DNA-strands
        total = []
        for col in res_cols:
            total += col

        # Convert list of values base 47 to binary string
        bits = base47_to_bin(total, codon_len=4, m=self.m)
        binary_string_to_file(bits, outputPath)

        # print(succes)

        return column_errors, row_errors, succes

    def set_column_redundancy(self, red):
        self.redB = red
        self.m = 46 - self.index_len - self.redB
        self.rs_col = galois.ReedSolomon(46, 46 - self.redB, field=self.GF)

    def set_row_redundancy(self, red):
        self.redA = red
        self.k = (47**2 - 1) - self.redA
        self.rs_row = galois.ReedSolomon(
            47**2 - 1, (47**2 - 1) - self.redA, field=self.GF2
        )

    def read_strands(self, inputPath):
        # Read in all the strands and generate the complementairy strands
        with open(inputPath, "r") as file:
            DNA_data = file.read().splitlines()
        file.close()

        DNA_strands = []
        for strand in DNA_data:
            compl_strand = complementairy_strand(strand)
            DNA_strands.append(strand)
            # DNA_strands.append(compl_strand)

        return DNA_strands

    def DNA_to_base47(self, DNA_data):
        # Convert DNA strand to a list of base 47 numbers and count number of false codons
        bases47 = []
        false_codons = 0
        for i in range(0, len(DNA_data), 3):
            codon = DNA_data[i : i + 3]
            if codon not in codons:
                bases47.append(0)
                false_codons += 1
                continue
            bases47.append(codon_dict[codon])

        return bases47, false_codons

    def RS_col_decoder(self, DNA_strands):
        # Decode the DNA strands in the column direction
        decoded_cols = []
        index_list = []
        errors = []
        for strand in DNA_strands:
            # If the DNA strands is missing one base, it is assumed that there is a deletion
            if len(strand) == 46 * 3 - 1:
                decoded_col, error, false_bases = self.del_error_correction(strand)
            # If the DNA strands is one base longer, it is assumed that there is an addition
            elif len(strand) == 46 * 3 + 1:
                decoded_col, error, false_bases = self.add_error_correction(strand)
            # If the DNA strands is the correct lenth, it is assumed that there is no deletion or addition
            elif len(strand) == 46 * 3:
                col, false_bases = self.DNA_to_base47(strand)
                decoded_col, error = self.rs_col.decode(col, errors=True)
            else:
                continue

            # if decoded_col[0] != 1:
            #     continue

            if false_bases > self.redB or error == -1:
                continue

            index = self.get_index_from_base47(decoded_col[: self.index_len])

            if index > len(DNA_strands) * 1.2:
                continue

            errors.append(error)
            decoded_cols.append(decoded_col[self.index_len :])
            index_list.append(index)

        return decoded_cols, index_list, errors

    def RS_row_decoder(self, columns):
        paired_cols = pair_columns(columns)

        encoded_rows = flip_matrix(paired_cols)
        result_rows, errors = self.rs_row.decode(encoded_rows, errors=True)
        decoded_cols = flip_matrix(result_rows)

        res_cols = separate_columns(decoded_cols)

        return res_cols, errors

    def get_index_from_base47(self, base47_list, index_len=4):
        id = 0

        for i in range(index_len):
            id += int(base47_list[i]) * (47 ** (index_len - 1 - i))

        return id

    def sort_columns_on_index(self, data, indexes):
        sorted_data = []

        for i in range(max(indexes)):
            data_index = indexes.index(i) if i in indexes else -1

            if data_index == -1:
                sorted_data.append([0] * self.m)
                continue

            sorted_data.append(data[data_index])
        return sorted_data

    def update_row_rs_field(self, nr_cols):
        red = int(round(nr_cols * self.red_frac_A))
        self.redA = red

        self.rs_row = galois.ReedSolomon(
            47**2 - 1, (47**2 - 1) - self.redA, field=self.GF2
        )

    def del_error_correction(self, strand):
        for i in range(45):
            added_strand = strand[: 3 * i] + "A" + strand[3 * i :]
            base_47, false_bases = self.DNA_to_base47(added_strand)
            col, err = self.rs_col.decode(base_47, errors=True)

            if err != -1:
                return col, err, false_bases
        return col, err, false_bases

    def add_error_correction(self, strand):
        for i in range(45):
            subtr_strand = strand[: 3 * i] + strand[3 * i + 1 :]
            base_47, false_bases = self.DNA_to_base47(subtr_strand)
            col, err = self.rs_col.decode(base_47, errors=True)

            if err != -1:
                return col, err, false_bases
        return col, err, false_bases
