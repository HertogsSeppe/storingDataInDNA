import galois

from util import file_to_binary_string, DNA_strand_to_file, flip_matrix
from config import codons

k = 40
m = 37
redA = 6
redB = 6


class Encoder:
    def __init__(self):
        print("Encoder")
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

    def encode(self, inputPath, outputPath):
        # Read in the binary string
        binaryString = file_to_binary_string(inputPath)

        # Convert the binary string to a list of values base 47
        base47_list = self.bits_to_base47(binaryString, 4)

        # Apply reed solomon error correction and cut up in lists of lenth k + redA
        encrypted_data = self.apply_reed_solomon(base47_list)

        dnaStrand = self.base47_to_DNA(encrypted_data)

        DNA_strand_to_file(dnaStrand, outputPath)

        print(
            len(binaryString),
            len(base47_list),
            len(binaryString) / (len(base47_list) * 3),
        )
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

        if wl == 16:
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
        columns = self.split_strands(raw_data)
        print("cols:")
        print(len(columns), len(columns[0]))

        paired_cols = []
        for i in range(0, len(columns) - 1, 2):
            # TODO: make work for odd nr cols
            new_col = []
            for j in range(len(columns[0])):
                new_col.append(columns[i][j] + 47 * columns[i + 1][j])
            paired_cols.append(new_col)
        if len(columns) % 2 == 1:
            paired_cols.append(columns[-1])

        print("pairedcols:")
        print(len(paired_cols), len(paired_cols[0]))

        # Reed Solomon decoding along rows
        self.update_row_rs_field(len(paired_cols))
        rows = flip_matrix(paired_cols)
        encoded_rows = self.rs_row.encode(rows)
        print("rows:")
        print(len(rows), len(rows[0]))
        print("encoded rows:")
        print(len(encoded_rows), len(encoded_rows[0]))

        cols = flip_matrix(encoded_rows)
        print("cols:")
        print(len(cols), len(cols[0]))
        indexed_cols = []
        for i, col in enumerate(cols):
            new_col1 = self.generate_index_base47(2 * i)
            new_col2 = self.generate_index_base47(2 * i + 1)

            for j in col:
                new_col1.append(j % 47)
                new_col2.append(j // 47)
            indexed_cols.append(new_col1)
            indexed_cols.append(new_col2)
        print("inexed:")
        print(len(indexed_cols), len(indexed_cols[0]))

        # Reed solomon in the column direction
        result = self.rs_col.encode(indexed_cols)
        print("res:")
        print(len(result), len(result[0]))

        return result

    def split_strands(self, raw_data):
        strands = []
        remainder = len(raw_data) % m

        for i in range(len(raw_data) // m):
            strand = raw_data[i * m : (i + 1) * m]
            strands.append(strand)
        if remainder == 0:
            return strands

        strand = raw_data[-remainder:] + (m - remainder) * [0]
        strands.append(strand)

        return strands

    def generate_index_base47(self, id=0):
        id_sequence = []

        for i in range(2, -1, -1):
            index = (id // (47**i)) % 47
            id_sequence.append(index)

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
