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

    def bits_to_base47(self, bits, nr_codons=3):
        if nr_codons == 3:
            wl = 16
        elif nr_codons == 4:
            wl = 22
        else:
            print("ERROR: This number of codons is not supported")
            return ""

        # Reserve space for the padding of the last column (calculated in the split_strands function)
        data_list = [0]

        if wl == 22:
            # Add zero paddig to the end of the binary string
            zeroPad = (wl - len(bits) % wl) % wl
            bits = bits + zeroPad * "0"

            # Set the second element of the list to the added padding
            data_list.append(zeroPad)

        # Go over the binary string in steps of "wl" and convert to a list of base 47 numbers
        for i in range(len(bits) // wl):
            word = int(bits[wl * i : wl * i + wl], 2)
            data_list += self.word_to_base47(word, nr_codons)

        return data_list

    def apply_reed_solomon(self, raw_data):
        # Split the strands into segments of length m
        columns = self.split_strands(raw_data)

        message_len = 2 * ((47**2 - 1) - self.redA)
        separated_cols = []

        for i in range(len(columns) // message_len + 1):
            # Pair every two columns for encoding with exension field
            paired_cols = pair_columns(columns[i * message_len : (i + 1) * message_len])

            # Reed Solomon encoding along rows
            # self.update_row_rs_field(len(paired_cols))
            rows = flip_matrix(paired_cols)
            encoded_rows = self.rs_row.encode(rows)

            # Separate the paired columns
            paired_encoded_cols = flip_matrix(encoded_rows)
            separated_cols += separate_columns(paired_encoded_cols)

        # Add index
        for i in range(len(separated_cols)):
            separated_cols[i] = self.generate_index_base47(i) + separated_cols[i]

        # Reed solomon in the column direction
        result = self.rs_col.encode(separated_cols)

        return result

    def split_strands(self, raw_data):
        strands = []
        # Calculate the amount of padding codons in the last strand
        remainder = len(raw_data) % self.m

        # Split the data into strands of the message size m
        for i in range(len(raw_data) // self.m):
            strand = raw_data[i * self.m : (i + 1) * self.m]
            strands.append(strand)

        # If there is no remainder, the strands are complete
        if remainder == 0:
            return strands

        # If there is a remainder, add padded list at the end
        strand = raw_data[-remainder:] + (self.m - remainder) * [0]
        strands.append(strand)

        # Store the remainder in the first codon
        strands[0][0] = self.m - remainder

        return strands

    def generate_index_base47(self, id):
        # Convert an id to an array of base 47 numbers based on the index length
        return self.word_to_base47(id, self.index_len)

    def word_to_base47(self, word, nr_codons):
        # Convert a integer to an array of base 47 numbers
        strand = []

        for i in range(nr_codons - 1, -1, -1):
            index = (word // (47**i)) % 47

            strand.append(index)

        return strand

    def base47_to_DNA(self, data):
        # Convert a list of base 47 codons to a DNA strand
        nucleotides = ""
        for strand in data:
            for index in strand:
                nucleotides += codons[int(index)]
            nucleotides += "\n"
        return nucleotides

    def update_row_rs_field(self, nr_cols):
        # Generate a redundancy based on an fraction
        red = int(round(nr_cols * self.red_frac_A))
        self.redA = red

        self.rs_row = galois.ReedSolomon(
            47**2 - 1, (47**2 - 1) - self.redA, field=self.GF2
        )
