import galois
from config import codons
from util import flip_matrix, b47_to_binary, binary_string_to_file
import random

k = 40
m = 37
redA = 6
redB = 6


class Decoder:
    def __init__(self):
        print("Decoder")
        self.GF = galois.GF(47)
<<<<<<< Updated upstream
        self.rs_row = galois.ReedSolomon(46, 46 - redA, field=self.GF)
        self.rs_col = galois.ReedSolomon(46, 46 - redB, field=self.GF)
    
=======
        self.rs_row = galois.ReedSolomon(46, 46 - redA, field=self.GF)
        self.rs_col = galois.ReedSolomon(46, 46 - redB, field=self.GF)

>>>>>>> Stashed changes
    def decode(self, inputPath, outputPath):
        print("Decoding...")

        # Read in the DNA strand
        with open(inputPath, "r") as file:
<<<<<<< Updated upstream
            DNA_strand = file.read()
        
        # Convert the DNA strand to base 47
        strand_base47 = self.DNA_to_base47(DNA_strand)

        # Reed Solomon decoding along column and index
        decoded_col = self.RS_col_decoder(strand_base47)

        # Sorting columns on index and make rows
        encoded_rows = self.col_to_row(decoded_col)
        
        # Reed Solomon decoding along rows
        decoded_rows = self.RS_row_decoder(encoded_rows)
=======
            DNA_data = file.readlines()
        file.close()

        # Convert the DNA strand to base 47
        DNA_strands = []
        for line in DNA_data:
            strand = self.DNA_to_base47(line)
            DNA_strands.append(strand)

        # Reed Solomon decoding along colums and index
        decoded_strands = []
        indexes = []
        for strand in DNA_strands:
            decoded_strand = self.rs_col.decode(strand)
            decoded_strands.append(decoded_strand[3:])
            indexes.append(self.get_index_from_base47(decoded_strand[:3]))

        sorted_data = self.sort_columns_on_index(decoded_strands, indexes)

        # Reed Solomon decoding along rows
        sorted_data = flip_matrix(sorted_data)
        result_rows = []
        for row in sorted_data:
            decoded_row = self.rs_col.decode(row)
            result_rows.append(decoded_row)

        result_cols = flip_matrix(result_rows)
>>>>>>> Stashed changes

        total = []
        for col in result_cols:
            total += col

        bits = b47_to_binary(total, 4)

        binary_string_to_file(bits, outputPath)
        # Convert list of values base 47 to binary string

<<<<<<< Updated upstream



        print("Decoding...")

    def DNA_to_base47(self, DNA_strand):
        bases47 = []
        for i in range(0, len(DNA_strand) -3, 3):
            codon = DNA_strand[i: i + 3]
            base47 = codons.index(codon) if codon in codons else None
            bases47.append(base47)
=======
    def DNA_to_base47(self, DNA_data):
        bases47 = []
        for i in range(0, len(DNA_data) - 3, 3):
            codon = DNA_data[i : i + 3]
            if codon not in codons:
                bases47.append(-1)
                print("False codon")
                continue

            # This can be optimised with a dict
            bases47.append(codons.index(codon))

>>>>>>> Stashed changes
        return bases47
    
    def RS_col_decoder(self, strand_base47):
        result_col = []
        decoded_col = self.rs_col.decode(strand_base47)
        result_col.append(decoded_col)
        return(result_col)
        
    def col_to_row(self, encoded_rows):
        print(encoded_rows)

    def RS_row_decoder(self, decoded_col):
        result_row = []
        decoded_row = self.rs_row.decode(decoded_col)
        result_row.append(decoded_row)
        return(result_row)

    def get_index_from_base47(self, base47_list):
        id = base47_list[0] * 47**2 + base47_list[1] * 47 + base47_list[2]
        return int(id)

    def sort_columns_on_index(self, data, indexes):
        return [x for _, x in sorted(zip(indexes, data))]
