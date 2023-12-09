import galois

from util import file_to_binary_string, DNA_srand_to_file, flip_matrix
from config import codons


class Encoder:
    def __init__(self):
        print("Encoder")
        self.GF = galois.GF(47)
        self.rs = galois.ReedSolomon(46, 40, field=self.GF)

    def encode(self, inputPath, outputPath):
        # Read in the binary string
        binaryString = file_to_binary_string(inputPath)

        # Convert the binary string to a list of values base 47
        dnaStrand = self.bits_to_base47(binaryString, 4)

        # Apply reed solomon error correction and cut up in lists of lenth 46
        encrypted_data = self.apply_reed_solomon(dnaStrand)
        # DNA_srand_to_file(dnaStrand, outputPath)
        # print(len(binaryString), len(dnaStrand), len(binaryString) / len(dnaStrand))

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
        strands = self.split_strands(raw_data)
        result = []

        # TODO: Refactor this code into sepparate functions

        for i in range(len(strands) // 40):
            # Flip rows and columns
            matrix = flip_matrix(strands[40 * i : 40 * (i + 1)])
            new_strands = []
            # Reed solomon in the row direction
            for row in matrix:
                codeword = self.rs.encode(row)
                new_strands.append(codeword)
            # Flip back
            new_strands = flip_matrix(new_strands)
            # Reed solomon in the column direction and index
            for j, strand in enumerate(new_strands):
                indexed = self.generate_index_base47(47 * i + j) + strand
                final = self.rs.encode(indexed)
                result.append(final)

        # TODO: deal with the remaining strands

        result = flip_matrix(flip_matrix(result))
        return result

        # print(len(flip_matrix(new_strands)), len(flip_matrix(new_strands)[0]))

    def split_strands(self, raw_data):
        strands = []
        remainder = len(raw_data) % 37

        for i in range(len(raw_data) // 37):
            strand = raw_data[i * 37 : (i + 1) * 37]
            strands.append(strand)
        if remainder == 0:
            return strands

        strand = raw_data[-remainder:] + (37 - remainder) * [0]
        strands.append(strand)

        return strands

    def generate_index_base47(self, id=0):
        id_secuence = []

        for i in range(2, -1, -1):
            index = (id // (47**i)) % 47
            id_secuence.append(index)

        return id_secuence

    def word_to_base47(self, word, nr_codons):
        strand = []

        for i in range(nr_codons - 1, -1, -1):
            index = (word // (47**i)) % 47

            strand.append(index)

        return strand
