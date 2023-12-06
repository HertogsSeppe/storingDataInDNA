from util import file_to_binary_string, DNA_srand_to_file
from config import codons


class Encoder:
    def __init__(self):
        print("Encoder")

    def encode(self, inputPath, outputPath):
        binaryString = file_to_binary_string(inputPath)
        dnaStrand = self.bits_to_bases(binaryString, 4)
        DNA_srand_to_file(dnaStrand, outputPath)
        # print(len(binaryString), len(dnaStrand), len(binaryString) / len(dnaStrand))

    def bits_to_bases(self, bits, nr_codons=3):
        if nr_codons == 3:
            wl = 16
        elif nr_codons == 4:
            wl = 22
        else:
            print("ERROR: Invalid sequence")
            return ""

        zeroPad = (wl - len(bits) % wl) % wl
        bits = bits + zeroPad * "0"

        id_secuence = self.generate_index(id=0, padding=zeroPad)

        strand = id_secuence

        for i in range(len(bits) // wl):
            word = int(bits[wl * i : wl * i + wl], 2)

            strand += self.word_to_bases(word, nr_codons)

        return strand

    def generate_index(self, id=0, padding=0):
        id_secuence = codons[padding % 47]

        for i in range(2, -1, -1):
            index = (id // (47**i)) % 47
            id_secuence += codons[index]

        return id_secuence

    def word_to_bases(self, word, nr_codons):
        strand = ""

        for i in range(nr_codons - 1, -1, -1):
            index = (word // (47**i)) % 47

            strand += codons[index]

        return strand
