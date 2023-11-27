from util import file_to_binary_string, bits_to_bases, DNA_srand_to_file


class Encoder:
    def __init__(self):
        print("Encoder")

    def encode(self, inputPath, outputPath):
        binaryString = file_to_binary_string(inputPath)
        print(binaryString)
        dnaStrand = bits_to_bases(binaryString, 3)
        print(dnaStrand)
        DNA_srand_to_file(dnaStrand, outputPath)
