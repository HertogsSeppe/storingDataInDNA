from util import file_to_binary_string, bits_to_bases, DNA_srand_to_file


class Encoder:
    def __init__(self):
        print("Encoder")

    def encode(self, inputPath, outputPath):
        binaryString = file_to_binary_string(inputPath)
        dnaStrand = bits_to_bases(binaryString, 3)
        print(len(binaryString), len(dnaStrand), len(binaryString) / len(dnaStrand))
        dnaStrand = bits_to_bases(binaryString, 4)
        print(len(binaryString), len(dnaStrand), len(binaryString) / len(dnaStrand))
        DNA_srand_to_file(dnaStrand, outputPath)
