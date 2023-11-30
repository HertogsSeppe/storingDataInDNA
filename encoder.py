from util import file_to_binary_string, bits_to_bases, DNA_strand_to_file


class Encoder:
    def __init__(self):
        print("Encoder")

    def encode(self, inputPath, outputPath):

        binaryString = file_to_binary_string(inputPath)

        dnaStrand3 = bits_to_bases(binaryString, 3)
        print(len(binaryString), len(dnaStrand3), len(binaryString) / len(dnaStrand3))

        dnaStrand4 = bits_to_bases(binaryString, 4)
        print(len(binaryString), len(dnaStrand4), len(binaryString) / len(dnaStrand4))

        DNA_strand_to_file(dnaStrand3, outputPath)
