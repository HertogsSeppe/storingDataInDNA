from util import bases_to_binary, binary_string_to_file

class Decoder:
    def __init__(self):
        print("Decoder")

    def decode(self, inputPath, outputPath):

        with open(inputPath, "r") as file:
            DNA_data = file.read()

        bit_data = bases_to_binary(DNA_data, 3)

        binary_string_to_file(bit_data, outputPath)

        print("Decoding...")
