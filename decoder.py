import galois
from config import codons

k = 40
m = 37
redA = 6
redB = 6

class Decoder:
    def __init__(self):
        print("Decoder")
        self.GF = galois.GF(47)
        self.rs_row = galois.ReedSolomon(k + redA, k, field=self.GF)
        self.rs_col = galois.ReedSolomon((m + 3 + redB), (m + 3), field=self.GF)
    
    def decode(self, inputPath, outputPath):
        # Read in the DNA strand
        with open(inputPath, "r") as file:
            DNA_data = file.read()
        
        # Convert the DNA strand to base 47
        data_base47 = self.DNA_to_base47(DNA_data)

        # Reed Solomon decoding along colums and index

        
        # Reed Solomon decoding along rows


        # Convert list of values base 47 to binary string




        print("Decoding...")

    def DNA_to_base47(self, DNA_data):
        bases47 = []
        for i in range(0, len(DNA_data), 3):
            codon = DNA_data[i: i + 3]
            base47 = codons.index(codon) if codon in codons else print("False codon")
            bases47.append(base47)
        return bases47

