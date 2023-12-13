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
        self.rs_row = galois.ReedSolomon(46, (46 - redA), field=self.GF)
        self.rs_col = galois.ReedSolomon(46, (46 - redB), field=self.GF)
    
    def decode(self, inputPath, outputPath):
        # Read in the DNA strand
        with open(inputPath, "r") as file:
            DNA_strand = file.read()
        
        # Convert the DNA strand to base 47
        strand_base47 = self.DNA_to_base47(DNA_strand)

        # Reed Solomon decoding along column and index
        decoded_col = self.RS_col_decoder(strand_base47)

        # Sorting columns on index and make rows
        encoded_rows = self.col_to_row(decoded_col)
        
        # Reed Solomon decoding along rows
        decoded_rows = self.RS_row_decoder(encoded_rows)

        # Convert list of values base 47 to binary string




        print("Decoding...")

    def DNA_to_base47(self, DNA_strand):
        bases47 = []
        for i in range(0, len(DNA_strand) -3, 3):
            codon = DNA_strand[i: i + 3]
            base47 = codons.index(codon) if codon in codons else None
            bases47.append(base47)
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

