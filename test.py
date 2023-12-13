from config import codons
import galois

GF = galois.GF(47)
rs_row = galois.ReedSolomon(46, 40, field=GF)
rs_col = galois.ReedSolomon(46, 40, field=GF)

def DNA_to_base47(DNA_data):
    bases47 = []
    for i in range(0, len(DNA_data), 3):
        codon = DNA_data[i: i + 3]
        base47 = codons.index(codon) if codon in codons else None
        bases47.append(base47)
    return bases47

def RS_decoder(strand_base47):
    result = []
    decoded_col = rs_col.decode(strand_base47)
    result.append(decoded_col)
    print(result)
    return(result)

def base47_to_DNA(data):
    nucleotides = ""
    # print(len(data), len(data[0]))
    for strand in data:
        for index in strand:
            nucleotides += codons[int(index)]
        nucleotides += "\n"
    return nucleotides

m = []
for i in range(40):
    m.append(i)

list_enc_base47 = []
enc_base47 = rs_col.encode(m)
list_enc_base47.append(enc_base47)
print(list_enc_base47)

enc_strand = base47_to_DNA(list_enc_base47)
print(enc_strand)
base47 = DNA_to_base47(enc_strand)
# base47.remove(None)
dec_strand = RS_decoder(base47)