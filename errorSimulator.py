import random

from util import DNA_strand_to_file


class ErrorSimulator:
    def __init__(self):
        print("ErrorSimulator")

    def induceErrors(
        self,
        inputPath,
        outputPath,
        dropout=0.01,
        insert_rate=0.00042,
        del_rate=0.00188,
        sub_rate=0.00407,
    ):
        self.dropout = dropout
        self.insert_rate = insert_rate
        self.del_rate = del_rate
        self.sub_rate = sub_rate

        # print("Inducing errors...")
        # Read in sequences from inputpath
        DNA_strands = self.read_strands(inputPath)

        # Use sim_seq_simple to induce errors in the strands
        error_strands = self.sim_seq_simple(
            DNA_strands,
            self.dropout,
            self.insert_rate,
            self.del_rate,
            self.sub_rate,
        )

        # Write the strands back into a file
        error_strands = "\n".join(error_strands)
        DNA_strand_to_file(error_strands, outputPath)

    def sim_seq_simple(
        self,
        sequences,
        dropout=0.01,
        insert_rate=0.00042,
        del_rate=0.00188,
        sub_rate=0.00407,
    ):
        # Takes sequences as a list of unique DNA sequences and simulates: synthesis, PCR and NGS
        # Dropout is the dropout rate (i.e. the number of sequences that are not recovered)
        # For baselevel errors we have insert_rate for insertions, del_rate for deletions and sub_rate for substitutions.
        # Default values for base errors are taken from https://www.nature.com/articles/nbt.4079
        s = sequences.copy()
        random.shuffle(s)
        for _ in range(0, int(dropout * len(s))):
            s.pop()

        for i in range(len(s)):
            seq_l = list(s[i])
            mod = 0  # Inserting messes with indexing so we skip over inserted bases using this counter
            for base_i in range(len(seq_l)):
                mutate_rand = random.random()
                if mutate_rand < insert_rate:
                    # Insert random nucleotide after this base
                    seq_l.insert(base_i + mod, random.choice(["A", "C", "T", "G"]))
                    mod += 1
                elif mutate_rand > insert_rate and mutate_rand < (
                    del_rate + insert_rate
                ):
                    # Delete this nucleotide (later)
                    seq_l[base_i + mod] = "_"
                elif mutate_rand > (insert_rate + del_rate) and mutate_rand < (
                    del_rate + insert_rate + sub_rate
                ):
                    # Substitute base
                    seq_l[base_i + mod] = random.choice(["A", "C", "T", "G"])
            if "_" in seq_l:
                for dels in range(seq_l.count("_")):
                    seq_l.remove("_")
            s[i] = "".join(seq_l)
        return s

    def read_strands(self, inputPath):
        with open(inputPath, "r") as file:
            DNA_data = file.read().splitlines()
        file.close()

        return DNA_data
