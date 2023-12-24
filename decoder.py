import galois
from config import codons, codon_dict
import numpy as np
from util import (
    flip_matrix,
    binary_string_to_file,
    base47_to_bin,
    GF_to_ints,
    complementary_strand,
    reversed_strand,
    pair_columns,
    separate_columns
)


class Decoder:
    def __init__(self):
        print("Decoder")

        self.k = 40
        self.m = 37
        self.red_frac_A = 0.1
        self.redB = 6

        self.GF = galois.GF(47)
        self.GF2 = galois.GF(47**2)

        self.complementary = False
        self.bases_per_b47_nr = 3
        self.corr_strand_length = self.k * self.bases_per_b47_nr + self.redB * self.bases_per_b47_nr

        self.rs_col = galois.ReedSolomon(46, 46 - self.redB, field=self.GF)
        self.rs_row = galois.ReedSolomon(
            47**2 - 1, int((47**2 - 1) * self.red_frac_A), field=self.GF2
        )

    def decode(self, inputPath, outputPath):
        print("Decoding...")

        # Read in the DNA strand and convert to strand with values base 47
        Base47_strands = self.read_strands(inputPath)

        if self.complementary == True:
            matched_strands, unmatched_strands = self.DNA_matching(Base47_strands)
            ins_del_corrected_strands, non_corrected_strands = self.insert_delet_compensation(matched_strands)

        # Reed Solomon decoding along colums and index
        decoded_cols, indexes, column_errors = self.RS_col_decoder(Base47_strands)

        # Sort the columns by index
        sorted_cols = self.sort_columns_on_index(decoded_cols, indexes)

        # Convert the galois field back to an array of integers
        sorted_cols = GF_to_ints(sorted_cols)

        # Reed Solomon decoding along rows
        res_cols, row_errors = self.RS_row_decoder(sorted_cols)

        # Concatenate all the DNA-strands
        total = []
        for col in res_cols:
            total += col

        # Convert list of values base 47 to binary string
        bits = base47_to_bin(total)
        binary_string_to_file(bits, outputPath)

        print("Done!")


    def read_strands(self, inputPath):
        with open(inputPath, "r") as file:
            DNA_data = file.read().splitlines()
        file.close()

        Base47_strands = []

        if self.complementary == False:
            for strand in DNA_data:
                if len(strand) != 46 * 3:
                        continue
                compl_strand = complementary_strand(strand) 
                # ---> Hier wordt de complementaire nu gesynthetiseed na sequencing ipv dat de complementaire daarvoor al bestaat?

                # Convert the DNA strand to base 47
                Base47_strands.append(self.DNA_to_base47(strand))
                Base47_strands.append(self.DNA_to_base47(compl_strand))
                Base47_strands.append(self.DNA_to_base47(reversed_strand(strand)))
                Base47_strands.append(self.DNA_to_base47(reversed_strand(compl_strand)))

            return Base47_strands

        elif self.complementary == True:
            return DNA_data
        

    def DNA_to_base47(self, DNA_data):
        bases47 = []
        for i in range(0, len(DNA_data), 3):
            codon = DNA_data[i : i + 3]
            if codon not in codons:
                bases47.append(0)
                # print("False codon")
                continue
            bases47.append(codon_dict[codon])

        return bases47
    

    def DNA_matching(self, DNA_strings):
        # Assumes DNA_strings is a list containing each sequenced nucleotide string. Can be changed if necessary.
    
        newly_sorted = []
        not_matched = []
        strings_copy = DNA_strings.copy()

        # Keep matching until either all complements are found or all (remaining) strands are determined unmatched for current error rate.
        while len(strings_copy) > 0:

            # Selects DNA string to find the complement to. Removes it from the list as to not perfomr unnecessary double checks.
            DNA_seq = strings_copy.pop(0)
        
            # Creates the exact complement to the DNA sequence.
            complement_seq = complementary_strand(DNA_seq)

            # Finds the (closest) complement in the DNA list.
            # --> Here does not yet specify a length for or a max error for the initial matching. Might need to be added.
            found_complement = self.find_complementary(complement_seq, strings_copy)

            # If a complement was found puts it with the DNA strand into a list and removes the complement from the still to be matched list.
            if found_complement != False:
                newly_sorted.append([DNA_seq, found_complement])
                strings_copy.remove(found_complement)
        
            # If no complement was found in the DNA list the DNA strand is added to a list for further inspection.
            else:
                not_matched.append(DNA_seq)

        no_existing_comp = []

        # If no complement was found during initial search it now does an addtional search on the remaining unmatched strands for
        # a higher allowed error margin.
        while len(not_matched) > 0:
        
            DNA_seq = not_matched.pop(0)
            complement_seq = complementary_strand(DNA_seq)

            # Here does a full (True) Levenshtein check between each strand with an allowed error margin of 20 (for now random value)
            found_complement = self.find_complementary(complement_seq, not_matched, True, 20)

            # If a match was found for these assumptions:
            if found_complement != False:
                newly_sorted.append([DNA_seq,found_complement])
                not_matched.remove(found_complement)
            
            # If still no match was found: (probably assume complementary strand was lost)
            else:
                no_existing_comp.append(DNA_seq)

        # Returns a list of lists with DNA strands and their complementary, and a list of unmatched sequences.
        return newly_sorted, no_existing_comp
    

    def find_complementary(self, DNA_strand, DNA_list, full_check=False, max_comp_error=7, comp_len=20):
        # Takes as inputs:
        # - DNA_strand to which to find the complement.
        # - DNA_list with all other DNA strands.
        # - comp_len the length of strand for first comparison.
        # - max_comp_error the max error allowed in the first comparison.

        comp_vals = [0]*len(DNA_list)

        # Performs the levenshtein distance calculation between the target DNA and the DNA list
        # where it only compares the first X nucleotides to save on amount of computation. 
        # ---> Standard set as first 20 bases to be compared. Is random atm, so can easily be changed.
        for i in range(0,len(DNA_list),1):

            comp_strand = DNA_list[i]

            if full_check ==True:
                comp_vals[i] = self.levenshtein_distance(DNA_strand,comp_strand)
            else:
                comp_vals[i] = self.levenshtein_distance(DNA_strand[:comp_len],comp_strand[:comp_len])

        min_diff = min(comp_vals)

        # Checks if the found minimal error is smaller than the given maximal allowed error.
        # If not, assumes no clear complement present at current time.
        if min_diff <= max_comp_error:

            # Checks if there is only one found complement. If not, performs levenshtein distance
            # calculation for all found complements with the minimal error to determine which of those has minimal error.    
            if comp_vals.count(min_diff) > 1 and full_check == False:

                # Find indexes of all found complements with minimal error.
                index_list = [i for i, val in enumerate(comp_vals) if val == min_diff]
                #https://www.geeksforgeeks.org/python-ways-to-find-indices-of-value-in-list/

                comp_vals_long = [0]*len(index_list)

                # Perform full levenshtein for all minimal error found complements
                for i in range(0,len(index_list),1):
                    comp_vals_long[i] = self.levenshtein_distance(DNA_strand,DNA_list[index_list[i]])

                # Find index of the strand with minimal error and select complement
                comp_index = index_list[comp_vals_long.index(min(comp_vals_long))]
                complement = DNA_list[comp_index]

            # If there is only one found complement with minimal error this is the complement.
            # If full_check = True and its still has two strands that have the same error rate it is assumed that two identical strands
            # were synthesized, so two identical complements exist and it doesn't matter which is selected here.
            else:
                complement = DNA_list[comp_vals.index(min_diff)]

        # If none of the DNA strands have an error smaller than the given max allowed error it is
        # assumed that the complementary strand is either very damaged or lost.
        else:
            # Assumes no complement can be found (at present time) due to either to much damage or the complement being lost.
            complement = False

        # Returns the closest match from the DNA list to the generated exact DNA string.
        return complement


    def levenshtein_distance(a, b, error_check=False):
        # Performs levenshtein distance calculation between two strings a and b.
        # If error_check is set to True, make sure that "a" is the string of correct length and "b" is the string of wrong legth.
        # Good explanation: https://blog.paperspace.com/measuring-text-similarity-using-levenshtein-distance/
    
        LS_matrix = np.zeros((len(b)+1,len(a)+1))
        len_a = len(a)
        len_b = len(b)
    
        for j in range(0, len_a+1, 1):
            LS_matrix[0,j] = j
            
        for k in range(0, len_b+1, 1):
            LS_matrix[k,0] = k
    
        for l in range(0, len_b, 1):
            for m in range(0, len_a, 1):
                vals = [LS_matrix[l,m],LS_matrix[l+1,m],LS_matrix[l,m+1]]
            
                if a[m] == b[l]:
                    LS_matrix[l+1,m+1] = vals[0]
                else:
                    LS_matrix[l+1,m+1] = min(vals) + 1

        # Performs an error correction for insertions and deletions if one string is of correct length and one string is of wrong length.
        if error_check == True:
            errors = []
            index_row = 1
            index_col = 1
            errors_present = True
        
            # Checks which error occurs at which position. At the moment also indicates a "subs" in cases no error occured
            # at a certain position. For the purposes of our use it is not important or necessary to fix this.
            while errors_present == True:
            
                index_a = len_a - index_col
                index_b = len_b - index_row
                error_check_vals = [LS_matrix[index_b+1, index_a], LS_matrix[index_b, index_a+1], LS_matrix[index_b, index_a]]
            
                if error_check_vals[2] == min(error_check_vals):
                    errors.append((index_b,"subs"))
                    index_col = index_col+1
                    index_row = index_row+1
                
                elif error_check_vals[1] == min(error_check_vals):
                    errors.append((index_b,"ins"))
                    index_row = index_row+1
            
                elif error_check_vals[0] == min(error_check_vals):
                    errors.append((index_b,"del"))
                    index_col = index_col+1
                
                # If the error value drops to zero there are no more errors present between the two strings,
                # so the loop can be ended early.
                if min(error_check_vals) == 0:
                    errors_present = False
        
            # Cycles through the errors to check which errors are insertions and which are deletions and counters them.
            for j in range(0,len(errors),1):
                error_nr = errors[j]

                # If the error found is an insertion, removes a base at the right index.
                if error_nr[2] == "ins":
                    split_b = [*b]
                    del_index = error_nr[0]
                    del split_b[del_index]
                    b = ''.join(split_b)

                # If the error found is an deletion, inserts an A at the right index. (can be any base, so can easily be changed)
                elif error_nr[2] == "del":
                    ins_index = error_nr[0]
                    part_1 = b[:ins_index+1]
                    part_2 = b[ins_index+1:]
                    b = part_1 + "A" + part_2
                    
            # If error_check = True, returns the corrected string b.
            return b
    
        # If error_check = False, Returns the Levenshtein distance between the two strings. 
        return int(LS_matrix[len(b),len(a)])
    

    def insert_delet_compensation(self, DNA_strands):
        # Partial insertion and deletion correction using Levenshtein.
        # - DNA_strands a list of lists of DNA_strands with matched complementary strands.
        
        cor_len = self.corr_strand_length
        corrected_strands_list = []
        uncorrected_strands_list = []

        for i in range(0,len(DNA_strands),1):
            strand_pair = DNA_strands[i]
            
            # If one of the two strands is of the correct length the other length can be adjusted at the right
            # locations where insertion or deletion occured.
            if len(strand_pair[0]) == cor_len and len(strand_pair[1]) != cor_len:
                strand_pair[1] = self.levenshtein(strand_pair[0],strand_pair[1],True)
                corrected_strands_list.append(strand_pair)

            elif len(strand_pair[1]) == cor_len and len(strand_pair[0]) != cor_len:
                strand_pair[0] = self.levenshtein(strand_pair[1],strand_pair[0],True)
                corrected_strands_list.append(strand_pair)

            # If both strands are of wrong length this method of error compensation cannot be applied
            elif len(strand_pair[0]) != cor_len and len(strand_pair[1]) != cor_len:
                uncorrected_strands_list.append(strand_pair)
            
            # If both strands are of correct length nothing has to be done as of yet.
            else:
                corrected_strands_list.append(strand_pair)

        return corrected_strands_list, uncorrected_strands_list

    def RS_col_decoder(self, strand_base47):
        decoded_cols = []
        index_list = []
        for col in strand_base47:
            if col[0] != 1:
                continue
            decoded_col, errors = self.rs_col.decode(col, errors=True)
            decoded_cols.append(decoded_col[3:])
            index_list.append(self.get_index_from_base47(decoded_col[:3]))
        return decoded_cols, index_list, errors


    def RS_row_decoder(self, columns):
        paired_cols = pair_columns(columns)

        self.update_row_rs_field(len(paired_cols))

        encoded_rows = flip_matrix(paired_cols)
        result_rows, errors = self.rs_row.decode(encoded_rows, errors=True)
        decoded_cols = flip_matrix(result_rows)

        res_cols = separate_columns(decoded_cols)

        return res_cols, errors

    def get_index_from_base47(self, base47_list):
        id = base47_list[1] * 47 + base47_list[2]
        return int(id)

    def sort_columns_on_index(self, data, indexes):
        return [x for _, x in sorted(zip(indexes, data))]

    def update_row_rs_field(self, nr_cols):
        red = int(round(nr_cols * self.red_frac_A))

        self.rs_row = galois.ReedSolomon(
            47**2 - 1, (47**2 - 1) - red, field=self.GF2
        )
