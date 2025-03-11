"""
This script calculates and prints the nucleotide composition of DNA sequences from given FASTA files,
to inform the development of a biomass composition.

It counts the occurrences of each nucleotide (A, T, G, C), and prints the counts and percentages
for each nucleotide in each sequence. It also calculates and prints the total
counts and percentages for all sequences combined.
Constants:
    SEQUENCES (list of str): List of file paths to the FASTA files containing the DNA sequences.
"""

from collections import Counter

SEQUENCES = ["data/genome/Rpom_chromosome__NC_003911.12.fasta",
             "data/genome/Rpom_megaplasmid__NC_006569.1.fasta"]


def main():
    total_G = 0
    total_C = 0
    total_A = 0
    total_T = 0
    for file in SEQUENCES:
        with open(file, "r") as f:
            # Skip header
            f.readline()

            # Read sequence 
            seq = f.read().replace("\n", "")
            print(f"Sequence {file} has length {len(seq)}")
            
            # Count nucleotides in both strands (sense and antisense)
            nucleotide_count = Counter(seq)
            G_nts = nucleotide_count["G"] + nucleotide_count["C"]
            C_nts = nucleotide_count["C"] + nucleotide_count["G"]
            A_nts = nucleotide_count["A"] + nucleotide_count["T"]
            T_nts = nucleotide_count["T"] + nucleotide_count["A"]
            print(f"G: {G_nts} ({100 * G_nts / (2 * len(seq)):.1f}%)")
            print(f"C: {C_nts} ({100 * C_nts / (2 * len(seq)):.1f}%)")
            print(f"A: {A_nts} ({100 * A_nts / (2 * len(seq)):.1f}%)")
            print(f"T: {T_nts} ({100 * T_nts / (2 * len(seq)):.1f}%)")
            
        # Update totals
        total_G += G_nts
        total_C += C_nts
        total_A += A_nts
        total_T += T_nts
    
    # Print totals
    total_nts = total_G + total_C + total_A + total_T
    print("=====================================")
    print(f"Total G: {total_G} ({100 * total_G / (total_nts):.2f}%)")
    print(f"Total C: {total_C} ({100 * total_C / (total_nts):.2f}%)")
    print(f"Total A: {total_A} ({100 * total_A / (total_nts):.2f}%)")
    print(f"Total T: {total_T} ({100 * total_T / (total_nts):.2f}%)")
    print(f"Total nucleotides: {total_nts}")
    print("=====================================")


if __name__ == "__main__":
    main()
