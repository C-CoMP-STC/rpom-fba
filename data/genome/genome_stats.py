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

    #the following are molecular weights of each nucelotide (masses in amu). Subrtact a pyrophosphate from each 
    
    mw_PPI = 176.9671
    mw_dATP = 487.15 - mw_PPI
    mw_dCTP = 463.13 - mw_PPI
    mw_dGTP = 503.15 - mw_PPI
    mw_dTTP = 478.14 - mw_PPI
    avogadro = 6.022e23
    
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
    #total counts of pyrophosphate should be all counts minus 1 
    total_nts = total_G + total_C + total_A + total_T
    total_ppi = total_G + total_C + total_A + total_T - 1
    
    #mass given in grams for each component
    mdATP = (total_A * (mw_dATP )) / avogadro
    mdCTP = (total_C * (mw_dCTP )) / avogadro
    mdGTP = (total_G * (mw_dGTP )) / avogadro
    mdTTP = (total_T * (mw_dTTP )) / avogadro
    mPPI = (total_ppi * (mw_PPI)) / avogadro
    total_mass = mdATP + mdCTP + mdGTP + mdTTP + mPPI + (mw_PPI/avogadro) #grams
   
    #in order to have 1g of DNA we need
    conversion_factor = 1/total_mass
    
    #following should be mmole/1gDNA 
    mmdATP = (mdATP * conversion_factor) / (mw_dATP*1000)
    mmdCTP = (mdCTP * conversion_factor) / (mw_dCTP*1000)
    mmdGTP = (mdGTP * conversion_factor) / (mw_dGTP*1000)
    mmdTTP = (mdTTP * conversion_factor) / (mw_dTTP*1000)
    mmPPI = (mPPI * conversion_factor) / (mw_PPI*1000)
    mDNA = 1/mmPPI #milligrams DNA per 1 mmole pyrophosphate

    # Print totals
    print("=====================================")
    print(f"Total G: {total_G} ({100 * total_G / (total_nts):.2f}%)")
    print(f"Total C: {total_C} ({100 * total_C / (total_nts):.2f}%)")
    print(f"Total A: {total_A} ({100 * total_A / (total_nts):.2f}%)")
    print(f"Total T: {total_T} ({100 * total_T / (total_nts):.2f}%)")
    print(f"Total nucleotides: {total_nts}")
    print("=====================================")

    print(f"mmol per 1 gram dna of dATP {mmdATP} mmol/g")
    print(f"mmol per 1 gram dna of dCTP  {mmdCTP} mmol/g")
    print(f"mmol per 1 gram dna of dGTP  {mmdGTP} mmol/g")
    print(f"mmol per 1 gram dna of dTTP  {mmdTTP} mmol/g")
    print(f"mmole pyrophosphate per 1 gDNA = {mmPPI}")
    print(f"grams DNA per mmole PPI = {mDNA/1000}")

if __name__ == "__main__":
    main()
