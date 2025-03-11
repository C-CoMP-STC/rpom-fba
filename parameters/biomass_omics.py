"""
"""
import requests
import concurrent.futures
import pandas as pd
import xmltodict

from tqdm import tqdm
from collections import Counter
from getpass import getpass


ORG_ID = "GCF_000011965"
DNA_SEQUENCES = {"NC_003911": "data/genome/Rpom_chromosome__NC_003911.12.fasta",
                 "NC_006569" : "data/genome/Rpom_megaplasmid__NC_006569.1.fasta"}
TRANSCRIPTOME = "data/clean/omics/rna-abs.csv"
PROTEOME = "data/clean/omics/prot.csv"

def dna():
    total_G = 0
    total_C = 0
    total_A = 0
    total_T = 0
    for file in DNA_SEQUENCES.values():
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

def get_rna_sequence(frame_id, username, password, replicon_sequences):
    s = requests.Session()
    r = s.post("https://websvc.biocyc.org/credentials/login/",
            data={"email": username, "password": password})
    r.raise_for_status()

    r = s.get(f"https://websvc.biocyc.org/getxml?id={ORG_ID}:{frame_id}&detail=full")
    r.raise_for_status()
    dat = xmltodict.parse(r.text)["ptools-xml"]

    # Get replicon, left and right end positions
    replicon = dat["Gene"]["replicon"]["Genetic-Element"]["@frameid"]
    left = int(dat["Gene"]["left-end-position"]["#text"])
    right = int(dat["Gene"]["right-end-position"]["#text"])

    # Get sequence
    seq = replicon_sequences[replicon][left - 1:right].replace("T", "U")
    return seq

def get_rna_sequences(frame_ids, username, password):
    results = {}

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = {
            executor.submit(get_rna_sequence, frame_id, username, password) : frame_id
            for frame_id in tqdm(frame_ids, desc='Submitting...')}

    with tqdm(total=len(frame_ids)) as pbar:
        for future in concurrent.futures.as_completed(futures):
                frame_id = futures[future]
                results[frame_id] = future.result()
                pbar.update(1)
    
    return results

def rna(abundances, frame_ids):
    pass


def protein(abundances, frame_ids):
    pass


def main():
    # Load DNA sequences
    replicon_sequences = {}
    for replicon, file in DNA_SEQUENCES.items():
        with open(file, "r") as f:
            # Skip header
            f.readline()

            # Read sequence
            seq = f.read().replace("\n", "")
            replicon_sequences[replicon] = seq

    dna()
    
    # Load omics data
    transcriptome = pd.read_csv(TRANSCRIPTOME)
    transcriptome = transcriptome[["DSS3_ac_mean_abund", "DSS3_glc_mean_abund", "SPO_ID (ACCESSION)", "frame_id"]]
    proteome = pd.read_csv(PROTEOME)
    proteome = proteome[["DSS3_ac_mean_abund", "DSS3_glc_mean_abund", "SPO_ID (ACCESSION)", "frame_id"]]

    # Prompt for Biocyc username and password
    u = input("Username: ")
    p = getpass("Password: ")

    # Make sure the u/p are correct
    while True:
        try:
            s = requests.Session()
            r = s.post("https://websvc.biocyc.org/credentials/login/",
                    data={"email": u, "password": p})
            r.raise_for_status()
            break
        except requests.exceptions.HTTPError:
            print("Invalid username or password. Please try again.")
            u = input("Username: ")
            p = getpass("Password: ")

    # Get rna and protein sequences
    get_rna_sequence(transcriptome["frame_id"][0], u, p, replicon_sequences)
    rna_sequences = get_rna_sequences(transcriptome["frame_id"], u, p)
    # protein_sequences = get_protein_sequences(proteome["frame_id"])

    rna_acetate = rna(transcriptome["DSS3_ac_mean_abund"], transcriptome["frame_id"])
    rna_glucose = rna(transcriptome["DSS3_glc_mean_abund"], transcriptome["frame_id"])

    protein_acetate = protein(proteome["DSS3_ac_mean_abund"], proteome["frame_id"])
    protein_glucose = protein(proteome["DSS3_glc_mean_abund"], proteome["frame_id"])


if __name__ == "__main__":
    main()
