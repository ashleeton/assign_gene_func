from helper_functions import global_alignment, scoring_function_blosum62, scoring_function_simple, local_alignment
# from Bio import SeqIO

# sars2 = SeqIO.read("data/sars_cov_2.fa", "fasta")
# print(len(sars2.seq))
# print(sars2.seq[:50])


# global_alignment("GTA", "TA", scoring_function_blosum62)
local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])

# scores = {}
# for name, accession in accession_codes.items():
#     handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
#     record = SeqIO.read(handle, "genbank")

# spike_dna = sars2.seq[ 21562 : 25384 ]
# print(len(spike_dna))
# print(spike_dna[:3])
# spike_protein = spike_dna.translate()
# print(len(spike_protein))
# print(spike_protein[-5:])