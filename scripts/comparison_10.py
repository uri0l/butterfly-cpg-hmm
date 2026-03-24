from Bio import SeqIO
import matplotlib.pyplot as plt

def detect_cpg_islands(sequence, window_size):
    cpg_count = 0
    non_cpg_count = 0

    for i in range(0, len(sequence), window_size):
        window = sequence[i:i+window_size]
        cg_count = window.count("CG")
        if cg_count >= window_size // 10:  # 10% threshold
            cpg_count += 1
        else:
            non_cpg_count += 1

    return cpg_count, non_cpg_count

# Load chromosome Z
with open('./chrZ.fasta', 'r') as fp:
    for record in SeqIO.parse(fp,"fasta"):
        seqZ = record.seq
z_cpg, z_non_cpg = detect_cpg_islands(seqZ, 300)
print("Chr Z → CpG:", z_cpg, "| Non-CpG:", z_non_cpg)

# Load chromosome W
with open('./chrW.fasta', 'r') as fp:
    for record in SeqIO.parse(fp,"fasta"):
        seqW = record.seq
w_cpg, w_non_cpg = detect_cpg_islands(seqW, 300)
print("Chr W → CpG:", w_cpg, "| Non-CpG:", w_non_cpg)

# Load chromosome 28
with open('./chr28.fasta', 'r') as fp:
    for record in SeqIO.parse(fp,"fasta"):
        seq28 = record.seq
c28_cpg, c28_non_cpg = detect_cpg_islands(seq28, 300)
print("Chr 28 → CpG:", c28_cpg, "| Non-CpG:", c28_non_cpg)


# COMPARISON OF CPG + NON-CPG REGIONS
chromosomes = ['Z', 'W', '28']
cpg_islands = [z_cpg, w_cpg, c28_cpg]
non_cpg_islands = [z_non_cpg, w_non_cpg, c28_non_cpg]

bar_width = 0.35
x = range(len(chromosomes))

plt.figure(figsize=(8, 6))
plt.bar(chromosomes, cpg_islands, width=bar_width, label='CpG Islands', color='cyan', zorder=2)
plt.bar(chromosomes, non_cpg_islands, width=bar_width, label='Non-CpG Regions', color='gray', zorder=1)

plt.yscale("log")

plt.xlabel('Chromosome')
plt.ylabel('Window Count (log scale)')
plt.title('CpG vs Non-CpG Islands in Chromosomes Z, W, 28 (10% Threshold)')
plt.xticks(x, chromosomes) 
plt.legend()
plt.savefig('cpg_island_comparison_10.png', dpi=300) 
plt.show()
