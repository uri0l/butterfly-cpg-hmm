from Bio import SeqIO
import matplotlib.pyplot as plt

# Assign hidden states to each nucleotide: 'C' for CpG, 'N' for non-CpG
def tag_hidden_states(sequence, cpg_nucleotide_positions):
    tags = ['N'] * len(sequence)
    for index in cpg_nucleotide_positions:
        tags[index] = 'C'
    return tags

# Count transitions between states
def compute_transition_counts(tags):
    transitions = {'C': {'C': 0, 'N': 0}, 'N': {'C': 0, 'N': 0}}
    for i in range(1, len(tags)):
        current_state = tags[i - 1]
        next_state = tags[i]
        transitions[current_state][next_state] += 1
    return transitions

# Normalize transition matrix
def normalize_transitions(transitions):
    for current_state, transitions_from_state in transitions.items():
        total = sum(transitions_from_state.values())
        if total > 0:
            for next_state in transitions_from_state:
                transitions_from_state[next_state] /= total
    return transitions

# Count nucleotide emissions per state
def compute_emissions(sequence, cpg_positions, non_cpg_positions, tags):
    emissions = {'C': {'A': 0, 'C': 0, 'G': 0, 'T': 0}, 'N': {'A': 0, 'C': 0, 'G': 0, 'T': 0}}
    for position in cpg_positions + non_cpg_positions:
        nucleotide = sequence[position]
        if nucleotide not in {'A', 'C', 'G', 'T'}:
            continue
        state = tags[position]
        emissions[state][nucleotide] += 1
    for state in emissions:
        total = sum(emissions[state].values())
        for nt in emissions[state]:
            emissions[state][nt] /= total
    return emissions

# Detect CpG islands using window and CG count threshold (10%)
def detect_cpg_islands(sequence, window_size):
    cpg_count = 0
    non_cpg_count = 0
    cpg_regions = []
    non_cpg_regions = []

    for i in range(0, len(sequence), window_size):
        window = sequence[i:i+window_size]
        cg_count = window.count("CG")
        if cg_count >= window_size // 10:
            cpg_regions.append((i, i+window_size, window))
            cpg_count += 1
        else:
            non_cpg_regions.append((i, i+window_size, window))
            non_cpg_count += 1

    return cpg_count, cpg_regions, non_cpg_count, non_cpg_regions

def get_nucleotide_positions(regions):
    positions = []
    for start, end, region_seq in regions:
        for i in range(len(region_seq)):
            absolute_position = start + i
            positions.append(absolute_position)
    return positions

# Load chromosome sequences
with open('chrZ.fasta', 'r') as file:
    sequenceZ = str(next(SeqIO.parse(file, "fasta")).seq).upper()

with open('chrW.fasta', 'r') as file:
    sequenceW = str(next(SeqIO.parse(file, "fasta")).seq).upper()

with open('chr28.fasta', 'r') as file:
    sequence28 = str(next(SeqIO.parse(file, "fasta")).seq).upper()

window_size = 500

# Z chromosome
Z_cpg_count, Z_cpg_regions, Z_non_count, Z_non_regions = detect_cpg_islands(sequenceZ, window_size)
Z_cpg_positions = get_nucleotide_positions(Z_cpg_regions)
Z_tags = tag_hidden_states(sequenceZ, Z_cpg_positions)
Z_transitions = normalize_transitions(compute_transition_counts(Z_tags))
Z_emissions = compute_emissions(sequenceZ, Z_cpg_positions, get_nucleotide_positions(Z_non_regions), Z_tags)
Z_percent_positions = [(p / len(sequenceZ)) * 100 for p in Z_cpg_positions]
Z_density = Z_cpg_count / (len(sequenceZ) / 1000)
print("Z transition matrix:")
print(Z_transitions)
print("Z emission matrix:")
print(Z_emissions)

# W chromosome
W_cpg_count, W_cpg_regions, W_non_count, W_non_regions = detect_cpg_islands(sequenceW, window_size)
W_cpg_positions = get_nucleotide_positions(W_cpg_regions)
W_tags = tag_hidden_states(sequenceW, W_cpg_positions)
W_transitions = normalize_transitions(compute_transition_counts(W_tags))
W_emissions = compute_emissions(sequenceW, W_cpg_positions, get_nucleotide_positions(W_non_regions), W_tags)
W_percent_positions = [(p / len(sequenceW)) * 100 for p in W_cpg_positions]
W_density = W_cpg_count / (len(sequenceW) / 1000)
print("W transition matrix:")
print(W_transitions)
print("W emission matrix:")
print(W_emissions)

# 28 chromosome
C28_cpg_count, C28_cpg_regions, C28_non_count, C28_non_regions = detect_cpg_islands(sequence28, window_size)
C28_cpg_positions = get_nucleotide_positions(C28_cpg_regions)
C28_tags = tag_hidden_states(sequence28, C28_cpg_positions)
C28_transitions = normalize_transitions(compute_transition_counts(C28_tags))
C28_emissions = compute_emissions(sequence28, C28_cpg_positions, get_nucleotide_positions(C28_non_regions), C28_tags)
C28_percent_positions = [(p / len(sequence28)) * 100 for p in C28_cpg_positions]
C28_density = C28_cpg_count / (len(sequence28) / 1000)
print("28 transition matrix:")
print(C28_transitions)
print("28 emission matrix:")
print(C28_emissions)

# Print densities
print("CpG Island Density per 1000 bp:")
print("Z density:", Z_density)
print("W density:", W_density)
print("28 density:", C28_density)

# Plot normalized positions of CpG islands
plt.figure(figsize=(15, 5))
plt.eventplot([Z_percent_positions, W_percent_positions, C28_percent_positions],
              lineoffsets=[2, 1, 0],
              colors=['blue', 'red', 'green'],
              linelengths=0.8)

plt.yticks([2, 1, 0], ['Chr Z', 'Chr W', 'Chr 28'])
plt.xlabel("Position in Chromosome (% of Length)", fontsize=12)
plt.title("Normalized Distribution of CpG Islands in Chromosomes Z, W, and 28 (10% Threshold)", fontsize=14)
plt.grid(True, axis='x', linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.show()

# Print all matrices again at the end
print("\nNormalized Transition Matrices:")
print("Chr Z:")
print(Z_transitions)
print("Chr W:")
print(W_transitions)
print("Chr 28:")
print(C28_transitions)

print("\nNormalized Emission Matrices:")
print("Chr Z:")
print(Z_emissions)
print("Chr W:")
print(W_emissions)
print("Chr 28:")
print(C28_emissions)

