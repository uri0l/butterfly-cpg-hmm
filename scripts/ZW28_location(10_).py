import matplotlib.pyplot as plt
from Bio import SeqIO
# Assign hidden states to each nucleotide: 'C' for CpG, 'N' for non-CpG
def tag_hidden_states(sequence, cpg_nucleotide_positions):
    #Fist we put everything as N 
    tags = ['N'] * len(sequence)
    # Replace with 'C' for positions that are CpG islands
    for index in cpg_nucleotide_positions:
        tags[index] = 'C'
    return tags

# Count the number of transitions between hidden states (C→C, C→N...)
def compute_transition_counts(tags):
    transitions = {'C': {'C': 0, 'N': 0}, 'N': {'C': 0, 'N': 0}}
    #Count how many change of states we have
    for i in range(1, len(tags)):
        current_state = tags[i - 1]
        next_state = tags[i]
        transitions[current_state][next_state] += 1
    return transitions

# Convert transition counts into probabilities
def normalize_transitions(transition_counts):
    for current_state, transitions_from_state in transition_counts.items():
        total = sum(transitions_from_state.values())  # Total transitions from the current state
        for next_state in transitions_from_state:
            transitions_from_state[next_state] /= total  # Convert to probability
    return transition_counts


# Count how often each nucleotide appears in each state
def compute_emissions(sequence, cpg_positions, non_cpg_positions, tags):
    emissions = {'C': {'A': 0, 'C': 0, 'G': 0, 'T': 0},
                 'N': {'A': 0, 'C': 0, 'G': 0, 'T': 0}}

    # Go through all nucleotide positions (CpG and non-CpG)
    for position in cpg_positions + non_cpg_positions:
        nucleotide = sequence[position]
        if nucleotide not in {'A', 'C', 'G', 'T'}:
            continue  # Skip any unexpected characters
        state = tags[position]  
        emissions[state][nucleotide] += 1  # Count the nucleotide for the current state

    return emissions

# Detect CpG islands using a fixed-size window and CG count threshold(which we will explain in the presentation why we choose those)
def detect_cpg_islands(sequence, window_size):
    cpg_regions = []
    non_cpg_regions = []

    # Slide through the genome in non-overlapping windows
    for i in range(0, len(sequence), window_size):
        window = sequence[i:i + window_size]
        cg_count = window.count("CG") 

        # Threshold(10% in this case)
        if cg_count >= window_size // 10:
            cpg_regions.append((i, i + window_size, window))
        else:
            non_cpg_regions.append((i, i + window_size, window))

    return cpg_regions, non_cpg_regions


# Get the list of absolute positions (in the full sequence) of nucleotides in each region
def get_nucleotide_positions(regions):
    positions = []
    for start, end, region_seq in regions:
        for i in range(len(region_seq)):
            absolute_position = start + i  
            positions.append(absolute_position)
    return positions

#CHROMOSOME Z
with open('./chrZ.fasta', 'rt') as fp:
    for record in SeqIO.parse(fp,"fasta"):
        seqZ = record.seq

window_size = 300

# Detect CpG and non-CpG regions using the windowed threshold method
Z_cpg_regions, Z_nocpg_regions = detect_cpg_islands(seqZ, window_size)
print("\nNumber of CpG islands (Z):", len(Z_cpg_regions))
print("Number of non-CpG islands(Z):", len(Z_nocpg_regions))

# Get all individual nucleotide positions for CpG and non-CpG regions
Z_cpg_positions = get_nucleotide_positions(Z_cpg_regions)
Z_nocpg_positions = get_nucleotide_positions(Z_nocpg_regions)

# Tag each nucleotide with a state: 'C' or 'N'
tags = tag_hidden_states(seqZ, Z_cpg_positions)

# Count and normalize transition probabilities between states
transitions = compute_transition_counts(tags)
normalized_transitions = normalize_transitions(transitions)

# Count nucleotide emissions in each state
emissions = compute_emissions(seqZ, Z_cpg_positions, Z_nocpg_positions, tags)

print("\nNormalized transitions(Z):")
print(normalized_transitions)

print("\nUnnormalized transition counts:")
print(transitions)

print("\nNormalized emissions (Z):")
normalized_emissions = {}
for state, emission_dict in emissions.items():
    total = sum(emission_dict.values())
    normalized_emissions[state] = {nt: count / total for nt, count in emission_dict.items()}

print(normalized_emissions)

print("\nUnnormalized emission counts:")
print(emissions)

print("\n----------------------------------------------------------------------------")

#CHROMOSOME W
with open('./chrW.fasta', 'rt') as fp:
    for record in SeqIO.parse(fp,"fasta"):
        seqW = record.seq

# Detect CpG and non-CpG regions using the windowed threshold method
W_cpg_regions, W_nocpg_regions = detect_cpg_islands(seqW, window_size)
print("\nNumber of CpG islands (W):", len(W_cpg_regions))
print("Number of non-CpG islands(W):", len(W_nocpg_regions))

# Get all individual nucleotide positions for CpG and non-CpG regions
W_cpg_positions = get_nucleotide_positions(W_cpg_regions)
W_nocpg_positions = get_nucleotide_positions(W_nocpg_regions)

# Tag each nucleotide with a state: 'C' or 'N'
tags = tag_hidden_states(seqW, W_cpg_positions)

# Count and normalize transition probabilities between states
transitions = compute_transition_counts(tags)
normalized_transitions = normalize_transitions(transitions)

# Count nucleotide emissions in each state
emissions = compute_emissions(seqW, W_cpg_positions, W_nocpg_positions, tags)

print("\nNormalized transitions(W):")
print(normalized_transitions)

print("\nUnnormalized transition counts:")
print(transitions)

print("\nNormalized emissions (W):")
normalized_emissions = {}
for state, emission_dict in emissions.items():
    total = sum(emission_dict.values())
    normalized_emissions[state] = {nt: count / total for nt, count in emission_dict.items()}

print(normalized_emissions)

print("\nUnnormalized emission counts:")
print(emissions)

print("\n----------------------------------------------------------------------------")

#CHROMOSOME 28
with open('./chr28.fasta', 'rt') as fp:
    for record in SeqIO.parse(fp,"fasta"):
        seq28 = record.seq

# Detect CpG and non-CpG regions using the windowed threshold method
C28_cpg_regions, C28_nocpg_regions = detect_cpg_islands(seq28, window_size)
print("\nNumber of CpG islands (C28):", len(C28_cpg_regions))
print("Number of non-CpG islands(C28):", len(C28_nocpg_regions))

# Get all individual nucleotide positions for CpG and non-CpG regions
C28_cpg_positions = get_nucleotide_positions(C28_cpg_regions)
C28_nocpg_positions = get_nucleotide_positions(C28_nocpg_regions)

# Tag each nucleotide with a state: 'C' or 'N'
tags = tag_hidden_states(seq28, C28_cpg_positions)

# Count and normalize transition probabilities between states
transitions = compute_transition_counts(tags)
normalized_transitions = normalize_transitions(transitions)

# Count nucleotide emissions in each state
emissions = compute_emissions(seq28, C28_cpg_positions, C28_nocpg_positions, tags)

print("\nNormalized transitions(C28):")
print(normalized_transitions)

print("\nUnnormalized transition counts:")
print(transitions)

print("\nNormalized emissions (C28):")
normalized_emissions = {}
for state, emission_dict in emissions.items():
    total = sum(emission_dict.values())
    normalized_emissions[state] = {nt: count / total for nt, count in emission_dict.items()}

print(normalized_emissions)

print("\nUnnormalized emission counts:")
print(emissions)

print("\n----------------------------------------------------------------------------")


#PLOT CpG ISLAND POSITIONS
offsets = {'Z': 0.5, 'W': 0, '28': -0.5}
colors = {'Z': 'blue', 'W': 'red', '28': 'green'}

plt.figure(figsize=(15, 4))

plt.eventplot(Z_cpg_positions, lineoffsets=offsets['Z'], colors=colors['Z'], linelengths=0.3, label='Chr Z')
plt.eventplot(W_cpg_positions, lineoffsets=offsets['W'], colors=colors['W'], linelengths=0.3, label='Chr W')
plt.eventplot(C28_cpg_positions, lineoffsets=offsets['28'], colors=colors['28'], linelengths=0.3, label='Chr 28')

plt.xlabel('Position in Chromosome', fontsize=14)
plt.yticks([offsets['28'], offsets['W'], offsets['Z']], ['Chr 28', 'Chr W', 'Chr Z'])
plt.title('CpG Island Distribution in Chromosomes Z, W, and 28 (10% Threshold)', fontsize=16)
plt.grid(True, linestyle='--', alpha=0.5)
plt.savefig('location_zw28_10.png', dpi=300)

plt.show()

