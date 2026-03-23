import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
from collections import defaultdict

class CpGHMM:
    """Complete Hidden Markov Model for CpG island detection"""
    
    def __init__(self):
        self.states = ['C', 'N']  # CpG island, Non-CpG
        self.nucleotides = ['A', 'C', 'G', 'T']
        self.transition_probs = None
        self.emission_probs = None
        self.initial_probs = None
    
    def train(self, sequence, cpg_positions):
        """Train HMM parameters from labeled data"""
        # Create state sequence
        states = ['N'] * len(sequence)
        for pos in cpg_positions:
            if pos < len(states):
                states[pos] = 'C'
        
        # Calculate transition probabilities
        self.transition_probs = self._calculate_transitions(states)
        
        # Calculate emission probabilities
        self.emission_probs = self._calculate_emissions(sequence, states)
        
        # Calculate initial state probabilities
        cpg_count = sum(1 for s in states if s == 'C')
        total = len(states)
        self.initial_probs = {
            'C': cpg_count / total,
            'N': (total - cpg_count) / total
        }
    
    def _calculate_transitions(self, states):
        """Calculate transition probability matrix"""
        transitions = defaultdict(lambda: defaultdict(int))
        
        for i in range(1, len(states)):
            prev_state = states[i-1]
            curr_state = states[i]
            transitions[prev_state][curr_state] += 1
        
        # Normalize to probabilities
        trans_probs = {}
        for from_state in self.states:
            total = sum(transitions[from_state].values())
            trans_probs[from_state] = {}
            for to_state in self.states:
                count = transitions[from_state][to_state]
                trans_probs[from_state][to_state] = count / total if total > 0 else 0.5
        
        return trans_probs
    
    def _calculate_emissions(self, sequence, states):
        """Calculate emission probability matrix"""
        emissions = defaultdict(lambda: defaultdict(int))
        
        for i, nucleotide in enumerate(sequence):
            if nucleotide in self.nucleotides:
                state = states[i]
                emissions[state][nucleotide] += 1
        
        # Normalize to probabilities
        emit_probs = {}
        for state in self.states:
            total = sum(emissions[state].values())
            emit_probs[state] = {}
            for nt in self.nucleotides:
                count = emissions[state][nt]
                emit_probs[state][nt] = count / total if total > 0 else 0.25
        
        return emit_probs
    
    def viterbi(self, sequence):
        """Viterbi algorithm for most likely state sequence"""
        n = len(sequence)
        
        # Initialize DP table
        dp = np.zeros((len(self.states), n))
        path = np.zeros((len(self.states), n), dtype=int)
        
        state_to_idx = {state: i for i, state in enumerate(self.states)}
        
        # Initialize first column
        for i, state in enumerate(self.states):
            if sequence[0] in self.nucleotides:
                emission_prob = self.emission_probs[state][sequence[0]]
                dp[i, 0] = np.log(self.initial_probs[state] + 1e-10) + np.log(emission_prob + 1e-10)
        
        # Fill DP table
        for t in range(1, n):
            if sequence[t] not in self.nucleotides:
                continue
                
            for curr_idx, curr_state in enumerate(self.states):
                emission_prob = self.emission_probs[curr_state][sequence[t]]
                
                best_prob = -np.inf
                best_prev = 0
                
                for prev_idx, prev_state in enumerate(self.states):
                    trans_prob = self.transition_probs[prev_state][curr_state]
                    prob = dp[prev_idx, t-1] + np.log(trans_prob + 1e-10) + np.log(emission_prob + 1e-10)
                    
                    if prob > best_prob:
                        best_prob = prob
                        best_prev = prev_idx
                
                dp[curr_idx, t] = best_prob
                path[curr_idx, t] = best_prev
        
        # Backtrack to find best path
        states_seq = ['N'] * n
        best_last_state = np.argmax(dp[:, -1])
        states_seq[-1] = self.states[best_last_state]
        
        for t in range(n-2, -1, -1):
            best_last_state = path[best_last_state, t+1]
            states_seq[t] = self.states[best_last_state]
        
        return states_seq

class CpGAnalyzer:
    """Complete CpG island analysis toolkit"""
    
    def __init__(self):
        self.hmm = CpGHMM()
    
    def sliding_window_detection(self, sequence, window_size=500, threshold=0.1):
        """Original sliding window method"""
        cpg_regions = []
        
        for i in range(0, len(sequence), window_size):
            window = sequence[i:i + window_size]
            cg_count = window.count("CG")
            cg_density = cg_count / len(window) if len(window) > 0 else 0
            
            if cg_density >= threshold:
                cpg_regions.append((i, i + window_size))
        
        return cpg_regions
    
    def get_cpg_positions(self, regions, sequence_length):
        """Convert regions to individual positions"""
        positions = []
        for start, end in regions:
            positions.extend(range(start, min(end, sequence_length)))
        return positions
    
    def compare_methods(self, sequence, window_size=500, threshold=0.1):
        """Compare sliding window vs HMM methods"""
        # Sliding window method
        sw_regions = self.sliding_window_detection(sequence, window_size, threshold)
        sw_positions = self.get_cpg_positions(sw_regions, len(sequence))
        
        # Train HMM on sliding window results
        self.hmm.train(sequence, sw_positions)
        
        # HMM prediction
        hmm_states = self.hmm.viterbi(sequence)
        hmm_positions = [i for i, state in enumerate(hmm_states) if state == 'C']
        
        return {
            'sliding_window': {
                'regions': sw_regions,
                'positions': sw_positions,
                'count': len(sw_regions)
            },
            'hmm': {
                'states': hmm_states,
                'positions': hmm_positions,
                'count': len([i for i in hmm_states if i == 'C'])
            }
        }
    
    def statistical_comparison(self, chr_data):
        """Statistical analysis of chromosome differences"""
        results = {}
        
        chromosomes = list(chr_data.keys())
        cpg_counts = [chr_data[chr]['sliding_window']['count'] for chr in chromosomes]
        chr_lengths = [len(chr_data[chr]['sequence']) for chr in chromosomes]
        
        # Normalize by chromosome length
        cpg_densities = [count/length * 1000000 for count, length in zip(cpg_counts, chr_lengths)]
        
        results['counts'] = dict(zip(chromosomes, cpg_counts))
        results['densities'] = dict(zip(chromosomes, cpg_densities))
        
        # Statistical tests
        if len(chromosomes) >= 3:
            # ANOVA for multiple comparisons
            f_stat, p_value = stats.f_oneway(*[chr_data[chr]['sliding_window']['count'] 
                                              for chr in chromosomes])
            results['anova'] = {'f_statistic': f_stat, 'p_value': p_value}
        
        # Pairwise comparisons
        results['pairwise'] = {}
        for i in range(len(chromosomes)):
            for j in range(i+1, len(chromosomes)):
                chr1, chr2 = chromosomes[i], chromosomes[j]
                # Chi-square test for independence
                observed = [[cpg_counts[i], chr_lengths[i] - cpg_counts[i]],
                           [cpg_counts[j], chr_lengths[j] - cpg_counts[j]]]
                chi2, p_val, _, _ = stats.chi2_contingency(observed)
                results['pairwise'][f"{chr1}_vs_{chr2}"] = {
                    'chi2': chi2, 'p_value': p_val
                }
        
        return results
    
    def create_comprehensive_plot(self, chr_data, stats_results):
        """Create comprehensive visualization"""
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        chromosomes = list(chr_data.keys())
        
        # 1. CpG island counts comparison
        counts = [chr_data[chr]['sliding_window']['count'] for chr in chromosomes]
        axes[0, 0].bar(chromosomes, counts, color=['blue', 'red', 'green'])
        axes[0, 0].set_title('CpG Island Counts by Chromosome')
        axes[0, 0].set_ylabel('Number of CpG Islands')
        
        # 2. CpG density comparison (normalized by length)
        densities = [stats_results['densities'][chr] for chr in chromosomes]
        axes[0, 1].bar(chromosomes, densities, color=['blue', 'red', 'green'])
        axes[0, 1].set_title('CpG Island Density (per Mb)')
        axes[0, 1].set_ylabel('CpG Islands per Million bp')
        
        # 3. Method comparison for one chromosome
        sample_chr = chromosomes[0]
        sample_data = chr_data[sample_chr]
        sw_count = sample_data['sliding_window']['count']
        hmm_count = sample_data['hmm']['count']
        
        methods = ['Sliding Window', 'HMM']
        method_counts = [sw_count, hmm_count]
        axes[1, 0].bar(methods, method_counts, color=['orange', 'purple'])
        axes[1, 0].set_title(f'Method Comparison - Chromosome {sample_chr}')
        axes[1, 0].set_ylabel('CpG Positions Detected')
        
        # 4. Statistical significance heatmap
        if 'pairwise' in stats_results:
            comparisons = list(stats_results['pairwise'].keys())
            p_values = [stats_results['pairwise'][comp]['p_value'] for comp in comparisons]
            
            # Create a simple visualization of p-values
            colors = ['red' if p < 0.05 else 'gray' for p in p_values]
            axes[1, 1].bar(range(len(comparisons)), [-np.log10(p) for p in p_values], 
                          color=colors)
            axes[1, 1].set_title('Statistical Significance (-log10 p-value)')
            axes[1, 1].set_ylabel('-log10(p-value)')
            axes[1, 1].set_xticks(range(len(comparisons)))
            axes[1, 1].set_xticklabels(comparisons, rotation=45)
            axes[1, 1].axhline(y=-np.log10(0.05), color='black', linestyle='--', 
                              label='p=0.05 threshold')
            axes[1, 1].legend()
        
        plt.tight_layout()
        return fig
    
    def generate_report(self, chr_data, stats_results):
        """Generate comprehensive analysis report"""
        report = []
        report.append("=== CpG Island Analysis Report ===\n")
        
        # Summary statistics
        report.append("SUMMARY STATISTICS:")
        for chr_name, data in chr_data.items():
            seq_length = len(data['sequence'])
            cpg_count = data['sliding_window']['count']
            density = stats_results['densities'][chr_name]
            
            report.append(f"Chromosome {chr_name}:")
            report.append(f"  - Length: {seq_length:,} bp")
            report.append(f"  - CpG Islands: {cpg_count}")
            report.append(f"  - Density: {density:.2f} islands/Mb")
        
        # Statistical analysis
        report.append("\nSTATISTICAL ANALYSIS:")
        if 'anova' in stats_results:
            anova = stats_results['anova']
            report.append(f"ANOVA F-statistic: {anova['f_statistic']:.4f}")
            report.append(f"ANOVA p-value: {anova['p_value']:.4e}")
            significance = "significant" if anova['p_value'] < 0.05 else "not significant"
            report.append(f"Overall difference is {significance}")
        
        report.append("\nPAIRWISE COMPARISONS:")
        for comparison, results in stats_results['pairwise'].items():
            chi2 = results['chi2']
            p_val = results['p_value']
            significance = "significant" if p_val < 0.05 else "not significant"
            report.append(f"{comparison}: χ² = {chi2:.4f}, p = {p_val:.4e} ({significance})")
        
        # Method comparison
        report.append("\nMETHOD COMPARISON:")
        for chr_name, data in chr_data.items():
            sw_count = data['sliding_window']['count']
            hmm_count = data['hmm']['count']
            agreement = abs(sw_count - hmm_count) / max(sw_count, hmm_count) * 100
            report.append(f"Chromosome {chr_name}: SW={sw_count}, HMM={hmm_count} "
                         f"(difference: {agreement:.1f}%)")
        
        return "\n".join(report)

# Example usage function
def analyze_chromosomes(chr_files):
    """
    Main analysis function
    chr_files should be a dict like {'Z': 'chrZ.fasta', 'W': 'chrW.fasta', '28': 'chr28.fasta'}
    """
    analyzer = CpGAnalyzer()
    chr_data = {}
    
    # Load and analyze each chromosome
    for chr_name, filename in chr_files.items():
        print(f"Analyzing chromosome {chr_name}...")
        
        # Load sequence (assuming FASTA format without headers in content)
        with open(filename, 'r') as f:
            sequence = ''.join(line.strip() for line in f if not line.startswith('>'))
        
        # Analyze with both methods
        results = analyzer.compare_methods(sequence)
        results['sequence'] = sequence
        chr_data[chr_name] = results
    
    # Statistical comparison
    stats_results = analyzer.statistical_comparison(chr_data)
    
    # Generate visualizations
    fig = analyzer.create_comprehensive_plot(chr_data, stats_results)
    
    # Generate report
    report = analyzer.generate_report(chr_data, stats_results)
    
    return chr_data, stats_results, fig, report

# To use this enhanced analysis:
# chr_files = {'Z': 'chrZ.fasta', 'W': 'chrW.fasta', '28': 'chr28.fasta'}
# data, stats, plot, report = analyze_chromosomes(chr_files)
# print(report)
# plot.show()
