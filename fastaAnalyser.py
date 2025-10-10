import numpy as np
from collections import defaultdict
from scipy.stats import gaussian_kde
import json

# name = "uniref50.fasta"
name = "cluster_reps.fasta"

def readFasta(filename):
    sequences = {}
    with open(filename, "r") as file:
        seq_id = None
        seq_list = []

        for line in file:
            line = line.strip()

            if not line:
                continue
            if line.startswith(">"):
                if seq_id:
                    sequences[seq_id] = "".join(seq_list)
                seq_id = line[1:].split()[0]
                seq_list = []
            else:
                seq_list.append(line)
        if seq_id:
            sequences[seq_id] = "".join(seq_list)
    return sequences

def proteinDistribution(sequences):
    distribution = defaultdict(int)
    positionalDistribution = defaultdict(lambda: defaultdict(int))
    lengthDistribution = defaultdict(int)

    for key in sequences:
        length = len(sequences[key])

        for position in range(length):
            letter = sequences[key][position]

            distribution[letter] += 1
            positionalDistribution[position][letter] += 1
        lengthDistribution[length] += 1
    
    return [lengthDistribution, distribution, positionalDistribution]

def normalize_distribution(dist):
    """Convert counts into percentages (0–100)."""
    total = sum(dist.values())
    if total == 0:
        return {k: 0.0 for k in dist}
    return {k: (v / total) * 100 for k, v in dist.items()}

def generate_kde_points(length_distribution, num_points=100, crop=True, threshold=0.01):
    lengths = np.array(list(length_distribution.keys()))
    counts = np.array(list(length_distribution.values()), dtype=float)
    
    # Stats calculation
    mean = np.average(lengths, weights=counts)
    variance = np.average((lengths - mean)**2, weights=counts)
    total_counts = counts.sum()

    stats = {
        "min": int(lengths.min()),
        "max": int(lengths.max()),
        "mean": float(mean),
        "variance": float(variance),
        "std_dev": float(np.sqrt(variance)),
    }

    # Always initialize properly
    length_list = []
    for length, count in length_distribution.items():
        fraction = count / total_counts
        if not crop or fraction >= threshold:
            length_list.extend([length] * int(count))

    # Fallback if crop removes everything
    if not length_list:
        for length, count in length_distribution.items():
            length_list.extend([length] * int(count))
    
    length_array = np.array(length_list)

    # Now gaussian_kde will always have enough samples
    kde = gaussian_kde(length_array)
    x_points = np.linspace(length_array.min(), length_array.max(), num_points * 5)
    y_points = kde(x_points)
    y_points_normalized = y_points / y_points.max()

    points = [(x, y) for x, y in zip(x_points, y_points_normalized)]
    return points, stats

print("Read Sequences")
sequences = readFasta(name)

print("")
print("Amount of sequences: " + str(len(sequences)))
print("")

print("Calculate Distribution")

protDist = proteinDistribution(sequences)
length_kde, length_stats = generate_kde_points(protDist[0])

# Convert all distributions into percentages
normalized_distribution = normalize_distribution(protDist[1])
normalized_positional = {
    str(pos+1): normalize_distribution(counts)
    for pos, counts in protDist[2].items()
}

data_to_save = {
    "num_sequences": len(sequences),
    "length_kde": length_kde,
    "length_stats": length_stats,
    "distribution": normalized_distribution,  # already %
    "positional_distribution": normalized_positional,  # already %
}

output_file = "distribution_AA.json"
with open(output_file, "w") as f:
    json.dump(data_to_save, f, indent=2)

print(f"Done! Results written to {output_file}")
