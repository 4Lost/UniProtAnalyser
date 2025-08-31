import numpy as np
from collections import defaultdict
from scipy.stats import gaussian_kde
import json

name = "uniref50.fasta"

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

def generate_kde_points(length_distribution):
    num_points=100
    length_list = []
    points = []

    for length, count in length_distribution.items():
        length_list.extend([length] * count)

    length_array = np.array(length_list)

    kde = gaussian_kde(length_array)

    x_points = np.linspace(length_array.min(), length_array.max(), num_points)
    y_points = kde(x_points)

    y_total = y_points.sum()
    y_points_normalized = y_points / y_points.max()

    for x, y in zip(x_points, y_points_normalized):
        points.append((x, y))

    # Stats calculation
    lengths = np.array(list(length_distribution.keys()))
    counts = np.array(list(length_distribution.values()), dtype=float)

    mean = np.average(lengths, weights=counts)
    variance = np.average((lengths - mean)**2, weights=counts)
    stats = {
        "min": int(lengths.min()),
        "max": int(lengths.max()),
        "mean": float(mean),
        "variance": float(variance),
        "std_dev": float(np.sqrt(variance))
    }

    return points, stats

print("Read Sequences")
sequences = readFasta(name)

print("")
print("Amount of sequences: " + str(len(sequences)))
print("")

print("Calculate Distribution")

protDist = proteinDistribution(sequences)
length_kde, length_stats = generate_kde_points(protDist[0])

data_to_save = {
    "num_sequences": len(sequences),
    "length_kde": length_kde,
    "length_stats": length_stats,
    "distribution": dict(protDist[1]),
    "positional_distribution": {str(pos+1): dict(counts) for pos, counts in protDist[2].items()}
}

output_file = "distribution_AA.json"
with open(output_file, "w") as f:
    json.dump(data_to_save, f, indent=2)

print(f"Done! Results written to {output_file}")