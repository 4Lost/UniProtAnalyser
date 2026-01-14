import json
from collections import defaultdict
import math
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

N_WORKERS = 12

FILE_NAME = "aa-uniref50.fasta"
#FILE_NAME = "exampleDatasetSmall.fasta"

# /-- aaindex Properties --/
AMINO_ORDER = "ARNDCQEGHILKMFPSTWYV"
FEATURE_SCALES = {
  "hydroscales": {"A":0.14358913378319937, "R":-0.22636644169555728, "N":-0.2092242521622241, "D":-0.26473494857761126, "C":0.32931700684958476, "Q":-0.22784290932653717, "E":-0.2809654937064575, "G":-0.05574240859408246, "H":-0.0005694593018983487, "I":0.7394661181561407, "L":0.6727226738424296, "K":-0.21494542772443354, "M":0.4557449156253738, "F":0.76326618347678, "P":0.17041140131272003, "S":-0.1371795756074432, "T":-0.005209995452296512, "W":0.7069106029546537, "Y":0.42808513792254455, "V":0.5120052529094391}, 
  "freeEnergie": {"A":0.4126725685975706, "R":0.3976573097883449, "N":0.5376667922318099, "D":0.5498713549024209, "C":0.42070550702821014, "Q":0.4305812297852952, "E":0.3931500062963445, "G":0.6419613997574174, "H":0.44649404914002194, "I":0.32993634469537864, "L":0.37253441201227444, "K":0.4297431453605854, "M":0.35660015135749124, "F":0.378033855997586, "P":1.0, "S":0.4656483843114838, "T":0.4145832124838414, "W":0.3879077058059451, "Y":0.38068244136095986, "V":0.3340384406240167}, 
  "stability": {"A":0.3480633629471504, "R":0.34399187846000145, "N":0.2073746228973096, "D":0.12041433302943207, "C":0.5407968873115502, "Q":0.28641797103574157, "E":0.27396796651284466, "G":0.038922249361399465, "H":0.33801074196451336, "I":0.8616284467863063, "L":0.7118262537208974, "K":0.38132182563329614, "M":0.5233023271400615, "F":0.8570013100720518, "P":0.34592244625918867, "S":0.28260349645510846, "T":0.4461459771266177, "W":0.8491195216807631, "Y":0.7288116232206955, "V":0.7828956685436782}, 
  "volume": {"A":0.36602623291450886, "R":0.8238919538564463, "N":0.5371914111807752, "D":0.4837632806855753, "C":0.48213006245979795, "Q":0.6464190909849685, "E":0.617329380563794, "G":0.231896382746156, "H":0.6882483635889526, "I":0.6993281023998736, "L":0.6977921936967111, "K":0.7286885100866508, "M":0.7189429802894034, "F":0.8328210272765146, "P":0.47900880512587585, "S":0.391636744280506, "T":0.5098512092286421, "W":0.9792804979957268, "Y":0.8766364342249389, "V":0.5940467373991879}, 
  "alpha-helix": {"A":0.7442034901304289, "R":0.49965192974183026, "N":0.2922649694733262, "D":0.39245804612261465, "C":0.3041156248788227, "Q":0.5808341080712011, "E":0.7501169192456335, "G":0.07290246184485562, "H":0.4961365750670697, "I":0.4896085869580447, "L":0.7216557552587273, "K":0.5799675060484489, "M":0.794644214563068, "F":0.5549011221246851, "P":-0.03535290560938736, "S":0.25080549451114953, "T":0.2969457164912867, "W":0.5408691741732791, "Y":0.28905098834816, "V":0.4830710076268101}, 
  "beta-sheet": {"A":0.34956591911092894, "R":0.460550677798243, "N":0.22256084553921518, "D":0.18452618634543372, "C":0.5099790895042989, "Q":0.4453882634893641, "E":0.11317872915496889, "G":0.3274437242442998, "H":0.3872233691659867, "I":0.9199532659318971, "L":0.5989115344257007, "K":0.27558096039030416, "M":0.6135421560885467, "F":0.7066307354658885, "P":0.09494029912392414, "S":0.3628039141726376, "T":0.5799367947829297, "W":0.6295171062240904, "Y":0.7336981819919552, "V":0.8903504557212697}, 
  "coil": {"A":0.09663551689300753, "R":0.2678770347284752, "N":0.6591846813607345, "D":0.4775261018204119, "C":0.27029350887007075, "Q":0.1909684984175682, "E":0.15314799173070323, "G":0.6665124673407554, "H":0.36393670117956045, "I":-0.08149077849430328, "L":-0.05275135589796503, "K":0.3039878566804798, "M":-0.1346607990280012, "F":0.0003032555159309797, "P":0.9318164083836635, "S":0.5001012727728252, "T":0.3513196863356348, "W":0.07066763248303932, "Y":0.18765062236002322, "V":-0.13569056048760558}, 
  "mutability": {"A":0.8394984282043829, "R":0.7005642770375807, "N":0.935471554028255, "D":0.8076649092189944, "C":0.46720358906849135, "Q":0.7524464644790849, "E":0.7617643564712345, "G":0.580491892462965, "H":0.7086617314344656, "I":0.7934216484762722, "L":0.5042077479375525, "K":0.6777650210486031, "M":0.7530836433744573, "F":0.42598335312599045, "P":0.5863325894099859, "S":0.9617377032363953, "T":0.8622880399307542, "W":0.18472971296482607, "Y":0.4128250719264721, "V":0.7553781011236005}
}


# /-- Import --/
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


# /-- Analyser --/
def analyseChunk(workerId):
    print(f"Worker {workerId + 1} started")
    count = 0
    sequences = loadForWorker(workerId).items()

    print(f"Worker {workerId + 1}: Sequences loaded")
    perSequenceResults = {feature: {} for feature in FEATURE_SCALES}
    aaResults = (defaultdict(int), defaultdict(lambda: defaultdict(int)), defaultdict(int)) # (distribution, positional distribution, length distribution)

    for seqId, seq in sequences:
        for feature in FEATURE_SCALES:
            value = calculateValueForFeature(seq.upper(), FEATURE_SCALES[feature])
            if value is None:
                continue
            perSequenceResults[feature][seqId] = value
        
        aaResults = calculateAAValues(seq, aaResults)
        count += 1
        
        if count % 2000 == 0: 
            print(f"Worker {workerId + 1}: {count}/{len(sequences)}")
            saveWorkerFeatureValues(workerId, perSequenceResults)
            perSequenceResults = {feature: {} for feature in FEATURE_SCALES}
            saveWorkerAAValues(workerId, aaResults)
            aaResults = (defaultdict(int), defaultdict(lambda: defaultdict(int)), defaultdict(int))

    saveWorkerFeatureValues(workerId, perSequenceResults)
    saveWorkerAAValues(workerId, aaResults)
    print(f"Worker {workerId + 1} finished")


# /-- Analyser -- Helper --/
def initWorker(sequences):
    items = list(sequences.items())
    chunk_size = math.ceil(len(items) / N_WORKERS)
    sequenceChunks = [items[i:i+chunk_size] for i in range(0, len(items), chunk_size)]

    for workerId in range(N_WORKERS):
        sequenceChunkDict = {}
        for seqId, seq in sequenceChunks[workerId]:
            sequenceChunkDict[seqId] = seq
        with open(f"worker_{workerId}.json", "w") as fout:
            json.dump(sequenceChunkDict, fout, indent=2)
        for feature in FEATURE_SCALES:
            with open(f"worker_{workerId}_{feature}.jsonl", "w") as f:
                pass


def loadForWorker(workerId):
    with open(f"worker_{workerId}.json", "r") as f:
        sequenceChunkDict = json.load(f)
    return sequenceChunkDict


def saveWorkerFeatureValues(workerId, perSequenceResults):
    for feature in FEATURE_SCALES:
        with open(f"worker_{workerId}_{feature}.jsonl", "a") as f:
            for seqId, value in perSequenceResults[feature].items():
                json.dump({seqId: value}, f)
                f.write("\n")


def saveWorkerAAValues(workerId, aaResults):
    (distribution, positionalDistribution, lengthDistribution) = aaResults
    
    with open(f"worker_{workerId}_distribution.jsonl", "a") as f:
        f.write(json.dumps(distribution) + "\n")
    
    with open(f"worker_{workerId}_positionalDistribution.jsonl", "a") as f:
        f.write(json.dumps(positionalDistribution) + "\n")
    
    with open(f"worker_{workerId}_lengthDistribution.jsonl", "a") as f:
        f.write(json.dumps(lengthDistribution) + "\n")


# /-- Analyser -- Scales --/
def calculateValueForFeature(sequence, dict):
    values = [dict[a] for a in sequence if a in dict]
    
    if not values:
        return None
    
    return sum(values) / len(values)


# /-- Analyser -- Amino Acid --/
def calculateAAValues(sequence, aaResults):
    (distribution, positionalDistribution, lengthDistribution) = aaResults
    length = len(sequence)

    for position in range(length):
        letter = sequence[position]

        distribution[letter] += 1
        positionalDistribution[position][letter] += 1
    lengthDistribution[length] += 1

    return (distribution, positionalDistribution, lengthDistribution)


# /-- Stats --/
# /-- Stats -- Feature --/
def processFeature(feature):
    print(f"Worker for {feature} started")
    data = loadFeatureValues(feature)
    
    print(f"{feature} loaded")
    values = list(data.values())
    arr = np.array(values, dtype=float)

    jsonData = {
        "stats": {
            "min": float(np.min(arr)),
            "max": float(np.max(arr)),
            "mean": float(np.mean(arr)),
            "std_dev": float(np.std(arr)),
            "median": float(np.median(arr)),
        },
        "perSequence": data,
    }

    print(f"Saving {feature}")
    saveFeatureValues(jsonData, feature)
    print(f"Saved {feature}")


def loadFeatureValues(feature):
    data = {}
    for workerId in range(N_WORKERS):
        with open(f"worker_{workerId}_{feature}.jsonl", "r") as f:
            for lineNumber, line in enumerate(f):
                line = line.strip()
                if line:
                    try:
                        d = json.loads(line)
                        data.update(d)
                    except json.JSONDecodeError as e:
                        print(f"Error in worker {workerId + 1}'s {feature} at line {lineNumber}: {line}")
                        raise e
    return data


def saveFeatureValues(jsonData, feature):
    with open(f"distribution_scales_{feature}.json", "w") as fout:
        json.dump(jsonData, fout, indent=2)


# /-- Stats -- Feature --/
def processAA():
    print(f"Worker for AA started")
    (distribution, positionalDistribution, lengthDistribution) = loadAAValues()

    for pos in positionalDistribution:
        positionalDistribution[pos] = normalizeDistribution(positionalDistribution[pos])

    print(f"AA loaded")

    jsonData = {
        "lengthKde": lengthDistribution,
        "lengthStats": getLengthStats(lengthDistribution),
        "distribution": normalizeDistribution(distribution),
        "positionalDistribution": positionalDistribution,
    }

    print(f"Saving AA")
    with open(f"distribution_AA.json", "w") as fout:
        json.dump(jsonData, fout, indent=2)
    print(f"Saved AA")


def loadAAValues():
    dist = defaultdict(int)
    posDist = defaultdict(lambda: defaultdict(int))
    lenDist = defaultdict(int)

    for workerId in range(N_WORKERS):
        with open(f"worker_{workerId}_distribution.jsonl", "r") as f:
            for lineNumber, line in enumerate(f):
                line = line.strip()
                if line:
                    try:
                        data = json.loads(line)
                        dist = assignDist(data, dist)
                    except json.JSONDecodeError as e:
                        print(f"Error in worker {workerId + 1}'s dirstibution at line {lineNumber}: {line}")
                        raise e
                    
    for workerId in range(N_WORKERS):
        with open(f"worker_{workerId}_positionalDistribution.jsonl", "r") as f:
            for lineNumber, line in enumerate(f):
                line = line.strip()
                if line:
                    try:
                        data = json.loads(line)
                        posDist = assignPosDist(data, posDist)
                    except json.JSONDecodeError as e:
                        print(f"Error in worker {workerId + 1}'s positional dirstibution at line {lineNumber}: {line}")
                        raise e
                    
    for workerId in range(N_WORKERS):
        with open(f"worker_{workerId}_lengthDistribution.jsonl", "r") as f:
            for lineNumber, line in enumerate(f):
                line = line.strip()
                if line:
                    try:
                        data = json.loads(line)
                        lenDist = assignLenDist(data, lenDist)
                    except json.JSONDecodeError as e:
                        print(f"Error in worker {workerId + 1}'s length dirstibution at line {lineNumber}: {line}")
                        raise e
                    
    return (dist, posDist, lenDist)


def assignDist(data, dist):
    for letter in data:
        dist[letter] += data[letter]
    return dist


def assignPosDist(data, posDist):
    for pos in data:
        for letter in data[pos]:
            posDist[pos][letter] += data[pos][letter]
    return posDist


def assignLenDist(data, lenDist):
    for length in data:
        lenDist[length] += data[length]
    return lenDist


def getLengthStats(length_distribution):
    lengths = np.array(
        [int(k) for k in length_distribution.keys()],
        dtype=np.float64
    )
    counts = np.array(list(length_distribution.values()), dtype=float)
    mean = np.average(lengths, weights=counts)
    variance = np.average((lengths - mean)**2, weights=counts)

    return {
        "min": int(np.min(lengths)),
        "max": int(np.max(lengths)),
        "mean": float(mean),
        "std_dev": float(np.sqrt(variance)),
    }


def normalizeDistribution(dist):
    total = sum(dist.values())
    if total == 0:
        return {k: 0.0 for k in dist}
    return {k: (v / total) * 100 for k, v in dist.items()}


# /-- Main --/
def main():
    print("Read Sequences")
    sequences = readFasta(FILE_NAME)

    print("\nCalculate Scales & AA Distribution\n")
    initWorker(sequences)
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futures = [ex.submit(analyseChunk, workerId) for workerId in range(N_WORKERS)]
        for future in as_completed(futures):
            future.result()

    print("\nGenerate Scales & AA Stats and Save Results\n")
    with ProcessPoolExecutor(max_workers=9) as ex:
        futures = [ex.submit(processFeature, feature) for feature in FEATURE_SCALES] + [ex.submit(processAA)]

        for future in futures:
            future.result()

    print(f"Done! All Results written!")


if __name__ == "__main__":
    main()

