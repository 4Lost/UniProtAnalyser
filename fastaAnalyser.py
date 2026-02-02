import glob
import json
from collections import defaultdict
import math
import os
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Queue

N_WORKERS = 12

FILE_NAME = "aa-uniref50.fasta"
#FILE_NAME = "exampleDatasetSmall.fasta"

# /-- aaindex Properties --/
AMINO_ORDER = "ARNDCQEGHILKMFPSTWYV"
FEATURE_SCALES = {
  "hydrophobicity": {"A":0.14358913378319937, "R":-0.22636644169555728, "N":-0.2092242521622241, "D":-0.26473494857761126, "C":0.32931700684958476, "Q":-0.22784290932653717, "E":-0.2809654937064575, "G":-0.05574240859408246, "H":-0.0005694593018983487, "I":0.7394661181561407, "L":0.6727226738424296, "K":-0.21494542772443354, "M":0.4557449156253738, "F":0.76326618347678, "P":0.17041140131272003, "S":-0.1371795756074432, "T":-0.005209995452296512, "W":0.7069106029546537, "Y":0.42808513792254455, "V":0.5120052529094391}, 
  "freeEnergie": {"A":0.4126725685975706, "R":0.3976573097883449, "N":0.5376667922318099, "D":0.5498713549024209, "C":0.42070550702821014, "Q":0.4305812297852952, "E":0.3931500062963445, "G":0.6419613997574174, "H":0.44649404914002194, "I":0.32993634469537864, "L":0.37253441201227444, "K":0.4297431453605854, "M":0.35660015135749124, "F":0.378033855997586, "P":1.0, "S":0.4656483843114838, "T":0.4145832124838414, "W":0.3879077058059451, "Y":0.38068244136095986, "V":0.3340384406240167}, 
  "stability": {"A":0.3480633629471504, "R":0.34399187846000145, "N":0.2073746228973096, "D":0.12041433302943207, "C":0.5407968873115502, "Q":0.28641797103574157, "E":0.27396796651284466, "G":0.038922249361399465, "H":0.33801074196451336, "I":0.8616284467863063, "L":0.7118262537208974, "K":0.38132182563329614, "M":0.5233023271400615, "F":0.8570013100720518, "P":0.34592244625918867, "S":0.28260349645510846, "T":0.4461459771266177, "W":0.8491195216807631, "Y":0.7288116232206955, "V":0.7828956685436782}, 
  "volume": {"A":0.36602623291450886, "R":0.8238919538564463, "N":0.5371914111807752, "D":0.4837632806855753, "C":0.48213006245979795, "Q":0.6464190909849685, "E":0.617329380563794, "G":0.231896382746156, "H":0.6882483635889526, "I":0.6993281023998736, "L":0.6977921936967111, "K":0.7286885100866508, "M":0.7189429802894034, "F":0.8328210272765146, "P":0.47900880512587585, "S":0.391636744280506, "T":0.5098512092286421, "W":0.9792804979957268, "Y":0.8766364342249389, "V":0.5940467373991879}, 
  "alphaHelix": {"A":0.7442034901304289, "R":0.49965192974183026, "N":0.2922649694733262, "D":0.39245804612261465, "C":0.3041156248788227, "Q":0.5808341080712011, "E":0.7501169192456335, "G":0.07290246184485562, "H":0.4961365750670697, "I":0.4896085869580447, "L":0.7216557552587273, "K":0.5799675060484489, "M":0.794644214563068, "F":0.5549011221246851, "P":-0.03535290560938736, "S":0.25080549451114953, "T":0.2969457164912867, "W":0.5408691741732791, "Y":0.28905098834816, "V":0.4830710076268101}, 
  "betaSheet": {"A":0.34956591911092894, "R":0.460550677798243, "N":0.22256084553921518, "D":0.18452618634543372, "C":0.5099790895042989, "Q":0.4453882634893641, "E":0.11317872915496889, "G":0.3274437242442998, "H":0.3872233691659867, "I":0.9199532659318971, "L":0.5989115344257007, "K":0.27558096039030416, "M":0.6135421560885467, "F":0.7066307354658885, "P":0.09494029912392414, "S":0.3628039141726376, "T":0.5799367947829297, "W":0.6295171062240904, "Y":0.7336981819919552, "V":0.8903504557212697}, 
  "coil": {"A":0.09663551689300753, "R":0.2678770347284752, "N":0.6591846813607345, "D":0.4775261018204119, "C":0.27029350887007075, "Q":0.1909684984175682, "E":0.15314799173070323, "G":0.6665124673407554, "H":0.36393670117956045, "I":-0.08149077849430328, "L":-0.05275135589796503, "K":0.3039878566804798, "M":-0.1346607990280012, "F":0.0003032555159309797, "P":0.9318164083836635, "S":0.5001012727728252, "T":0.3513196863356348, "W":0.07066763248303932, "Y":0.18765062236002322, "V":-0.13569056048760558}, 
  "mutability": {"A":0.8394984282043829, "R":0.7005642770375807, "N":0.935471554028255, "D":0.8076649092189944, "C":0.46720358906849135, "Q":0.7524464644790849, "E":0.7617643564712345, "G":0.580491892462965, "H":0.7086617314344656, "I":0.7934216484762722, "L":0.5042077479375525, "K":0.6777650210486031, "M":0.7530836433744573, "F":0.42598335312599045, "P":0.5863325894099859, "S":0.9617377032363953, "T":0.8622880399307542, "W":0.18472971296482607, "Y":0.4128250719264721, "V":0.7553781011236005}
}


# /-- Import --/
def streamFasta(filename):
    with open(filename, "r") as file:
        seq_id = None
        seq_parts = []

        for line in file:
            line = line.strip()

            if not line:
                continue
            if line.startswith(">"):
                if seq_id:
                    yield seq_id, "".join(seq_parts)
                seq_id = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if seq_id:
            yield seq_id, "".join(seq_parts)

def producer(queue):
    for seqId, seq in streamFasta(FILE_NAME):
        queue.put((seqId, seq))

    # stop signals
    for _ in range(N_WORKERS):
        queue.put(None)


# /-- Analyser --/
def analyseWorker(workerId, queue):
    print(f"Worker {workerId + 1} started")
    perSequenceResults = {feature: {} for feature in FEATURE_SCALES}
    aaResults = (defaultdict(int), defaultdict(lambda: defaultdict(int)), defaultdict(int)) # (distribution, positional distribution, length distribution)
    count = 0
    
    while True:
        item = queue.get()

        if item is None:
            break

        _, seq = item

        for feature in FEATURE_SCALES:
            value = calculateValueForFeature(seq.upper(), FEATURE_SCALES[feature])
            if value is None:
                continue
            perSequenceResults[feature][seq] = value
        
        aaResults = calculateAAValues(seq, aaResults)
        count += 1
        
        if count % 2000 == 0: 
            print(f"Worker {workerId + 1}: {count}")
            saveWorkerFeatureValues(workerId, perSequenceResults)
            perSequenceResults = {feature: {} for feature in FEATURE_SCALES}
            saveWorkerAAValues(workerId, aaResults)
            aaResults = (defaultdict(int), defaultdict(lambda: defaultdict(int)), defaultdict(int))

    saveWorkerFeatureValues(workerId, perSequenceResults)
    saveWorkerAAValues(workerId, aaResults)
    print(f"Worker {workerId + 1} finished")


# /-- Analyser -- Helper --/
def saveWorkerFeatureValues(workerId, perSequenceResults):
    for feature in FEATURE_SCALES:
        with open(f"worker_{workerId}_{feature}.jsonl", "a") as f:
            for seq, value in perSequenceResults[feature].items():
                json.dump({seq: value}, f)
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
    
    n = len(values)
    min_val = np.min(arr)
    max_val = np.max(arr)
    stdDev = np.std(arr)

    bandwidth = 0.75 * stdDev * n**(-1/5)

    kdePoints = []
    valRange = max_val - min_val
    
    for i in range(200):
        x = min_val + (i / 200) * valRange
        y = 0
        for value in values:
            y += math.exp(-((x - value) ** 2) / (2 * bandwidth ** 2))
        
        y /= n * bandwidth * math.sqrt(2 * math.pi)
        kdePoints.append((float(x), float(y)))
    
    sumKDE = sum(point[1] for point in kdePoints) * (valRange / 200)
    normalizedKDE = [(point[0], point[1] / sumKDE) for point in kdePoints]

    mean = np.mean(arr)

    jsonData = {
        "stats": {
            "min": float(min_val),
            "max": float(max_val),
            "mean": float(mean),
            "stdDev": float(stdDev),
            "median": float(np.median(arr)),
        },
        "values": normalizedKDE
    }

    print(f"Saving {feature}")
    saveFeatureValues(jsonData, feature)
    print(f"Saved {feature}")

    return mean, stdDev


def loadFeatureValues(feature):
    data = {}
    for workerId in range(N_WORKERS):
        with open(f"worker_{workerId}_{feature}.jsonl", "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    d = json.loads(line)
                    data.update(d)
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

    mean, stdDev, stats = getLengthStats(lengthDistribution)

    jsonData = {
        "lengthKde": calculateLengthKDE(lengthDistribution),
        "lengthStats": stats,
        "distribution": normalizeDistribution(distribution),
        "positionalDistribution": positionalDistribution,
    }

    print(f"Saving AA")
    with open(f"distribution_AA.json", "w") as fout:
        json.dump(jsonData, fout, indent=2)
    print(f"Saved AA")

    return mean, stdDev


def loadAAValues():
    dist = defaultdict(int)
    posDist = defaultdict(lambda: defaultdict(int))
    lenDist = defaultdict(int)

    for workerId in range(N_WORKERS):
        with open(f"worker_{workerId}_distribution.jsonl", "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    data = json.loads(line)
                    dist = assignDist(data, dist)
                    
    for workerId in range(N_WORKERS):
        with open(f"worker_{workerId}_positionalDistribution.jsonl", "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    data = json.loads(line)
                    posDist = assignPosDist(data, posDist)
                    
    for workerId in range(N_WORKERS):
        with open(f"worker_{workerId}_lengthDistribution.jsonl", "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    data = json.loads(line)
                    lenDist = assignLenDist(data, lenDist)
                    
    return (dist, posDist, lenDist)


def calculateLengthKDE(lenDist):
    lenDist = dict(sorted(lenDist.items()))
    amount = sum(lenDist.values())
    x = 0.05

    max_key = max((k for k, v in lenDist.items() if v / amount > x), default=None)

    if max_key == None:
        max_key = max(lenDist.keys())

    cutoff = {}
    for k in lenDist:
        if float(k) <= float(max_key):
            cutoff[k] = lenDist[k]

    values = []

    for key, value in cutoff.items():
        for i in range(value):
            values.append(float(key))

    valMin = values[0]
    valMax = values[len(values) - 1]
    valRange = valMax - valMin
    bandwidth = 20

    kdePoints = []
    for i in range(200):
        x = valMin + (i / 200) * valRange
        y = 0
        for value in values:
            y += math.exp(-math.pow(x - value, 2) / (2 * bandwidth * bandwidth))
      
        y /= len(values) * bandwidth * math.sqrt(2 * math.pi)
        kdePoints.append((x, y))
    
    sumKDE = sum(point[1] for point in kdePoints) * (valRange / 200)
    normalizedKDE = [(point[0], point[1] / sumKDE) for point in kdePoints]
    return normalizedKDE


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

    stdDev = np.sqrt(variance)

    return mean, stdDev, {
        "min": int(np.min(lengths)),
        "max": int(np.max(lengths)),
        "mean": float(mean),
        "stdDev": float(stdDev),
    }


def normalizeDistribution(dist):
    total = sum(dist.values())
    if total == 0:
        return {k: 0.0 for k in dist}
    return {k: (v / total) * 100 for k, v in dist.items()}


def cleanup_worker_files():
    for file in glob.glob("worker_*_*.jsonl"):
        os.remove(file)
        print(f"Removed {file}")


# /-- Surprise Metric --/
def featureProducer(feature, queue):
    for seqId, seq in streamFeature(feature):
        queue.put((seqId, seq))

    # stop signals
    for _ in range(N_WORKERS):
        queue.put(None)


def streamFeature(feature):
    for workerId in range(N_WORKERS):
        with open(f"worker_{workerId}_{feature}.jsonl", "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    d = json.loads(line)
                    seq_id, value = next(iter(d.items()))
                    yield seq_id, value


def featureProducer(feature, queue):
    for seqId, seq in streamFeature(feature):
        queue.put((seqId, seq))

    # stop signals
    queue.put(None)


def surpriseMetric(queues, means, stdDevs):
    run = True
    seqSave = ''
    perSequenceResults = {}
    classStats = {'': 0, 'extremly surprising': 0, 'highly surprising': 0, 'surprising': 0, 'slightly surprising': 0, 'ordinary': 0}
    count = 0

    while run:
        values = []
        for i, queue in enumerate(queues):
            item = queue.get()

            if item is None:
                run = False
                break

            seq, value = item
            
            if i == 0:
                seqSave = seq

            if seqSave != seq:
                print('Error: ${seqSqve} != ${seq}')

            values.append(value)

        if run == False:
            break
        
        values.append(len(seq))
        count += 1

        surpriseClass, classification = classifySequence(values, means, stdDevs)
        perSequenceResults[seq] = classification

        classStats[surpriseClass] += 1

        if count % 2000 == 0: 
            print(f"SurpriseMetric: {count}")
            saveSurpriseMetric(perSequenceResults)
            perSequenceResults = {}

    saveSurpriseMetric(perSequenceResults)
    saveMetricStats(classStats)


def saveSurpriseMetric(perSequenceResults):
    with open(f"surpriseMetric.jsonl", "a") as f:
        for seq, value in perSequenceResults.items():
            json.dump({seq: value}, f)
            f.write("\n")


def saveMetricStats(surpriseMetricStats):
    with open(f"surpriseMetric.json", "w") as f:
        json.dump(surpriseMetricStats, f, indent=2)


def classifySequence(values, means, stdDevs):
    z_scores = {'class': '', 'surpriseFactor': 0}
    totalFactor = 0

    if len(values) < len(means):
        return '', z_scores

    for i, feature in enumerate(FEATURE_SCALES):
        val = abs(means[i] - values[i]) / stdDevs[i]
        z_scores[feature] = val
        totalFactor += val
    
    val = abs(means[8] - values[8]) / stdDevs[8]
    z_scores['length'] = val
    totalFactor += val

    totalFactor = totalFactor / 9
    z_scores['class'] = totalFactor

    if totalFactor >= 5.0:
        z_scores['class'] = 'extremly surprising'
    elif totalFactor >= 4:
        z_scores['class'] = 'highly surprising'
    elif totalFactor >= 3:
        z_scores['class'] = 'surprising'
    elif totalFactor >= 2.5:
        z_scores['class'] = 'slightly surprising'
    else:
        z_scores['class'] = 'ordinary'

    return z_scores['class'], z_scores


# /-- Main --/
def main():
    from multiprocessing import Process, Queue
    print("Start streaming analysis")
    queue = Queue(maxsize=5000)  # important: backpressure

    # start workers
    workers = []
    for workerId in range(N_WORKERS):
        p = Process(target=analyseWorker, args=(workerId, queue))
        p.start()
        workers.append(p)

    # run producer in main process
    producer(queue)

    # wait for workers
    for p in workers:
        p.join()

    means = []
    stdDevs = []

    print("\nGenerate Scales & AA Stats and Save Results\n")
    for feature in FEATURE_SCALES:
        mean, stdDev = processFeature(feature)
        means.append(mean)
        stdDevs.append(stdDev)
    mean, stdDev = processAA()
    means.append(mean)
    stdDevs.append(stdDev)

    print("\nGenerate Surprise Metric\n")

    queues = [Queue(maxsize=2000), Queue(maxsize=2000), Queue(maxsize=2000), Queue(maxsize=2000), Queue(maxsize=2000), Queue(maxsize=2000), Queue(maxsize=2000), Queue(maxsize=2000)]
    workers = []
    for number, feature in enumerate(FEATURE_SCALES):
        p = Process(target=featureProducer, args=(feature, queues[number]))
        p.start()
        workers.append(p)
    workers.append(p)

    # run producer in main process
    surpriseMetric(queues, means, stdDevs)

    # wait for workers
    for p in workers:
        p.join()


    print(f"Done! All Results written!")
    cleanup_worker_files()


if __name__ == "__main__":
    main()

