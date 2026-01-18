import json
from aaindex import aaindex1

# /-- aaindex Properties --/
AMINO_ORDER = "ARNDCQEGHILKMFPSTWYV"
HYDROPHOBICITY_SCALES = [
        "MIYS850101", "ROBB790101", "ARGP820101", "CIDH920101", "CIDH920102",
        "CIDH920103", "CIDH920104", "CIDH920105", "EISD840101", "EISD860103",
        "FAUJ830101", "GOLD730101", "JOND750101", "JURD980101", "KYTJ820101",
        "NADH010102", "NADH010103", "NADH010104", "PONP930101", "RADA880108",
        "SWER830101", "WILM950101", "WILM950102", "ZIMJ680101", "BROC820101",
        "GUOD860101", "MEEJ800101", "MEEJ800102", "MEEJ810101", "MEEJ810102",
        "BROC820102", "EISD860101", "NOZY710101", "LAWE840101", "RADA880101",
        "RADA880102", "SIMZ760101", "WIMW960101",
    ]
FREE_ENERGIE_SCALES = [
        "MUNV940101", "MUNV940102", "MUNV940103", "MUNV940104", "MUNV940105",
    ]
STABILITY_SCALES = [
        "BASU050102", "TAKK010101", "ZHOH040101", "ZHOH040102", "OOBM850101",
        "PTIO830102", "FAUJ880107",
    ]
VOLUME_SCALES = [ # "LINS030101", 
        "CHOC760101", "RADA880106", "ROSG850101", "BIGC670101",
        "CHOC750101", "DAWD720101", "FAUJ880103", "GOLD730102", "GRAR740103",
        "HARY940101", "KRIW790103", "LEVM760106", "PONJ960101", "TSAJ990101",
        "TSAJ990102", "FASG760101",
    ]
ALPHA_HELIX_SCALES = [
        "AURR980109", "AURR980111", "AURR980112", "AURR980113", "BEGF750101",
        "BURA740101", "CHOP780201", "GEIM800101", "GEIM800104", "ISOY800101",
        "KANM800101", "KANM800103", "LEVM780101", "LEVM780104", "MAXF760101",
        "NAGK730101", "PALJ810101", "PALJ810102", "PALJ810109", "PRAM900102",
        "QIAN880105", "QIAN880109", "QIAN880110", "RACS820108", "RICJ880112",
        "ROBB760101", "TANS770101", "AURR980110", "CRAJ730101", "QIAN880106",
        "QIAN880107", "QIAN880108", "RICJ880109", "ROBB760103", "AURR980108",
        "QIAN880104",
    ]
BETA_SHEET_SCALES = [
        "RACS820105", "BEGF750102", "CHOP780202", "CHOP780208", "CHOP780209",
        "CRAJ730102", "KANM800102", "LEVM780102", "LEVM780105", "NAGK730102",
        "PALJ810103", "PALJ810104", "PALJ810110", "PALJ810111", "PALJ810112",
        "PRAM900103", "QIAN880118", "QIAN880119", "QIAN880120", "QIAN880121",
        "ROBB760106",
    ]
COIL_SCALES = [
        "GEIM800108", "GEIM800111", "CHAM830101", "GEIM800110", "NAGK730103",
        "QIAN880130", "QIAN880131", "QIAN880132", "QIAN880133", "QIAN880134",
        "QIAN880135", "ROBB760112", "CHOP780207",
    ]
MUTABILITY_SCALES = [
        "DAYM780201", "JOND920102", "WEBA780101",
    ]
SCALE_IDS = {
        "hydroscales": HYDROPHOBICITY_SCALES,
        "freeEnergie": FREE_ENERGIE_SCALES,
        "stability": STABILITY_SCALES,
        "volume": VOLUME_SCALES,
        "alpha-helix": ALPHA_HELIX_SCALES,
        "beta-sheet": BETA_SHEET_SCALES,
        "coil": COIL_SCALES,
        "mutability": MUTABILITY_SCALES,
    }

# /-- Formating --/
class CompactDict(dict):
    pass


class CustomEncoder(json.JSONEncoder):
    def iterencode(self, obj, _one_shot=False):
        for chunk in super().iterencode(obj, _one_shot=_one_shot):
            yield chunk

    def default(self, obj):
        if isinstance(obj, CompactDict):
            return json.loads(json.dumps(obj, separators=(",", ":")))
        return super().default(obj)


# /-- Scale Calculation --/
def calculateFeatureScale(feature):
    scales = []

    for scale in SCALE_IDS[feature]:
        scales.append(normalizeScale(getDict(scale)))

    values = {}
    for aa in AMINO_ORDER:
        aaSum = 0
        for scale in scales:
            aaSum += scale[aa]
        values.update({aa: aaSum/len(scales)})

    return CompactDict({aa: float(values[aa]) for aa in AMINO_ORDER})

def normalizeScale(scale):
    maxValue = max(abs(v) for v in scale.values())

    return {aa: v / maxValue for aa, v in scale.items()}

def getDict(id: str):
    values = aaindex1[id].values
    
    return {aa: float(values[aa]) for aa in AMINO_ORDER}


# /-- Main --/
def main():
    featureScales = {}
    for feature in SCALE_IDS:
        featureScales.update({feature: calculateFeatureScale(feature)})

    outputFile = "scales.json"
    with open(outputFile, "w") as f:
        f.write("{\n")

        items = list(featureScales.items())
        for i, (scale, aa_map) in enumerate(items):
            compact_inner = json.dumps(aa_map, separators=(", ", ":"))
            comma = ", " if i < len(items) - 1 else ""
            f.write(f'  "{scale}": {compact_inner}{comma}\n')

        f.write("}\n")


    print(f"Done! All Results written!")

if __name__ == "__main__":
    main()

