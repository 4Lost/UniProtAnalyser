import json

FILE_NAME='SCL-surpriseMetric'
#FILE_NAME='UniProt50-surpriseMetric'
#FILE_NAME='SCL-UniProt50-surpriseMetric'

with open(f'{FILE_NAME}.jsonl', 'r') as infile, open(f'{FILE_NAME}Filtered.jsonl', 'w') as outfile:
    for line in infile:
        data = json.loads(line.strip())
        if list(data.values())[0].get('class') == 'extremly surprising':
            json.dump(data, outfile)
            outfile.write('\n')
        print(list(data.values())[0])
