#!/usr/bin/env python3

import json
import csv
import gzip
import sys

pext_file = csv.DictReader(gzip.open(sys.argv[1],"rt"),delimiter="\t")
formated = open(sys.argv[2],"w")

for line in pext_file:

    line_dict = json.loads(line['tx_annotation'])

    for gene in line_dict:

        end = int(line['pos']) + 1

        formated.write(f"{line['chrom']}\t{line['pos']}\t{end}\t{line['ref']}\t{line['alt']}\t{gene['mean_proportion']}\t{gene['Brain_Cerebellum']}\t{gene['ensg']}\n")

formated.close()
