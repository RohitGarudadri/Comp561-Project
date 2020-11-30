import os
import math
import json
from collections import defaultdict
__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))


def parseGenome():
    probability_file = "Probabilities.txt"
    genome_file = "Genome.txt"
    if os.path.exists(os.path.join(__location__, "Genome_probabilities.json")):
        with open(os.path.join(__location__, "Genome_probabilities.json"), "r") as my_file:
            probability_dictionary = json.load(my_file)
        with open(os.path.join(__location__, genome_file), "r") as my_file:
            genome = my_file.readline().strip()
        return genome, probability_dictionary

    # Extract data from files
    with open(os.path.join(__location__, genome_file),  "r") as my_file:
        genome = my_file.readline().strip()

    with open(os.path.join(__location__, probability_file), "r") as my_file:
        probabilities = my_file.readline().strip().split(" ")

    probability_dictionary = {}
    nucleotides = set(("A", "C", "G", "T"))
    # Convert data to dictionary of probabilities based on index value
    for i in range(len(genome)):
        probability = float(probabilities[i])
        remaining_prob = (1-probability)/3
        nucleotide = genome[i].upper()
        prob_dict = {}
        prob_dict[nucleotide] = probability
        for nucleotide in nucleotides.difference(set((nucleotide))):
            prob_dict[nucleotide] = remaining_prob
        probability_dictionary[str(i)] = prob_dict

    # Convert dictionary to json for easy transfer
    prob_dict_json = json.dumps(probability_dictionary)
    with open(os.path.join(__location__, "Genome_probabilities.json"), "w") as out_file:
        out_file.write(prob_dict_json)
    print("JSON file of probabilistic genome created at: ",
          os.path.join(__location__, "Genome_probabilities.json"))

    return genome, probability_dictionary


def preprocess():
    if os.path.exists(os.path.join(__location__, "Genome_8mers.json")):
        with open(os.path.join(__location__, "Genome_8mers.json"), "r") as my_file:
            size_8_mers = json.load(my_file)
        with open(os.path.join(__location__, "Genome_11mers.json"), "r") as my_file:
            size_11_mers = json.load(my_file)
        with open(os.path.join(__location__, "Genome_15mers.json"), "r") as my_file:
            size_15_mers = json.load(my_file)
        return size_8_mers, size_11_mers, size_15_mers
    size_15_mers = defaultdict(list)
    size_11_mers = defaultdict(list)
    size_8_mers = defaultdict(list)
    for i in range(len(genome)-8):
        seq8 = genome[i:i+8]
        max_prob8 = 0
        for j in range(i, i+8):
            max_prob8 += probability_dictionary[str(j)][genome[j]]
        max_prob8 /= 8
        size_8_mers[seq8].append((i, max_prob8))
        if (i < len(genome) - 11):
            seq11 = genome[i:i+11]
            max_prob11 = 0
            for j in range(i, i+11):
                max_prob11 += probability_dictionary[str(j)][genome[j]]
            max_prob11 /= 11
            size_11_mers[seq11].append((i, max_prob11))
        if (i < len(genome) - 15):
            seq15 = genome[i:i+15]
            max_prob15 = 0
            for j in range(i, i+15):
                max_prob15 += probability_dictionary[str(j)][genome[j]]
            max_prob15 /= 15
            size_15_mers[seq15].append((i, max_prob15))

    with open(os.path.join(__location__, "Genome_8mers.json"), "w") as out_file:
        out_file.write(json.dumps(size_8_mers))
    with open(os.path.join(__location__, "Genome_11mers.json"), "w") as out_file:
        out_file.write(json.dumps(size_11_mers))
    with open(os.path.join(__location__, "Genome_15mers.json"), "w") as out_file:
        out_file.write(json.dumps(size_15_mers))

    return (size_8_mers, size_11_mers, size_15_mers)


genome, probability_dictionary = parseGenome()
size_8_mers, size_11_mers, size_15_mers = preprocess()
