import os
import math
import json
__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))


def parseGenome(genome_file, probability_file):
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
    return probability_dictionary


probability_file = input("Please enter the probability file name: ")
genome_file = input("Please enter the genome file name: ")

parseGenome(genome_file, probability_file)
