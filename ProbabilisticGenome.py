import os
import operator
import json
from itertools import repeat
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
        return defaultdict(list, size_8_mers), defaultdict(list, size_11_mers), defaultdict(list, size_15_mers)
    size_15_mers = defaultdict(list)
    size_11_mers = defaultdict(list)
    size_8_mers = defaultdict(list)
    for i in range(len(genome)-8):
        seq8 = genome[i:i+8]
        max_prob8 = 1
        score8 = 0
        for j in range(i, i+8):
            val = probability_dictionary[str(j)][genome[j]]
            max_prob8 *= val
            score8 += val
        size_8_mers[seq8].append([i, max_prob8, score8])
        if (i < len(genome) - 11):
            seq11 = genome[i:i+11]
            max_prob11 = 1
            score11 = 0
            for j in range(i, i+11):
                val = probability_dictionary[str(j)][genome[j]]
                max_prob11 *= val
                score11 += val
            size_11_mers[seq11].append([i, max_prob11, score11])
        if (i < len(genome) - 15):
            seq15 = genome[i:i+15]
            max_prob15 = 1
            score15 = 0
            for j in range(i, i+15):
                val = probability_dictionary[str(j)][genome[j]]
                max_prob15 *= val
                score15 += val
            size_15_mers[seq15].append([i, max_prob15, score15])

    with open(os.path.join(__location__, "Genome_8mers.json"), "w") as out_file:
        out_file.write(json.dumps(size_8_mers))
    with open(os.path.join(__location__, "Genome_11mers.json"), "w") as out_file:
        out_file.write(json.dumps(size_11_mers))
    with open(os.path.join(__location__, "Genome_15mers.json"), "w") as out_file:
        out_file.write(json.dumps(size_15_mers))

    return (size_8_mers, size_11_mers, size_15_mers)


def shortPerfectMatch(query, w=11):
    if w == 8:
        genome_matches = size_8_mers
    elif w == 15:
        genome_matches = size_15_mers
    else:
        genome_matches = size_11_mers
    matches = defaultdict(list)
    for i in range(len(query)-w):
        seq = query[i:i+w]
        if seq in matches.keys():
            matches[seq][0].append(i)
        elif len(genome_matches[seq]) > 0:
            matches[seq] = [[i], genome_matches[seq]]
    return matches


def score(nucleotide, index):
    vals = probability_dictionary[str(index)]
    best_nucleotide = max(vals.items(), key=operator.itemgetter(1))[0]
    if nucleotide == best_nucleotide:
        return vals[nucleotide], vals[nucleotide]
    else:
        return vals[nucleotide] - vals[best_nucleotide], vals[nucleotide]


def ungappedExtension(query, matches):
    # Need another threshold value for determining what moves on to gapped phase
    ungapped_matches = defaultdict(list)
    for k, v in matches.items():
        for match in v[1]:
            for start_index in v[0]:
                right_index, right_genome_index, right_probability, right_maxima = extension(
                    k, start_index, query, match[0], 0)
                left_index, left_genome_index, left_probability, left_maxima = extension(k,
                                                                                         start_index, query, match[0], 1)
                extended_sequence = query[left_index: right_index+1]
                score = right_maxima + match[-1] + left_maxima
                probability = left_probability * match[1] * right_probability
                # if score > threshold_t:
                info = [left_index, left_genome_index, probability, score]
                if all(map(lambda x, y: x != y, [v[1] for v in ungapped_matches[extended_sequence]], repeat(left_genome_index))):
                    ungapped_matches[extended_sequence].append(info)
    return ungapped_matches


def extension(sequence, start, query, genome_start, dir):
    # dir variable for determining right vs left extension (0 = right, 1 = left)
    if dir == 0:
        genome_start += len(sequence) - 1
        start += len(sequence) - 1
    # Need a threshold value still. Use 5 for now as a default
    threshold = 5
    drop = 0
    max_score = 0
    max_index = start
    max_genome_index = genome_start
    cur_score = 0
    cur_prob = 1
    prob_at_max = 1
    while drop < threshold:
        if dir == 0:
            start += 1
            genome_start += 1
            if start > len(query)-1 or genome_start > len(genome) - 1:
                break
        else:
            start -= 1
            genome_start -= 1
            if start < 0 or genome_start < 0:
                break

        vals = score(query[start], genome_start)
        cur_score += vals[0]
        cur_prob *= vals[1]
        if cur_score > max_score:
            max_score = cur_score
            max_index = start
            max_genome_index = genome_start
            prob_at_max = cur_prob
            drop = 0
        else:
            drop = max_score - cur_score

    return max_index, max_genome_index, prob_at_max, max_score


# Possible probability normalization via taking nth root
#  of probability where n is the length of the sequence
# seems to make a difference, and math checks out for doing it only at end


genome, probability_dictionary = parseGenome()
size_8_mers, size_11_mers, size_15_mers = preprocess()
# og = "TAAGATGTTGGTATAAATACTCTGGAGCAA"
# query = "TAAGA_GTAGGTATAAATACTCAGGAGAAA"
# stripped_seq = "".join(query.split("_"))
# matches = shortPerfectMatch(stripped_seq)
# print(matches)
# ungapped_matches = ungappedExtension(stripped_seq, matches)
# print(ungapped_matches)
