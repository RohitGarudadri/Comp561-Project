from operator import indexOf
import os
import operator
import json
import datetime
from math import log10, isinf, ceil
from itertools import repeat
from collections import defaultdict
__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))
# Very small valued number for log(0) equivalent
zero = 10**-200


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
    total_score = 0
    # Convert data to dictionary of probabilities based on index value
    for i in range(len(genome)):
        probability = float(probabilities[i])
        total_score += probability
        remaining_prob = (1-probability)/3
        nucleotide = genome[i].upper()
        prob_dict = {}
        prob_dict[nucleotide] = probability
        for nucleotide in nucleotides.difference(set((nucleotide))):
            prob_dict[nucleotide] = remaining_prob
        probability_dictionary[str(i)] = prob_dict
    gap_score = total_score / len(genome)
    probability_dictionary["GAP"] = gap_score

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
        max_prob8 = 0
        score8 = 0
        for j in range(i, i+8):
            val = probability_dictionary[str(j)][genome[j]]
            max_prob8 += Log(val)
            score8 += val
        size_8_mers[seq8].append([i, max_prob8, score8])
        if (i < len(genome) - 11):
            seq11 = genome[i:i+11]
            max_prob11 = 0
            score11 = 0
            for j in range(i, i+11):
                val = probability_dictionary[str(j)][genome[j]]
                max_prob11 += Log(val)
                score11 += val
            size_11_mers[seq11].append([i, max_prob11, score11])
        if (i < len(genome) - 15):
            seq15 = genome[i:i+15]
            max_prob15 = 0
            score15 = 0
            for j in range(i, i+15):
                val = probability_dictionary[str(j)][genome[j]]
                max_prob15 += Log(val)
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


def Log(x):
    if x <= 0:
        return log10(zero)
    else:
        return log10(x)


def score(nucleotide, index):
    vals = probability_dictionary[str(index)]
    best_nucleotide = max(vals.items(), key=operator.itemgetter(1))[0]
    if nucleotide == best_nucleotide:
        return vals[nucleotide], Log(vals[nucleotide])
    else:
        return vals[nucleotide] - vals[best_nucleotide], Log(vals[nucleotide])


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
                probability = left_probability + match[1] + right_probability
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
    cur_prob = 0
    prob_at_max = 0
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
        cur_prob += vals[1]
        if cur_score > max_score:
            max_score = cur_score
            max_index = start
            max_genome_index = genome_start
            prob_at_max = cur_prob
            drop = 0
        else:
            drop = max_score - cur_score

    return max_index, max_genome_index, prob_at_max, max_score


def gappedExtension(query, ungapped_matches):
    # Need to perform NW variant on the right and left sides of the remaining query sequence.
    alignment_info = []
    for k, v in ungapped_matches.items():
        for i, match in enumerate(v):
            # print(k, i)
            left_index, genome_index, cur_prob, cur_score = tuple(match)
            left_query = query[:left_index]
            left_genome_index = genome_index - 1
            # if k == "CCTCATCCCAC":
            #     print(match)
            #     print(left_query)
            #     print(len(left_query))
            left_info = align(left_query, left_genome_index, True)
            if isinf(left_info[1]):
                continue
            right_query = query[left_index + len(k):]
            right_genome_index = genome_index + len(k)
            right_info = align(right_query, right_genome_index)
            if isinf(right_info[1]):
                continue

            total_score = left_info[2] + cur_score + right_info[2]
            total_prob = left_info[1] + cur_prob + right_info[1]
            genome_start_index = genome_index - len(left_info[0])
            alignment = left_info[0] + k + right_info[0]
            alignment_info.append(
                (alignment, total_score, total_prob, genome_start_index))
    return sorted(alignment_info, key=lambda info: info[1], reverse=True)[:5]


def align(query, index, reverse=False):
    n = len(query)
    if reverse:
        s1 = query[::-1]
        new_index = max(index -
                        (2*ceil(float(probability_dictionary["GAP"]) * n) + 1), 0)
        s2 = genome[index: new_index: -1]
        index = new_index
    else:

        s1 = query
        s2 = genome[index: index + 2 *
                    ceil(float(probability_dictionary["GAP"]) * n)+1]
    m = len(s2)
    if m < n:
        print("Not possible")
        return ("", float("-inf"), 0)
    # if query == "GAGGGGAGAGCGGGCCCGGCGGGCTCCGGGAGGAGGTGG":
    #     print(s2)
    # Create empty dp array. Structured as (score, probability, prev_x, prev_y)
    M = [[(0, 0, 0, 0) for _ in range(m + 1)] for _ in range(n + 1)]
    # Can't align query to gaps in genome so those indices are marked as impossible
    for i in range(1, n + 1):
        for j in range(i):
            M[i][j] = (float("-inf"), 0,  0, 0)
    for i in range(m - n):
        for j in range(i+n+1, m + 1):
            M[i][j] = (float("-inf"), 0, 0, 0)

    for i in range(1, n + 1):
        for j in range(i, m-n+i+1):
            s_diag = score(s1[i-1], index+j-1)
            diagonal = M[i-1][j-1][0] + s_diag[0]
            left = M[i][j-1][0] - probability_dictionary["GAP"]
            if diagonal >= left:
                M[i][j] = (diagonal, M[i-1][j-1][1] + s_diag[1], i-1, j-1)
            else:
                # Temporarily set probability to s_diag[1], but it needs to change
                M[i][j] = (left, M[i][j-1][1] + s_diag[1], i, j-1)

    # Now we need the highest score in the bottom row.
    index = M[-1].index(max(M[-1]))
    max_score = M[-1][index][0]
    max_prob = M[-1][index][1]
    s1_final = []
    i = n
    j = index
    while i > 0 and j > 0:
        vals = M[i][j]
        if vals[-2] == i-1 and vals[-1] == j-1:
            if reverse:
                s1_final.append(s1[i-1])
            else:
                s1_final.insert(0, s1[i-1])
            i -= 1
            j -= 1
        else:
            if reverse:
                s1_final.append("_")
            else:
                s1_final.insert(0, "_")
            j -= 1
    s1_final = "".join(s1_final)
    return (s1_final, max_prob, max_score)


print(datetime.datetime.now())
genome, probability_dictionary = parseGenome()
size_8_mers, size_11_mers, size_15_mers = preprocess()
# print(len(probability_dictionary))
sequence_files = ["query_seq1000.txt"]
for file in sequence_files:
    new_file = file.split(".")[0]+".json"
    with open(os.path.join(__location__, "Query Seqs", file), "r") as seq_file:
        print(file)
        seqs = []
        for line in seq_file:
            seqs.append(line.strip().split(" "))
        seq_alignments = {}
        for i in range(len(seqs)):
            print(datetime.datetime.now())

            print("i= ", i)
            alignments = {}
            for j in range(3):
                print("j= ", j)
                query = seqs[i][j+1]
                stripped_seq = "".join(query.split("_"))
                matches = shortPerfectMatch(stripped_seq)
                # # print(matches)
                ungapped_matches = ungappedExtension(stripped_seq, matches)
                # print(ungapped_matches)
                final_matches = gappedExtension(stripped_seq, ungapped_matches)
                # print(final_matches)
                alignments[query] = final_matches
            seq_alignments[seqs[i][0]] = alignments
    with open(os.path.join(__location__, "Query Predictions", new_file), "w") as align_file:
        align_file.write(json.dumps(seq_alignments))
# query = "GTCTCTAGCTGTCCGAGGCTTGGGCAACCACCTAGCCTCTATGAGCCTTAGGTGACACATAGGTAGAAATCATAGTACTT_CCTGATCTAGGGTTGAGATTAA_AATGCCTGGAATACCTGATC_GGCAGATGGAGAATTCCACACAGTCCTTGGCACACAGTAGG_CTTCAGCGGATATTATCCGTGTAGATTCATATTTCCTGGCACAAGTTCAGTGTCTCCACCTCATCCCACAGTTGACTCAGGATCTGATGATGCCAGTTCACAGTTCTCTAG_CCTCTTTTCTACCCAAACCCCAAACCTCGTTTAGGACCTTTCATTGCTCATAATGCAGGTCGCCACTTCATTGGCAGCACCCTTAGGAAGATTTGTTGAGGGAAGGTTGAAATTTATAGGAAAATTGACGTGTCT_CAGTGTCTGTGGAGCTGTCAGGTCCTCTGGCCAC_GGTGGAGGATGGGGTGTGCTTTCGCCTGGCCCGGAGCCCAGGCCTGCCGTCATGAGAAGATGGAGT"
# stripped_seq = "".join(query.split("_"))
# matches = shortPerfectMatch(stripped_seq)
# # # print(matches)
# ungapped_matches = ungappedExtension(stripped_seq, matches)
# # print(ungapped_matches)
# final_matches = gappedExtension(stripped_seq, ungapped_matches)
# print(final_matches)
print(datetime.datetime.now())
