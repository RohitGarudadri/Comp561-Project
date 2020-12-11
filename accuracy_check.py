import os
import math
import json
from collections import defaultdict
__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))
with open(os.path.join(__location__, 'Genome_probabilities.json')) as json_file:
    data = json.load(json_file)
    f = open(os.path.join(__location__, "Accuracies",
                          "accuracies500.txt"), "a")
    with open(os.path.join(__location__, "Query Predictions", 'query_seq500.json')) as json_file:
        seqs = json.load(json_file)

        base_seq_num = 0

        for og_seq in seqs:  # og_seq will hold the original sequence, so loops 4 times per file
            mod_seqs = seqs[og_seq]
            base_seq_num = base_seq_num + 1
            mod_num = 0

            for modded in mod_seqs:  # modded holds each modded sequence as a string loops 3 times
                list_of_objects = mod_seqs[modded]
                mod_num = mod_num + 1
                alignment_num = 0
                mod_avg = 0
                seq_count = 0
                for object in list_of_objects:  # loops 5 times
                    alignment_num = alignment_num + 1
                    predicted_alignment = object[0]
                    predicted_index = object[3]
                    total = 0
                    accurate = 0
                    i = str(predicted_index)
                    # end of preprocessing

                    for nuc in predicted_alignment:  # loops through predicted nucs
                        curr = data[i]
                        keep = 0
                        hold = 'A'

                        for curr_nuc in curr:  # for loop going through curr index to find appropriate nuc
                            # curr_nuc is the letter of Nuc and curr[curr_nuc] gives the prob of tht nuc
                            val = curr[curr_nuc]
                            if (val > keep):
                                keep = val
                                hold = curr_nuc
                                # swaps to largest

                        # check if you're correct, update total and accurate as appropriate
                        if(hold == nuc):
                            accurate = accurate + 1
                        elif(nuc in ['A', 'C', 'G', 'T']):
                            # gives probability of the predicted nuc @ index in genome
                            predicted = curr[nuc]
                            predicted = predicted / keep
                            accurate = accurate + predicted  # add this in instead of 0 as it is still possible

                        total = total + 1

                        # update
                        predicted_index = predicted_index + 1
                        i = str(predicted_index)

                    fin_accuracy = (accurate / total)*100
                    seq_count += 1
                    mod_avg = mod_avg + (fin_accuracy/100)
                    f.write('%d.' % base_seq_num)
                    f.write('%d.' % mod_num)
                    f.write('%d Accuracy is: ' % alignment_num)
                    f.write('%d' % fin_accuracy)
                    f.write('%\n')
                mod_avg = (mod_avg / seq_count) * 100
                f.write('The average accuracy for %d.' % base_seq_num)
                f.write('%d is: ' % mod_num)
                f.write('%d' % mod_avg)
                f.write('%\n')
    f.close()
