import os
import math
import json
from collections import defaultdict

with open('Genome_probabilities.txt') as json_file:
    data = json.load(json_file)
    f = open("accuracies100_8.txt", "a")
    with open('query_seq100_8.json') as json_file:
        seqs = json.load(json_file)

        base_seq_num = 0

        for og_seq in seqs:   #og_seq will hold the original sequence, so loops 4 times per file
            mod_seqs = seqs[og_seq]
            base_seq_num = base_seq_num + 1
            mod_num = 0

            for modded in mod_seqs:   #modded holds each modded sequence as a string loops 3 times
                if(modded == 'Correct'):
                    f.write('The original index is : %d\n\n' % mod_seqs[modded])
                    continue
                list_of_objects = mod_seqs[modded]
                mod_num = mod_num + 1
                alignment_num = 0
                mod_avg = 0


                for object in list_of_objects:  # loops 5 times
                    alignment_num = alignment_num + 1
                    predicted_alignment = object[0]
                    predicted_index = object[3]
                    total = 0
                    accurate = 0
                    i = str(predicted_index)
                    f.write('The predicted index was: %d\n' % predicted_index)
                    ##end of preprocessing

                    for nuc in predicted_alignment: #loops through predicted nucs
                        curr = data[i]
                        keep = 0
                        hold = 'A'

                        for curr_nuc in curr:  # for loop going through curr index to find appropriate nuc
                            val = curr[curr_nuc]  # curr_nuc is the letter of Nuc and curr[curr_nuc] gives the prob of tht nuc
                            if (val > keep):
                                keep = val
                                hold = curr_nuc
                                # swaps to largest

                        ##check if you're correct, update total and accurate as appropriate
                        if(hold == nuc):
                            accurate = accurate + 1
                        elif(nuc in ['A','C','G','T']):
                            predicted = curr[nuc]  ##gives probability of the predicted nuc @ index in genome
                            predicted = predicted / keep
                            accurate = accurate + predicted ##add this in instead of 0 as it is still possible

                        total = total + 1

                        ##update
                        predicted_index = predicted_index + 1
                        i = str(predicted_index)

                    fin_accuracy = (accurate / total)*100
                    mod_avg = mod_avg + (fin_accuracy/100)
                    f.write('%d.' % base_seq_num)
                    f.write('%d.' % mod_num)
                    f.write('%d Accuracy is: ' % alignment_num)
                    f.write('%d' % fin_accuracy)
                    f.write('%\n\n')
                mod_avg = (mod_avg / 5) * 100
                f.write('The average accuracy for %d.' % base_seq_num)
                f.write('%d is: ' % mod_num)
                f.write('%d' % mod_avg)
                f.write('%\n\n')
    f.close()