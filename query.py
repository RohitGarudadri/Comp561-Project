import os
import math
import json
import random
from collections import defaultdict

changelog = []
with open('Genome_probabilities.txt') as json_file:
    data = json.load(json_file)
    random.seed(50)
    i=0
    f = open("query_seq100.txt", "a")
    cl = open("change_log100.txt", "a")
    while(i<4):
        r_int = random.randint(0,len(data))
        start = r_int
        #print(r_int)
        s = str(r_int)
        query = []
        j = 0

        while(j<100): #results in query seq stored in query list.
            curr = data[s] #gives the index of the nucleotide and stores the list of nucleotides @ index
            keep = 0
            hold = 'A'

            for curr_nuc in curr: #for loop going through curr index to find appropriate nuc
                val = curr[curr_nuc] #curr_nuc is the letter of Nuc and curr[curr_nuc] gives the prob of tht nuc
                if(val>keep):
                    keep = val
                    hold = curr_nuc
                    #swaps to largest
            query.append(hold)
            j= j+ 1
            r_int = r_int +1
            s = str(r_int)
        #i = i+1
        #record base query seq
        for nuc in query:
            f.write('%s' % nuc)
        f.write(' ')

        #make copies for modification
        copy_one = []
        copy_two = []
        copy_three = []
        for nuc in query:
            copy_one.append(nuc)
            copy_two.append(nuc)
            copy_three.append(nuc)

        #first copy
        log = 0
        k=0
        for nuc in copy_one:
            if(random.random()>0.99424):
                copy_one[k] = '_'

                if(i==0):
                    cl.write('1.1: ')
                    cl.write('%d to _\n' % k)
                elif(i==1):
                    cl.write('2.1: ')
                    cl.write('%d to _\n' % k)
                elif(i==2):
                    cl.write('3.1: ')
                    cl.write('%d to _\n' % k)
                elif(i==3):
                    cl.write('4.1: ')
                    cl.write('%d to _\n' % k)

            elif(random.random()>0.99 and copy_one[k] != 'A'):
                copy_one[k] = 'A'

                if(i==0):
                    cl.write('1.1: ')
                    cl.write('%d to A\n' % k)
                elif(i==1):
                    cl.write('2.1: ')
                    cl.write('%d to A\n' % k)
                elif(i==2):
                    cl.write('3.1: ')
                    cl.write('%d to A\n' % k)
                elif(i==3):
                    cl.write('4.1: ')
                    cl.write('%d to A\n' % k)

            elif(random.random()>0.989876 and copy_one[k] != 'G'):
                copy_one[k] = 'G'

                if(i==0):
                    cl.write('1.1: ')
                    cl.write('%d to G\n' % k)
                elif(i==1):
                    cl.write('2.1: ')
                    cl.write('%d to G\n' % k)
                elif(i==2):
                    cl.write('3.1: ')
                    cl.write('%d to G\n' % k)
                elif(i==3):
                    cl.write('4.1: ')
                    cl.write('%d to G\n' % k)

            elif(random.random()>0.983609 and copy_one[k] != 'T'):
                copy_one[k] = 'T'

                if(i==0):
                    cl.write('1.1: ')
                    cl.write('%d to T\n' % k)
                elif(i==1):
                    cl.write('2.1: ')
                    cl.write('%d to T\n' % k)
                elif(i==2):
                    cl.write('3.1: ')
                    cl.write('%d to T\n' % k)
                elif(i==3):
                    cl.write('4.1: ')
                    cl.write('%d to T\n' % k)

            elif(random.random()>0.979863and copy_one[k] !='C'):
                copy_one[k] = 'C'

                if(i==0):
                    cl.write('1.1: ')
                    cl.write('%d to C\n' % k)
                elif(i==1):
                    cl.write('2.1: ')
                    cl.write('%d to C\n' % k)
                elif(i==2):
                    cl.write('3.1: ')
                    cl.write('%d to C\n' % k)
                elif(i==3):
                    cl.write('4.1: ')
                    cl.write('%d to C\n' % k)

            elif(random.random()>0.976235):
                copy_one.insert(k,'A')

                if(i==0):
                    cl.write('1.1: ')
                    cl.write('%d to +\n' % k)
                elif(i==1):
                    cl.write('2.1: ')
                    cl.write('%d to +\n' % k)
                elif(i==2):
                    cl.write('3.1: ')
                    cl.write('%d to +\n' % k)
                elif(i==3):
                    cl.write('4.1: ')
                    cl.write('%d to +\n' % k)
            k = k + 1

        for nuc in copy_one:
            f.write('%s' % nuc)
        f.write(' ')

        #second copy
        k=0
        for nuc in copy_two:
            if(random.random()>0.99346):
                copy_two[k] = '_'

                if(i==0):
                    cl.write('1.2: ')
                    cl.write('%d to _\n' % k)
                elif(i==1):
                    cl.write('2.2: ')
                    cl.write('%d to _\n' % k)
                elif(i==2):
                    cl.write('3.2: ')
                    cl.write('%d to _\n' % k)
                elif(i==3):
                    cl.write('4.2: ')
                    cl.write('%d to _\n' % k)

            elif(random.random()>0.987341 and copy_two[k] != 'A'):
                copy_two[k] = 'A'

                if(i==0):
                    cl.write('1.2: ')
                    cl.write('%d to A\n' % k)
                elif(i==1):
                    cl.write('2.2: ')
                    cl.write('%d to A\n' % k)
                elif(i==2):
                    cl.write('3.2: ')
                    cl.write('%d to A\n' % k)
                elif(i==3):
                    cl.write('4.2: ')
                    cl.write('%d to A\n' % k)

            elif(random.random()>0.984218 and copy_two[k] != 'G'):
                copy_two[k] = 'G'

                if(i==0):
                    cl.write('1.2: ')
                    cl.write('%d to G\n' % k)
                elif(i==1):
                    cl.write('2.2: ')
                    cl.write('%d to G\n' % k)
                elif(i==2):
                    cl.write('3.2: ')
                    cl.write('%d to G\n' % k)
                elif(i==3):
                    cl.write('4.2: ')
                    cl.write('%d to G\n' % k)

            elif(random.random()>0.979543 and copy_two[k] != 'T'):
                copy_two[k] = 'T'

                if(i==0):
                    cl.write('1.2: ')
                    cl.write('%d to T\n' % k)
                elif(i==1):
                    cl.write('2.2: ')
                    cl.write('%d to T\n' % k)
                elif(i==2):
                    cl.write('3.2: ')
                    cl.write('%d to T\n' % k)
                elif(i==3):
                    cl.write('4.2: ')
                    cl.write('%d to T\n' % k)

            elif(random.random()>0.978453 and copy_two[k] !='C'):
                copy_two[k] = 'C'

                if(i==0):
                    cl.write('1.2: ')
                    cl.write('%d to C\n' % k)
                elif(i==1):
                    cl.write('2.2: ')
                    cl.write('%d to C\n' % k)
                elif(i==2):
                    cl.write('3.2: ')
                    cl.write('%d to C\n' % k)
                elif(i==3):
                    cl.write('4.2: ')
                    cl.write('%d to C\n' % k)

            elif(random.random()>0.97242):
                copy_two.insert(k,'G')

                if(i==0):
                    cl.write('1.2: ')
                    cl.write('%d to +\n' % k)
                elif(i==1):
                    cl.write('2.2: ')
                    cl.write('%d to +\n' % k)
                elif(i==2):
                    cl.write('3.2: ')
                    cl.write('%d to +\n' % k)
                elif(i==3):
                    cl.write('4.2: ')
                    cl.write('%d to +\n' % k)
            k = k + 1

        for nuc in copy_two:
            f.write('%s' % nuc)
        f.write(' ')

        #third copy
        k = 0
        for nuc in copy_three:
            if (random.random() > 0.99236):
                copy_three[k] = '_'

                if(i==0):
                    cl.write('1.3: ')
                    cl.write('%d to _\n' % k)
                elif(i==1):
                    cl.write('2.3: ')
                    cl.write('%d to _\n' % k)
                elif(i==2):
                    cl.write('3.3: ')
                    cl.write('%d to _\n' % k)
                elif(i==3):
                    cl.write('4.3: ')
                    cl.write('%d to _\n' % k)

            elif (random.random() > 0.989346 and copy_three[k] != 'A'):
                copy_three[k] = 'A'

                if(i==0):
                    cl.write('1.3: ')
                    cl.write('%d to A\n' % k)
                elif(i==1):
                    cl.write('2.3: ')
                    cl.write('%d to A\n' % k)
                elif(i==2):
                    cl.write('3.3: ')
                    cl.write('%d to A\n' % k)
                elif(i==3):
                    cl.write('4.3: ')
                    cl.write('%d to A\n' % k)

            elif (random.random() > 0.9849345 and copy_three[k] != 'G'):
                copy_three[k] = 'G'

                if(i==0):
                    cl.write('1.3: ')
                    cl.write('%d to G\n' % k)
                elif(i==1):
                    cl.write('2.3: ')
                    cl.write('%d to G\n' % k)
                elif(i==2):
                    cl.write('3.3: ')
                    cl.write('%d to G\n' % k)
                elif(i==3):
                    cl.write('4.3: ')
                    cl.write('%d to G\n' % k)

            elif (random.random() > 0.9773123 and copy_three[k] != 'T'):
                copy_three[k] = 'T'

                if(i==0):
                    cl.write('1.3: ')
                    cl.write('%d to T\n' % k)
                elif(i==1):
                    cl.write('2.3: ')
                    cl.write('%d to T\n' % k)
                elif(i==2):
                    cl.write('3.3: ')
                    cl.write('%d to T\n' % k)
                elif(i==3):
                    cl.write('4.3: ')
                    cl.write('%d to T\n' % k)

            elif (random.random() > 0.97431 and copy_three[k] != 'C'):
                copy_three[k] = 'C'

                if(i==0):
                    cl.write('1.3: ')
                    cl.write('%d to C\n' % k)
                elif(i==1):
                    cl.write('2.3: ')
                    cl.write('%d to C\n' % k)
                elif(i==2):
                    cl.write('3.3: ')
                    cl.write('%d to C\n' % k)
                elif(i==3):
                    cl.write('4.3: ')
                    cl.write('%d to C\n' % k)

            elif (random.random() > 0.9716789):
                copy_three.insert(k, 'T')

                if(i==0):
                    cl.write('1.3: ')
                    cl.write('%d to +\n' % k)
                elif(i==1):
                    cl.write('2.3: ')
                    cl.write('%d to +\n' % k)
                elif(i==2):
                    cl.write('3.3: ')
                    cl.write('%d to +\n' % k)
                elif(i==3):
                    cl.write('4.3: ')
                    cl.write('%d to +\n' % k)
            k = k + 1
        for nuc in copy_three:
            f.write('%s' % nuc)

        f.write(' %d' % start)
        f.write(' %d' % r_int)
        f.write('\n')
        i = i + 1
    cl.write('seeded w 50')

    f.close()
    cl.close()