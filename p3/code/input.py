#!/usr/bin/python
#You can use the script to generate input instances for the problem
# Use: ./input.py <input_length_desired> <output_file_name>

import random
import sys

if len(sys.argv) != 3:
    print "Usage: ", sys.argv[0], " <length> <output_file>"
    exit()

#output file
f = open(sys.argv[2], 'w')

#length of the list needed
l = sys.argv[1]

f.write(l + "\n")
l = int(l)

for i in range(l):
    r = random.randint(1, 2*l)
    f.write(str(r) + "\n")

f.close()
