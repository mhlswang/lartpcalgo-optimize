'''
make plots for fftw
'''


import numpy as np
import matplotlib.pyplot as plt

nwires = 5000 # or something
ntcks = 4096

line_num = 0
wire_num = 0
expected_data = [[]]
computed_data = [[]]

filename   = 'fftw_output.csv.txt'
#filename   = 'mkl_output.csv.txt'

f = open(filename, "r")

for line in f:

	l = line.split(',')

	expected_data[wire_num].append(l[0])
	computed_data[wire_num].append(l[1])

	if (line_num % ntcks) == (ntcks-1):
		expected_data.append([])
		computed_data.append([])
		wire_num += 1

	line_num += 1

f.close()

