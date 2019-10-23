'''
make plots for fftw
'''


import numpy as np
import matplotlib.pyplot as plt

nwires = 5000 # or something
ntcks = 4096

line_num = 0
wire_num = 0
error_data = []

# filename   = 'fftw_output.csv'
filename   = 'mkl_output.csv'

f = open(filename, "r")
nnan = 0
for line in f:

	n = float(line.strip())
	if not np.isnan(n):
		error_data.append(n)
	else:
		nnan += 1

f.close()

print "number of nan = %d"%nnan
print "number good   = %d"%len(error_data)


num_bins = 50
plt.figure()
plt.hist(error_data, num_bins, log = True)
plt.xlabel('error')
plt.ylabel('count')
plt.title('Histogram of mkl errors')
plt.show()
plt.savefig('mkl.png')
