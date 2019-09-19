'''
ok so
fprintf(fout,"%d %d %d %lf %lf %lf %lf %lf\n",
od_vec[iv][ic_min].n,od_vec[iv][ic_min].imh,od_vec[iv][ic_min].ipp,
rd_vec[iv].simtck,rd_vec[iv].rectck,rd_vec[iv].rms,
od_vec[iv][ic_min].mytck,od_vec[iv][ic_min].mysigma)


rd_vec[iv].simtck is the simulated time of the hit
rd_vec[iv].rectck is the reconstructed time of the hit using larsoft
rd_vec[iv].rms is the spread of the hit using larsoft (not sure if it's simulated or reconstructed - but I can check)

od_vec[iv][ic_min].mytck is the  reconstructed time of the hit using  our fitter
od_vec[iv][ic_min].mysigma is the reconstructed spread of the gaussian using our fitter


so we can compare `rd_vec[iv].rms` to `od_vec[iv][ic_min].mysigma`
and `rd_vec[iv].simtck` compared to `rd_vec[iv].rectck` tells you how good larsoft is at finding the time
which then we can compare to our versions of the fitterr
by looking at `rd_vec[iv].rms` compared to `od_vec[iv][ic_min].mytck` for both marqfit and MKL

'''
import numpy as np
import matplotlib.pyplot as plt

NAMES = {
  "od_vec.n": 0,
  "od_vec.imh": 1,
  "od_vec.ipp": 2,
  "rd_vec.simtck": 3,
  "rd_vec.rectck": 4,
  "rd_vec.rms": 5,
  "od_vec.mytck": 6,
  "od_vec.mysigma": 7
}

class ResultData(object):
	"""docstring for ResultData"""
	def __init__(self, filename):
		
		self.filename = filename
		self.data_li = []
		self._load_data()


	def _load_data(self):

		f = open(self.filename, "r")

		for line in f:
			self.data_li.append([float(x) for x in line.split()])

		f.close()

	def get_column(self,col_name):
		try:
			li = [x[NAMES[col_name]] for x in self.data_li]	
		except:
			li = []
			print "ERROR: get_column; no such name"	

		return li


def fraction_different(self,old, new):
	for old_row,new_row in zip(old.data_li,new.data_li):
		out = ""
		i = 0
		for o,n in zip(old_row,new_row):
			if i > 2:
				out = out + str((o-n)/o) + ','
			i = i + 1 

		print out


if __name__ == '__main__':

	old = ResultData("marquart_results_512.txt")
	# old = ResultData("og_result.txt")
	
	new = ResultData("mkl_result.txt")


	# basic comparison of the results from our implementations
	old_sigma  = old.get_column("od_vec.mysigma")
	old_rms    = old.get_column("rd_vec.rms")
	old_mytck  = old.get_column("od_vec.mytck")

	new_sigma  = new.get_column("od_vec.mysigma")
	new_rms    = new.get_column("rd_vec.rms")
	new_mytck  = new.get_column("od_vec.mytck")

 	diff_sigma = [o-n for o,n in zip(old_sigma,new_sigma)]
 	diff_mytck = [o-n for o,n in zip(old_mytck,new_mytck)]

	x = range(len(old_sigma))

	diff_bins = 10 #len(old_sigma)/100

	plt.figure()
	plt.plot(x, new_sigma, 'bs',  ms=0.7)
	plt.plot(x, old_sigma, 'r^', ms=0.2)
	plt.xlabel('Line of Result file')
	plt.ylabel('od_vec.mysigma')
	plt.title('Comparison of MKL (Blue) and Marquardt (red) Sigma')
	plt.savefig('lart_sigma.png')

	plt.figure()
	plt.plot(x, new_mytck, 'bs',  ms=0.7)
	plt.plot(x, old_mytck, 'r^', ms=0.2)
	plt.xlabel('Line of Result file')
	plt.ylabel('od_vec.mytck')
	plt.title('Comparison of MKL (Blue) and Marquardt (red) mytck')
	plt.savefig('lart_mytck.png')
 
	plt.figure()
	plt.plot(x, diff_mytck, 'bs',  ms=0.7)
	plt.xlabel('Line of Result file')
	plt.ylabel('Difference')
	plt.title('Marquardt - MKL for mytck')
	plt.savefig('lart_diff_mytck.png')

	plt.figure()
	plt.hist(diff_mytck, diff_bins, log=True)
	plt.xlabel('Line of Result file')
	plt.ylabel('count')
	plt.title('Histogram of mytck difference values')
	plt.savefig('lart_hist_mytck.png')

	plt.figure()
	plt.plot(x, diff_sigma,   'r^', ms=0.2)
	plt.xlabel('Line of Result file')
	plt.ylabel('Difference')
	plt.title('Marquardt - MKL for sigma')
	plt.savefig('lart_diff_sigma.png')

	plt.figure()
	plt.hist(diff_mytck, diff_bins, log=True)
	plt.xlabel('Line of Result file')
	plt.ylabel('count')
	plt.title('Histogram of sigma difference values')
	plt.savefig('lart_hist_sigma.png')


	# comparison to LArSoft
	old_simtck = old.get_column("rd_vec.simtck")
	old_rectck = old.get_column("rd_vec.rectck")
	new_simtck = new.get_column("rd_vec.simtck")
	new_rectck = new.get_column("rd_vec.rectck")

	# rd_vec[iv].simtck-rd_vec[iv].rectck=resolution of the time of 
    # the fitter in larsoft compared to the truth
	old_larsoft_accuracy = [s-r for s,r in zip(old_simtck,old_rectck)]
	new_larsoft_accuracy = [s-r for s,r in zip(new_simtck,new_rectck)]

	acc_bins = 10 #len(old_sigma)/100
	
	plt.figure()
	plt.plot(x, new_larsoft_accuracy,   'bs', ms=0.2)
	plt.xlabel('Line of Result file')
	plt.ylabel('Difference')
	plt.title('LArTPC Accuracy (rd_vec.simtck-rd_vec.rectck)')
	plt.savefig('lart_new_lart_accuracy.png')

	plt.figure()
	plt.plot(x, old_larsoft_accuracy,   'r^', ms=0.2)
	plt.xlabel('Line of Result file')
	plt.ylabel('Difference')
	plt.title('LArTPC Accuracy (rd_vec.simtck-rd_vec.rectck)')
	plt.savefig('lart_old_lart_accuracy.png')


	plt.figure()
	plt.hist(old_larsoft_accuracy, acc_bins, range=(-500,500), log=True)
	plt.xlabel('histogram bins')
	plt.ylabel('counts')
	plt.title('LArTPC Accuracy Histogram')
	plt.savefig('lart_hist_lart_acc.png')

	# rd_vec[iv].simtck-od_vec[iv][ic_min].mytck = resolution of the
	# marqfit and MKL versions compared to truth
	marquat_accuracy  = [s-r for s,r in zip(old_simtck,old_mytck)]
	mkl_accuracy      = [s-r for s,r in zip(new_simtck,new_mytck)]


	plt.figure()
	plt.plot(x, marquat_accuracy,   'bs', ms=0.2)
	plt.xlabel('Line of Result file')
	plt.ylabel('Difference')
	plt.title('Marquardt Accuracy (rd_vec.simtck-od_vec.mytck)')
	plt.savefig('lart_maquat_accuracy.png')

	plt.figure()
	plt.hist(marquat_accuracy, acc_bins, range=(-500,500), log=True)
	plt.xlabel('histogram bins')
	plt.ylabel('counts')
	plt.title('Marquardt Accuracy Histogram')
	plt.savefig('lart_hist_maquat_acc.png')

	plt.figure()
	plt.plot(x, mkl_accuracy,   'r^', ms=0.2)
	plt.xlabel('Line of Result file')
	plt.ylabel('Difference')
	plt.title('MKL Accuracy (rd_vec.simtck-od_vec.mytck)')
	plt.savefig('lart_mkl_accuracy.png')

	plt.figure()
	plt.hist(mkl_accuracy, acc_bins, range=(-500,500), log=True)
	plt.xlabel('histogram bins')
	plt.ylabel('counts')
	plt.title('MKL Accuracy Histogram')
	plt.savefig('lart_hist_mkl_acc.png')

	# od_vec.sigma/rd_vec.rms = some other accuracy thing
	marquat_ratio  = [s/r for s,r in zip(old_sigma,old_rms)]
	mkl_ratio      = [s/r for s,r in zip(new_sigma,new_rms)]

	plt.figure()
	plt.plot(x, marquat_ratio,   'bs', ms=0.2)
	plt.xlabel('Line of Result file')
	plt.ylabel('Ratio')
	plt.title('Marquardt Ratio (od_vec.sigma/rd_vec.rms)')
	plt.savefig('lart_maquat_ratio.png')

	plt.figure()
	plt.plot(x, mkl_ratio,   'r^', ms=0.2)
	plt.xlabel('Line of Result file')
	plt.ylabel('Difference')
	plt.title('MKL Ratio (od_vec.sigma/rd_vec.rms)')
	plt.savefig('lart_mkl_ratio.png')


	plt.figure()
	colors = ['gray', 'violet', 'blue']
	labels = ['LArTPC', 'Marquardt', 'MKL']
	data = (old_larsoft_accuracy,marquat_accuracy,mkl_accuracy)
	plt.hist(data, acc_bins, color=colors, label=labels, log=True)
	plt.xlabel('histogram bins')
	plt.ylabel('counts')
	plt.legend(loc='upper left')
	plt.title('Combined Accuracy Histograms')

	plt.savefig('lart_comb_hist_acc.png')





