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

	old = ResultData("og_result.txt")
	new = ResultData("mkl_results.txt")


	old_rms   = old.get_column("rd_vec.rms")
	old_mytck = old.get_column("od_vec.mytck")
	new_rms   = new.get_column("rd_vec.rms")
	new_mytck = new.get_column("od_vec.mytck")
 	diff_rms   = [o-n for o,n in zip(old_rms,new_rms)]
 	diff_mytck = [o-n for o,n in zip(old_mytck,new_mytck)]

	x = range(len(old_rms))

	plt.figure(1)
	plt.plot(x, new_rms, 'bs',  ms=0.7)
	plt.plot(x, old_rms, 'r^', ms=0.2)
	plt.savefig('pos0.png')

	plt.figure(2)
	plt.plot(x, new_mytck, 'bs',  ms=0.7)
	plt.plot(x, old_mytck, 'r^', ms=0.2)
	plt.savefig('pos1.png')
 
	plt.figure(3)
	plt.plot(x, diff_mytck, 'bs',  ms=0.7)
	plt.plot(x, diff_rms,   'r^', ms=0.2)
	plt.savefig('pos2.png')





