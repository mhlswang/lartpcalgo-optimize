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



if __name__ == '__main__':
	
	f_known = open("og_result.txt", "r")
	f_new   = open("mkl_results.txt", "r")

	known_li = []
	new_li   = []

	for line in f_known:
		known_li.append([float(x) for x in line.split()])

	for line in f_new:
		new_li.append([float(x) for x in line.split()])

	for old,new in zip(known_li,new_li):
		out = ""
		i = 0
		for o,n in zip(old,new):
			if i > 2:
				out = out + str((o-n)/o) + ','
			i = i + 1 

		print out

	f_new.close()
	f_known.close()



