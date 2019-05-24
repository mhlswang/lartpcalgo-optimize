



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
		for o,n in zip(old,new):
			out = out + str(o-n) + ','
		print out

	f_new.close()
	f_known.close()



