'''
data order:
x (time in tcks)
y (points to fit)
f (initial guess)
f (final fuess)
'''
import numpy as np
import matplotlib.pyplot as plt

class EventData(object):
	"""docstring for EventData"""
	def __init__(self):
		self.data = {'x':[],'y':[],'fi':[],'fo':[]}
		

if __name__ == '__main__':


	filename   = 'fit.txt'
	num_ev     = 10
	event_data = [] # list of event data
	line_cnt   = 0
	line_li    = ['x','y','fi','fo']
	event_cnt  = 0

	for i in range(num_ev):
		event_data.append(EventData())

	f = open(filename, "r")

	for line in f:
		if (event_cnt < num_ev) and (len(line.split(',')) > 1):
			event_data[event_cnt].data[line_li[line_cnt]] = [float(x) for x in line.split(',')[:-1]]
			line_cnt = (line_cnt+1) % 4
			if line_cnt == 0:
				event_cnt += 1

	f.close()

	for i in range(num_ev):
		plt.figure()
		plt.plot(event_data[i].data['x'], event_data[i].data['y'],  'bo',  label="data to fit to")
		plt.plot(event_data[i].data['x'], event_data[i].data['fi'], 'g*-', label="initial guess")
		plt.plot(event_data[i].data['x'], event_data[i].data['fo'], 'r-',  label="final fit")
		plt.title("Fit of one Event for LArTPC")
		plt.xlabel("time in simulation")
		plt.ylabel("y value")
		plt.legend(loc='upper left')
		plt.savefig('fit_%d.png' % i)






