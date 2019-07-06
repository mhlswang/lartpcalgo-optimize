import os

nrep = 5

os.environ['OMP_NESTED'] = 'TRUE'
print(os.environ['OMP_NESTED'])

os.environ['KMP_AFFINITY'] = 'scatter'
#os.environ['KMP_AFFINITY'] = 'compact'
print(os.environ['KMP_AFFINITY'])

vnthevt1  = ['1,1','1,2','1,5','1,10','1,20','1,40','1,60','1,80','1,100','1,120']
vnthroi1  = ['2,1','5,1','10,1','20,1','40,1','60,1','80,1','100,1','120,1']
vnthevt5  = ['5,2','5,4','5,8','5,12','5,16','5,20','5,24']
vnthevt10 = ['10,2','10,4','10,6','10,8','10,10','10,12']
vnthevt20 = ['20,2','20,3','20,4','20,5','20,6']

vnthstr = []
vnthstr.append(vnthevt1)
vnthstr.append(vnthroi1)
vnthstr.append(vnthevt5)
vnthstr.append(vnthevt10)
vnthstr.append(vnthevt20)

print(vnthstr)

for series in vnthstr:
    nte = []
    ntw = []
    tim = []
    for thstr in series:
        os.environ['OMP_NUM_THREADS'] = thstr
        avg = 0.
        for it in range(0,nrep):
            os.system('./gshf-mrqdt3 | grep without >> logtmp.txt')
        with open('logtmp.txt') as fp:
            for line in fp:
                time = float(line.split('=')[1].replace('\n',''))
                avg = avg + time
        os.system('rm logtmp.txt')
        avg = avg/float(nrep)
        print(os.environ['OMP_NUM_THREADS'],avg)
        nte.append(int(thstr.split(',')[0]))
        ntw.append(int(thstr.split(',')[1]))
        tim.append(float('%.2f'%(avg)))
    print('nte = ',nte)
    print('ntw = ',ntw)
    print('tim = ',tim)
