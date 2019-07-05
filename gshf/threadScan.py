import os

nrep = 5

os.environ['OMP_NESTED'] = 'TRUE'
print(os.environ['OMP_NESTED'])

os.environ['KMP_AFFINITY'] = 'scatter'
#os.environ['KMP_AFFINITY'] = 'compact'
print(os.environ['KMP_AFFINITY'])

vnthevt1  = ['1,1','1,2','1,4','1,6','1,10','1,15','1,20','1,25','1,30','1,40','1,50','1,60','1,70','1,80','1,90','1,100','1,150','1,200']
vnthroi1  = ['2,1','4,1','6,1','10,1','15,1','20,1','25,1','30,1','40,1','50,1','60,1','70,1','80,1','90,1','100,1','150,1','200,1']
vnthevt5  = ['5,2','5,4','5,6','5,10','5,13','5,15','5,16','5,20','5,30','5,40']
vnthevt10 = ['10,2','10,4','10,6','10,8','10,10','10,15','10,20']
vnthevt20 = ['20,2','20,3','20,4','20,5','20,7','20,10']

vnthstr = []
vnthstr.extend(vnthevt1)
vnthstr.extend(vnthroi1)
vnthstr.extend(vnthevt5)
vnthstr.extend(vnthevt10)
vnthstr.extend(vnthevt20)

print(vnthstr)

nte = []
ntw = []
tim = []

for thstr in vnthstr:
    os.environ['OMP_NUM_THREADS'] = thstr
    avg = 0.
    for it in range(0,nrep):
        os.system('./gshf-mrqdt3 | grep without >> logtmp.txt')
    with open('logtmp.txt') as fp:
        for line in fp:
            time = float(line.split('=')[1].replace('\n',''))
            #print(time)
            avg = avg + time
    os.system('rm logtmp.txt') 
    avg = avg/float(nrep)
    print(os.environ['OMP_NUM_THREADS'],avg)
    nte.append(thstr.split(',')[0])
    ntw.append(thstr.split(',')[1])
    tim.append(avg)

print(nte)
print(ntw)
print(tim)
