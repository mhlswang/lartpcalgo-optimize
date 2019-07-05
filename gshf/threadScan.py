import os

nrep = 5

os.environ['OMP_NESTED'] = 'TRUE'
print(os.environ['OMP_NESTED'])

os.environ['KMP_AFFINITY'] = 'scatter'
#os.environ['KMP_AFFINITY'] = 'compact'
print(os.environ['KMP_AFFINITY'])

vnthevt1  = ['1,1','1,2','1,4','1,6','1,10','1,15','1,20','1,25','1,30','1,40','1,50','1,60','1,80','1,100']
vnthroi1  = ['2,1','4,1','6,1','10,1','15,1','20,1','25,1','30,1','40,1','50,1','60,1','80,1','100,1']
vnthevt5  = ['5,2','5,4','5,6','5,10','5,15','5,20']
vnthevt10 = ['10,2','10,4','10,6','10,10']
vnthevt20 = ['20,2','20,3','20,4','20,5']
addtnl = ['1,150','1,200','150,1','200,1','5,30','5,40','10,15','10,20','20,7','20,10']
addtnl2 = ['70,1','90,1','1,70','1,90','5,13','5,16','10,8']
cmpct = ['10,8','1,80','80,1','1,10','10,1','5,6','20,3']

vnthstr = []
vnthstr.extend(vnthevt1)
vnthstr.extend(vnthroi1)
vnthstr.extend(vnthevt5)
vnthstr.extend(vnthevt10)
vnthstr.extend(vnthevt20)
vnthstr.extend(addtnl2)
vnthstr.extend(cmpct)

print(vnthstr)

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
