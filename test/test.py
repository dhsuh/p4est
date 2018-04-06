import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt

outfile = open('temp','rw+')
final = []
proc = []
width = 0.35
for i in range(1,10):
    o = i*16
    proc.append(o)
    temp = "mpirun -np " + str(o) + " " + "./p4est_test_fusion"
    subprocess.call(temp, shell=True, stdout = outfile)
    outfile = open('temp','rw+')
    last = outfile.readlines()[-2]
    last = last[19:-3].split()
    last = [float(i) for i in last]
    final.append(last)

ind = range(1,10)
final = np.array(final).T.tolist()
clr = ['blue', 'green', 'yellow', 'magenta', 'white', 'cyan', 'red']
print(final)
#for i in range(len(final)):
p1 = plt.bar(ind, final[0], width, align='center', color=clr[0])
p2 = plt.bar(ind, final[1], width, align='center', color=clr[1], bottom = final[0])
p3 = plt.bar(ind, final[2], width, align='center', color=clr[2], bottom = [final[0][j] + final[1][j] for j in range(len(final[0]))])
p4 = plt.bar(ind, final[3], width, align='center', color=clr[3], bottom = [final[0][j] + final[1][j] + final[2][j] for j in range(len(final[0]))])
p5 = plt.bar(ind, final[4], width, align='center', color=clr[4], bottom = [final[0][j] + final[1][j] + final[2][j] + final[3][j] for j in range(len(final[0]))])
p6 = plt.bar(ind, final[5], width, align='center', color=clr[5], bottom = [final[0][j] + final[1][j] + final[2][j] + final[3][j] + final[4][j] for j in range(len(final[0]))])
p7 = plt.bar(ind, final[6], width, align='center', color=clr[6], bottom = [final[0][j] + final[1][j] + final[2][j] + final[3][j] + final[4][j] + final[5][j] for j in range(len(final[0]))])

    
plt.xticks(ind, [str(i) for i in proc])
plt.yticks(np.arange(0,1.5,0.1))
plt.ylabel('Seconds')
plt.xlabel('Procs')
plt.legend((p1[0],p2[0],p3[0],p4[0],p5[0],p6[0],p7[0]), ('Full Loop', 'Local Coarsening', 'Local Refinement', 'Balance', 'Partition', 'Ghost', 'Optimized'), loc='upper left')
plt.show()
