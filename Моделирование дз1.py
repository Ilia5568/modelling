from math import*
import os
import matplotlib.pyplot as plt

#вариант 8

#создание директории
if not os.path.isdir("results"):
     os.mkdir("results")

#программа

A=9.66459
def y(x):
    return -abs(sin(x)*cos(A)*exp(abs(1-sqrt(x**2+A**2)/pi)))
X=[]
Y=[]


os.chdir("results")
f=open('res1.xml','w')
f.write('<?xml version="1.1" encoding="UTF-8" ?>\n')
f.write('<data>\n\t')
f.write('<xdata>\n\t\t')


i=-10.0
while (i<10.1):
    
    X.append(i)
    Y.append(y(i))
    
    i=i+0.1
for i in X:
         f.write('<x>{:.1f}</x>\n\t\t'.format(i))
f.write('</xdata>\n\t')
f.write('<ydata>\n\t')
for i in Y:
         f.write('<y>{:.3f}</y>\n\t\t'.format(i))

f.write('</ydata>\n\t')
f.write('</data>\n')



f.close()
plt.plot(X,Y)
plt.show()
