import urllib.request
import math
import os
import scipy.special as scp
import matplotlib.pyplot as plt
if __name__ == '__main__':
    class Rcs:
        def __init__(self,r,fmin,fmax):
            self.r=r
            self.fmin=fmin
            self.fmax=fmax
            self.mass_f=[]
            self.mass_rcs=[]
        def calc(self):
            f=fmin
            def hn(n,x):
                return complex(scp.spherical_jn(n,x),scp.spherical_yn(n,x))
            while f<=fmax:
                l=3*10**8/f
                k=2*math.pi/l
                s=0
                for n in range(1,20): 
                   an=scp.spherical_jn(n,k*r)/hn(n,k*r)
                   bn=(k*r*scp.spherical_jn(n-1,k*r)-n*scp.spherical_jn(n,k*r))/(k*r*hn(n-1,k*r)-n*hn(n,k*r))
                   s+=((-1)**n)*(n+0.5)*(bn-an)
                self.mass_f.append(f)
                self.mass_rcs.append((l**2)*((abs(s))**2)/math.pi)
                f+=1000000
            plt.plot(self.mass_f,self.mass_rcs)
            plt.show()
    class Output:
        def __init__(self,mass_f,mass_rcs):
            self.mass_f=mass_f
            self.mass_rcs=mass_rcs
        def output(self):
            #создание директории
            if not os.path.isdir("results"):
                os.mkdir("results")
            os.chdir("results")
            f=open('res2.xml','w')
            f.write('<data>\n\t')
            f.write('<frequencydata>\n\t\t')
            for i in self.mass_f:
                f.write('<f>{:e}</f>\n\t\t'.format(i))
            f.write('</frequencydata>\n\t')
            f.write('<lambdadata>\n\t\t')
            for i in self.mass_f:
                l=3*10**8/i
                f.write('<lambda>{:e}</lambda>\n\t\t'.format(l))
            f.write('</lambdadata>\n\t')
            f.write('<rcsdata>\n\t\t')
            for i in self.mass_rcs:
                f.write('<rcs>{:e}</rcs>\n\t\t'.format(i))
            f.write('</rcsdata>\n\t')
            f.write('</data>\n')
            f.close()
    #Загрузка входного файла
    url='https://jenyay.net/uploads/Student/Modelling/task_rcs_02.txt'
    urllib.request.urlretrieve(url,"input.txt")
    #извлекаем информацию из файла
    F=open("input.txt","r")
    for i in F:
         if i.split()[0]=="8":
            string=i.split()
            break
    D=float(string[1])
    fmin=float(string[2])
    fmax=float(string[3])
    r=D/2

    sph=Rcs(r,fmin,fmax)
    sph.calc()
    o=Output(sph.mass_f,sph.mass_rcs)
    o.output()
    F.close()
    



