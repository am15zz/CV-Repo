import numpy as np
import matplotlib.pyplot as plot
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d


sample=['Brock cell _CSA 0 ticks_0_25.1.txt',
'Brock cell _CSA 100 ticks_0_25.txt',
'Brock cell _CSA 10 ticks_0_25.1.txt',
'Brock cell _CSA 20 ticks_0_25.txt',
'Brock cell _CSA 30 ticks_0_25.1.txt',
'Brock cell _CSA 40 ticks_0_25.txt',
'Brock cell _CSA 50 ticks_0_25.txt',
'Brock cell _CSA 60 ticks_0_25.txt',
'Brock cell _CSA 70 ticks_0_25.txt',
'Brock cell _CSA 80 ticks_0_25.txt',
'Brock cell _CSA 90 ticks_0_25.txt'
]
base='Brock cell _CSA 50 ticks_0_25.0.txt'


def fun(x,a,b):
       return(a*x+b)



n=len(sample)
lamda=np.empty(152)
AC_DC=np.empty([152,n])
baselineAC_DC=np.empty(152)
HT=np.empty([152,n])



baseline=np.loadtxt(base,skiprows=41)
count=0
cutoffs=[]
for i in range(len(sample)):
      first=True
      for k in range(0,152):
            if count==0:
                lamda[k]=baseline[k][0]
                baselineAC_DC[k]=baseline[k][1]  
            current_sample=np.loadtxt(sample[i],skiprows=41)
            AC_DC[k,i]=(current_sample[k][1]-baselineAC_DC[k])
            HT[k,i]=current_sample[k][2]
            if first==True and HT[k,i]<-400 and k!=0:
                  cutoffs.append(k)
                  first=False


colors = plot.cm.nipy_spectral(np.linspace(0,1,n))

dyn_der_pool=[]
for i in range(len(sample)):
      der=savgol_filter(np.transpose(AC_DC)[i],11,1,deriv=1,delta=-1)
      f = interp1d(lamda, der, bounds_error=False, kind='cubic')
      dyn_der_pool.append(f(lamda))


der_pool=np.array(dyn_der_pool,dtype=object)





# current working secion

CD=np.transpose(AC_DC)
all_turning_points=[]
count=1
max_turn=[0,0]
for j in range(len(sample)):
      start=cutoffs[j]
      value=der_pool[j]
      turning_points=[]
      for k in range(20,int(start)):
            if k+2>start:
                pass
            elif np.sign(value[k])!=np.sign(value[k+1]):
                   if np.sign(value[k])==np.sign(value[k-1])==np.sign(value[k-2]) and np.sign(value[k+1])==np.sign(value[k+2]) and abs(value[k+2]-value[k-2])/5>0.01:
                        if (value[k+2]-value[k-2])<0 and CD[j][k]>0:
                            pass
                        else:
                            turning_points.append(k) 
                   elif np.sign(value[k])==np.sign(value[k-1])==np.sign(value[k-2]) and np.sign(value[k+1])==np.sign(value[k+3]) and abs(value[k+2]-value[k-2])/5>0.01:
                        if (value[k+2]-value[k-2])<0 and CD[j][k]>0:
                            pass
                        
                        else:
                            turning_points.append(k)
      if max_turn[0]<len(turning_points):
            max_turn=[len(turning_points),j]
       
      all_turning_points.append(turning_points)



for i in range(len(sample)):
      
      plot.figure()
      plot.plot(lamda[1:cutoffs[i]],np.transpose(AC_DC)[i][1:cutoffs[i]],linestyle='-')
      plot.plot(lamda[1:cutoffs[i]],der_pool[i][1:cutoffs[i]],linestyle='--')
      plot.plot([170,320],[0,0],'k-')
      x=[]
      y=[]
      for j in all_turning_points[i]:
            x.append(lamda[j+1])
            y.append(np.transpose(AC_DC)[i][j+1])
      plot.plot(x,y,'kx')
      plot.show()
      plot.close()
            
ranges=[]
for items in all_turning_points[max_turn[1]]:
    ranges.append(np.linspace(items-10,items+10,21))

temps=Temps=[0,100,10,20,30,40,50,60,70,80,90]

total_plot = [[] for _ in range(max_turn[0])]
for item in range(len(all_turning_points)):
    for entry in range(len(all_turning_points[item])): 
        for subset in range(len(ranges)):
            if all_turning_points[item][entry] in ranges[subset]:
                total_plot[entry].append([lamda[all_turning_points[item][entry]],temps[item],CD[item][all_turning_points[item][entry]]])



dyn_interp=[]
for i in range(max_turn[0]):
    dyn_interp.append(np.transpose(np.array(total_plot[i])))
count=0

for i in (dyn_interp):
      popt,popc=curve_fit(fun,i[1],i[0])
      x=np.linspace(0,100,10)
      y=fun(x,round(popt[0],3),round(popt[1],3))
      plot.figure()
      plot.plot(x,y,label='y='+str(round(popt[0],2))+'x + '+str(round(popt[1],2)))

      plot.legend()
      plot.xlabel('Ticks above minimum')
      plot.ylabel('Peak wavelength $\lambda$ (nm)')
      plot.plot(i[1],i[0],'kx')
      plot.savefig('peakshift_pathlength'+str(count)+'.png')
      plot.show()
      count+=1



dyn_interp=[]
for i in range(max_turn[0]):
    dyn_interp.append(np.transpose(np.array(total_plot[i])))
count=0
for i in (dyn_interp):
      popt,popc=curve_fit(fun,i[1],i[2])
      x=np.linspace(0,100,10)
      y=fun(x,round(popt[0],3),round(popt[1],3))
      plot.figure()
      plot.plot(x,y,label='y='+str(round(popt[0],2))+'x + '+str(round(popt[1],2)))
      plot.legend()
      plot.xlabel('Ticks above minimum')
      plot.ylabel('CD peak maximum')
      plot.plot(i[1],i[2],'kx')
      plot.savefig('peakshift_CD'+str(count)+'.png')
      plot.show()
      count+=1
