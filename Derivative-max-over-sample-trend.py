import numpy as np
import matplotlib.pyplot as plot
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

sample=['5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _24.6 .txt','5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _29.9 .txt','5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _35 .txt'
,'5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _40 .txt'
,'5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _45 .txt'
,'5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _50.0999980.txt'
,'5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _55.099998 .txt'
,'5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _60.099998 .txt'
,'5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _65.099998 .txt'
,'5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _70.199997 .txt'
,'5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _75.400002 .txt'
,'5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _80.400002 .txt'
,'5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _84.699997 .txt'
,'5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _80.4000020.txt'
,'5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _75.099998 .txt'
,'5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _70.099998 .txt'
,'5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _65 .txt'
,'5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _60 .txt'
,'5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _55 .txt'
,'5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _50.099998 .txt'
,'5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _45.0 .txt'
,'5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _40.099998 .txt'
,'5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _35.0 .txt'
,'5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _30.1 .txt'
,'5 mM dCG9, 4M NaCl, thermal _25C - 85C, 5C step, go back _0 _25.1 .txt'
]
base='4M NaCl _baseline for thermal scan _0.0001 _30.1 .txt'
def fun(x,a,b):
       return(a*x+b)


n=len(sample)
lamda=np.empty(152)
AC_DC=np.empty([152,n])
baselineAC_DC=np.empty(152)
HT=np.empty([152,n])
HT_reals=np.empty(n)



baseline=np.loadtxt(base,skiprows=41)
count=0
for i in range(len(sample)):
      first=True
      for k in range(152):
            if count==0:
                lamda[k]=baseline[k][0]
                baselineAC_DC[k]=baseline[k][1]  
            current_sample=np.loadtxt(sample[i],skiprows=41)
            AC_DC[k,i]=(current_sample[k][1]-baselineAC_DC[k])
            HT[k,i]=(current_sample[k][2])
            if k==151 and first==True:
                HT_reals[i]=135
            if first==True and HT[k,i]<-400 and k!=0:
                if k<135:
                       HT_reals[i]=135
                else:
                       HT_reals=k
                first=False
      count+=1

print(HT_reals)
der_pool=np.empty([n,152])
for i in range(len(sample)):
      der=savgol_filter(np.transpose(AC_DC)[i],11,1,deriv=1,delta=-1)
      der_pool[i]=der


all_turning_points=[]
count=1
for j in range(len(sample)):
      start=HT_reals[j]
      value=der_pool[j]
      turning_points=[]
      for k in range(10,int(start)):
            if k+2>start:
                pass
            elif np.sign(value[k])!=np.sign(value[k+1]):
                   if np.sign(value[k])==np.sign(value[k-1])==np.sign(value[k-2]) and np.sign(value[k+1])==np.sign(value[k+2]):
                           turning_points.append(k) 
                   elif np.sign(value[k])==np.sign(value[k-1])==np.sign(value[k-2])==np.sign(value[k-3]) and np.sign(value[k+1])==np.sign(value[k+3]):
                            turning_points.append(k)

      all_turning_points.append(turning_points)
     
      

'''
for i in range(len(sample)):
      f = interp1d(lamda, der_pool[i], bounds_error=False, kind='cubic')
      plot.figure()
      plot.plot(lamda,np.transpose(AC_DC)[i],linestyle='-')
      plot.plot(lamda,der_pool[i],linestyle='--')
      plot.plot(lamda,f(lamda))
      plot.plot([170,320],[0,0],'k-')
      x=[]
      y=[]
      for j in all_turning_points[i]:
            x.append(np.average([lamda[j],lamda[j+1]]))
            y.append(np.average([np.transpose(AC_DC)[i][j],np.transpose(AC_DC)[i][j+1]]))
      plot.plot(x,y,'kx')

            
      plot.show()
            
 '''
shift=np.zeros([len(all_turning_points[0]),len(all_turning_points)])
for l in range(len(all_turning_points[0])):
      for p in range(len(all_turning_points)):
            shift[l,p]=np.average([lamda[int(all_turning_points[p][l])],lamda[int(all_turning_points[p][l])+1]])



Temps=[25,30,35,40,45,50,55,60,65,70,75,80,85,80,75,70,65,60,55,50,45,40,35,30,25]
 


for i in range(len(sample)):
      popt,popc=curve_fit(fun,Temps,shift[i])
      x=np.linspace(20,90,10)
      y=fun(x,popt[0],popt[1])
      plot.figure()
      plot.plot(x,y)
      plot.plot(Temps,shift[i],'kx')
      plot.savefig('peakshift_temps'+str(i)+'.png')
      plot.show()
                  
            










