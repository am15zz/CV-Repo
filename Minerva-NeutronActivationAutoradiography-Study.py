#############################
# This is a python 3 script #
#############################

########
# Imports 
#####

import math
import numpy as np 
import matplotlib.pyplot as plot
from scipy.integrate import simps
from  PIL import Image
from scipy.optimize import nnls
import scipy
from decimal import Decimal as D



###
# Create pigment composition
##########



######## 
# Create the Colour palette
########


# define the path to the colours on Rembrandt's palette and place them in an iterable list 


palette_names=['Yellow-Ochre.jpeg','Azurite.jpeg','Umber.jpeg','Lead-White.jpeg','Malachite.jpeg','Madder-Lake.jpeg']

palette_list=[]    # creates a list of lists of RGB values for each colour swatch with a 4th check parameter for the linear combination 
counter=[]         # establishes a list of a list of pairs for each colour swatch [name,occurrences]  

# swatch_colour takes in the path to an image file (string) and outputs a list of the average RBG value on the image with a check parameter for the linear combination.

#swatch_colour: str->list(float,float,float,float)

def swatch_colour(image): # take in a string with the file name and extension can be in gif jpeg or png format
    rgb_image=Image.open(image).convert('RGB') 
    arr=np.asarray(rgb_image) 
    colour=image.split('.')[0]
    r=0
    g=0
    b=0
    n_pixels=rgb_image.size[0]*rgb_image.size[1] # calculate the size of the image array
    rgb_image.close()
    for i in arr:    
        for j in i:  # sum over all RGB values in the array
            r+=j[0]
            g+=j[1]
            b+=j[2]
    r_bar=r/n_pixels # calculate the average value for each colour
    g_bar=g/n_pixels
    b_bar=b/n_pixels
    return([r_bar,g_bar,b_bar,0.01]) # return the average RGB colours with a check variable 0.01 -helps set limit for nnls which doesn't take many external factors. 0.01 should eliminate all trivial solutions for black.


for swatch in palette_names:    # interate swatch_colour over all swatches in Fembrandt's palette
    counter.append([swatch.split('.')[0],0])
    palette_list.append(swatch_colour(swatch))

x=np.array(palette_list)  # create an array of the RGB colour values on the palette


colour_map=[]  # create list to update with the linear combination values for the swatch colours 
FILENAME= 'minerva.jpeg' #image can be in gif jpeg or png format  # read in the 4K minerva Image 
image=Image.open(FILENAME).convert('RGB')
arr=np.asarray(image)
image.close()

# iterate over all pixels in the minerva image and do a linear combination to find the most accurate combination of colours to recreate that RGB value using only the RGB values from the palette

for i in arr :
    colour_list=[]
    coefficients=[]
    for j in i :
        coefficients=nnls(x.T,np.array((j.tolist())+[0.01])) # make the RGB minerva value the same length as the others with the same 4th fitting parameter 
        colour_list.append(coefficients[0].tolist()) # add the coefficients to a list 
    colour_map.append(colour_list) update the list to get the correct dimension array (needed to recreate the image hurts memory but allows for the recreation)

colour_profile=np.array(colour_map)
arr=0
colour_list=0

# iterate over the values from the linear combination RGB array to find the percent of each pixel belongs to what colour. 

for i in colour_profile: 
    for j in i : 
        k=0 
        summ=sum(j)
        for p in j:  
            (counter[k])[1]+=p/summ # add the percent composition of each pixel to the counter for that colour 
             k+=1
                 
width,height=image.size
pixel_count=width*height
   
for i in colour_profile:         #iterate over the values of percent pixel composition summs to add a value of 1 per pixel to lead white to account for the layer of lead white paint underneath.
    if i[0]=='Lead-White': 
         i[1]+=pixel_count         
                
####
# Convert the linear combination RGB colour array back to RGB colour array
####

for i in counter:     # print out the percent values for each colour 
    print(i[0]+":"+str(i[1]*100/*2))+"%") # pixel count is multiplied by 2 becaue of the excess added to lead white to account for the under layer 
counter=0
newpicinset=[]
newpic=[]

####
# Recreate the Image 
####

for i in colour_profile: # iterate over the linear combination array
    newpicinset=[]
    for j in i: 
        k=0 # keeps track of palette placement
        r_value=0
        g_value=0
        b_value=0
        for p in j: #p is the coefficient of an element in the palette 
            r_value+=(palette_list[k])[0]*(p)
            g_value+=(palette_list[k])[1]*(p)
            b_value+=(palette_list[k])[2]*(p)
            k+=1
            if r_value > 255:  # cases when the maximized linear combination goes over the max RGB value. 
               r_value=255
            if g_value>255: 
               g_value=255
            if b_value>255: 
               b_value=255
        newpicinset.append([r_value,g_value,b_value]) # create new RGB array
    newpic.append(newpicinset)
newpicinset=-0
picture=np.array(newpic,dtype=np.uint8) # set up for image writing 

####
# Save the Digital recreation of the image 
####

scipy.misc.imsave('try2.png', picture)





####
#Create cross section plots
########

# list of all of the files

lofiles=['1H.txt','2H.txt','12C.txt','13C.txt','16O.txt','17O.txt','18O.txt','27Al.txt','55Mn.txt',
'54Fe.txt','56Fe.txt','57Fe.txt','58Fe.txt','63Cu.txt','65Cu.txt',
'204Pb.txt','206Pb.txt','207Pb.txt','208Pb.txt']

for i in lofiles:  # read all files and plot the cross sections as a function of energy in ev
    file1=open(i,'r')
    E=[]
    sigma=[]
    lines=file1.readlines()[:-1]
    for line in lines:
        energy,cross=line.split(' ')
        E.append((float(energy)))
        sigma.append((float(cross)))
    plot.figure()
    plot.title('Cross section for '+i[:-4])
    plot.xlabel('Energy (eV)')
    plot.ylabel('Cross section (b)')
    plot.loglog(E,sigma,'g-')
    plot.savefig(i[:-4]+'.png')
    plot.show()
    
###
# Numeric integration over cross sections
#########    


########
# Constants
#####

kb=1.38*10**(-23)
mn=1.674929*10**(-27)

#
#velocitymaker: float-> float converts the electron volt energy into velocities
#

def velocitymaker(energy):
    v=((energy*(2.0)*(1.6022*10**(-19)))/((mn)))**(1.0/2.0)
    return(v)

#
# barntom2: float->float  converts the cross section from barns to m^2
#

def barntom2(area):
    m=(area*(10.0**(-28.0)))
    return(m)

#
# listmaker str->listoflists : opens a given file within the directory and splits the energy and crosssection compoenets and converts them to velocity and crosssection in m^2
#

def listmaker(filename):
      file1=open(filename,'r')
      v=[]
      sigma_v=[]
      lines=file1.readlines()[:-1]
      for line in lines:
          energy,cross=line.split(' ')
          v.append(velocitymaker(float(energy)))
          sigma_v.append(barntom2(float(cross)))
      return([v,sigma_v])

#
# maxwellboltsmann: float+float->float takes in a velocity and a temperature and spits out the value of the maxwell boltzmann distribution at that point
#

def maxwellboltsmann(v,T):
    f=(((4.0)*(np.pi)*v**(2.0))*((mn)/(4.0*np.pi*kb*T))**(3.0/2.0))*(np.exp(((-1)*(mn)*v**(2.0))/(2.0*kb*T)))
    return(f) 
    
inputs=input("Input file name and extension: ") # user input for the file name and extensions

dataset=listmaker(inputs) 
velocity=np.array(dataset[0])
cros=np.array(dataset[1])
d=velocity[0]

y1=np.array([])
counter=0

for i in velocity: 
    y1=np.append(y1,float(i*maxwellboltsmann(i,20)*(cros[counter])*1e9*1.10062594e22))  # calculate the value of the numeric integral for each element in velocity
    counter+=1

area=simps(y1,velocity) # calculate the simpson integral

# plot the integration functions. 
plot.figure()
plot.loglog(velocity,y1,'b-')
plot.title(inputs[:-4]+'something')
plot.xlabel('velocity')
plot.ylabel('something')
plot.annotate('Area Under the Curve:\n'+ str(area), xy=(0.9, .5),xycoords='axes fraction',horizontalalignment='right', verticalalignment='bottom')
plot.savefig(inputs[:-4]+'.png')
plot.show()



###
# Activitiy and Exposure
########

    
    
# t changes with each simulation specifying the time region you want to look at in the activity graphs

t=np.linspace(0,60,100)

# list of all the elements for naming purposes

list_of_elements=[r'$^{3}H$',r'$^{14}C$',r'$^{19}O$',r'$^{28}Al$',r'$^{56}Mn$',r'$^{55}Fe$',r'$^{59}Fe$',r'$^{64}Cu$',r'$^{66}Cu$',r'$^{205}Pb$',r'$^{209}Pb$']

# list of lambda valuues in years

list_of_lambdasy=[0.056262,1.21604e-4,826559.1545,3.710965,2354.488137,0.252604657,5.686003372,478.0704542,1.624562054,4e-8,1866.575415]

list_of_lambdass=[]

for i in list_of_lambdasy: # convert lambdas from years to seconds (so the activities will be in Bq
    list_of_lambdass.append(i/31557600)

list_of_N0s=[4.31759438807,11.6999503588,1.32879141049,1997.80350455,113586.607053,
19265.3907216,9875.45564172,38247.844942,18385.7303831,5653.55870412,1.98931256685] # a list of the initial number of each element per cm^2

list_of_gammas=[1,1,0.4838,0.8372,0.8580,1,0.624,0.1087,0.0458,0.5355,1]   
# list of the gamma factors (exposure rate constants) 
# if an exposure rate constant was not availible the number 1 was used

time=43200     # total time the painting is being bombarded with neutrons

for i in range(len(list_of_N0s)):
    list_of_N0s[i]=list_of_N0s[i]*time

activated=sum(list_of_N0s)         # calculate the total ammount of altered particles

print(activated)

Act=[]

#
# l1activity float+float+float->float determines the activity of a given quantity of a given isotope at a time t
#

def l1activity(lambdas,No,t):               
    A=lambdas*No/3.7e10*np.exp(-1*lambdas*t)
    return(A)

# loop over all unstable isotopes present after neutron bombardment and calculate activity over different temperature ranges
c=0
for i in range(0,len(list_of_lambdass)):
    As=[]
    for j in t:
        As.append(l1activity(list_of_lambdass[i],list_of_N0s[i],j))
    if As[-1]>=1e-80:
       Act.append(As)
       c+=1
    else: 
       del list_of_elements[c]
NUM_COLORS = c
cm = plot.get_cmap('gist_rainbow')
fig = plot.figure()
ax = fig.add_subplot(111)
ax.set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
for i in range(NUM_COLORS):
    ax.loglog(t,Act[i],label=r'Activity'+list_of_elements[i])
ax.legend(loc=0)
plot.title('Activity of Nucleides')
plot.xlabel('time in seconds')
plot.ylabel('Activity in Curie')
plot.savefig('act1.png')
plot.show()    
	
# 
# lactivity float+float+float->Decimal same as l1activity but returns a more accurate decimal
#

def lactivity(lambdas,No,t):
    A=D(lambdas)*D(No)*D.exp(D(-1*D(lambdas)*D(t)))
    return(A)
#################
# exposure
#################

times=[86400,604800,2592000,2592000*12*50]

overallexposure1d=0
overallexposure1w=0
overallexposure1m=0
overallexposure1y=0

for i in range(len(list_of_lambdass)):
    for j in times: 
        if j==86400:
           overallexposure1d+=lactivity(list_of_lambdass[i],list_of_N0s[i],j)*D(list_of_gammas[i])
        elif j==604800:
           overallexposure1w+=lactivity(list_of_lambdass[i],list_of_N0s[i],j)*D(list_of_gammas[i])
        elif j==2592000:
         overallexposure1m+=lactivity(list_of_lambdass[i],list_of_N0s[i],j)*D(list_of_gammas[i])
        else:
         overallexposure1y+=lactivity(list_of_lambdass[i],list_of_N0s[i],j)*D(list_of_gammas[i])

# these values are the exposure rate in terms of R*m^2/hr

print(overallexposure1d)
print(overallexposure1w)
print(overallexposure1m)
print(overallexposure1y)

# adjust those values to give Sv for a person standing 2m away

overallexposure1d=overallexposure1d*D(0.00877/(2)**2)
overallexposure1w=overallexposure1w*D(0.00877/(2)**2)
overallexposure1m=overallexposure1m*D(0.00877/(2)**2)
overallexposure1y=overallexposure1y*D(0.00877/(2)**2)
print(overallexposure1d)
print(overallexposure1w)
print(overallexposure1m)
print(overallexposure1y)

    
    
    
    
    
    
    
    
    
    
    
    
