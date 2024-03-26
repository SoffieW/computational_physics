import math
import scipy.stats 
import matplotlib.pyplot as plt
import numpy as np


#-------ANALYTIC INVERSION AND REJECT ACCEPT METHOD FOR PDF OF SINE----

def inverseCDF(x): #define inverse of cumulative distribution function
    return math.acos(1-2*x)

def inversionmethod(N): #use the inversion method to plot prob dist func of sinx
    xgen=np.random.uniform(0,1,N) #generate random numbers in a uniform distribution from 0 to 1
    x_req=[] #initialize x_req
    for i in range(len(xgen)): 
        xreq=inverseCDF(xgen[i]) #compute the required x value for the ith random number contained in xgen using inverse of CDF
        x_req.append(xreq) #append the value of xreq to x_req, the empty list
    x_req=np.array(x_req) #convert list x_req to a numpy array to allow for the .size attribute later on
    return x_req #returns the array x_req containing all the x values transformed from uniform distribution function to pdf of sin(x)

def desiredPDF(x): #returns sin(x) 
    return np.sin(x)

def rejaccmethod(N): #use the reject-accept method to plot prob dist func of sinx
    xgen=np.random.uniform(0,np.pi,N) #generate random numbers in a uniform distribution from 0 to pi
    ygen=np.random.uniform(0,1,N) #generate random numbers in uniform distribution from 0 to 1
    newx=[] #initialize empty list for required x values
    for i in range(len(xgen)): 
        if ygen[i]>desiredPDF(xgen[i]):
            continue # reject required x if the generated y value is greater than the PDF function of the corresponding generated x value
        else:
            newx.append(xgen[i]) #accept the required x value if the generated y value is less than the func of the corresponding generated x value
        
    newx=np.array(newx) #convert list newx to numpy array to allow for .size attribute later on
    return newx #returns the array newx containing all the x values transformed from uniform distribution function to pdf of sin(x)

def sin():  #returns sin(x) from 0 to pi using built-in function, for comparison
    x = np.arange(0,np.pi,0.01)
    return x, np.sin(x) 

#---------PLOT THE PDF USING ANALYTIC INVERSION METHOD---------

def plotanalytic():
    sine=sin()
    N=1*10**5 #number of samples
    x_reqnew=inversionmethod(N) #generates N random values from the Inversion method
    weight = np.empty_like(x_reqnew) 
    weight.fill(2* 50 / (np.pi)/x_reqnew.size) #normalizes histogram to have area 2 (the integral of sin(x) from 0 to pi is 2)
    count, bins, ignored = plt.hist(x_reqnew, 50,edgecolor='black', weights=weight) #plots histogram for N values generated from the Inversion method 
    plt.plot(sine[0],sine[1],'r') #plots sine line graph for comparison of PDF to true function
    plt.title('Analytic Inversion Method')
    plt.xlabel('θ')
    plt.ylabel('sin(θ)')
    plt.show()
    print('Mean is',np.mean(x_reqnew)) #prints the mean of the PDF
    print('Variance is',scipy.stats.tvar(x_reqnew)) #prints variance of PDF
    print('Skew is',scipy.stats.skew(x_reqnew)) #prints skew of PDF

    return
#---------PLOT THE PDF USING REJECT-ACCEPT METHOD---------

def plotrejacc(): 
    sine=sin() 
    N=1*10**5 #number of samples
    newx=rejaccmethod(N)   
    weight = np.empty_like(newx)
    weight.fill(2* 50 / (np.pi)/newx.size)
    count, bins, ignored = plt.hist(newx, 50,edgecolor='black', weights=weight)
    plt.plot(bins, np.ones_like(bins), linewidth=0)
    plt.plot(sine[0],sine[1],'r')
    plt.title('Reject-Accept Method')
    plt.xlabel('θ')
    plt.ylabel('sin(θ)')
    plt.show()
    print('Mean is',np.mean(newx))
    print('Variance is',scipy.stats.tvar(newx))
    print('Skew is',scipy.stats.skew(newx))
    return

#--GENERATE X AND Y POSITIONS FOR PARTICLES HITTING DETECTOR SCREEN--

def decayhitmap():
    N=5*10**6 #number of decays
    tau=550*10**-6 #mean lifetime before decay
    vz=2000 #speed of particles in beam
    distance=tau*vz #mean distance = speed x mean lifetime
    theta=inversionmethod(N) #generate theta angles using inversion method, for angles perpendicular to detector plane for N particles
    phi = 2*np.pi*np.random.random(N) #generate phi angles from 0 to 2pi, for angles in the detector plane for N particles
    theta -= np.pi/2 
    dist=[] 
    for i in range(0,N):
        s=np.random.exponential(scale=distance) #generate random distance from an exponential distribution
        dist.append(s) #append distance values to array
    x=[]
    y=[]
    rx=[]
    ry=[]
    for i in range(0,N): 
        if 2-dist[i]>0: #distance values generated beyond the detector screen will not be considered
            xpos=(2-dist[i])*np.tan(phi[i]) #find the x position using trig
            ypos= (2-dist[i])*np.tan(theta[i]) / np.cos(phi[i]) #find the y position using trig
            x.append(xpos) 
            y.append(ypos)
            res=resolution(xpos,ypos) #use the resolution function to produce a smearing of each x and y position
            rx.append(res[0])
            ry.append(res[1])
    histrange=[-2,2] # defines size of the detector screen
    bins=100 #defines number of bins
    h=plt.figure(figsize=(7.4, 6 )) #defines the dimensions of the figure
    h=plt.hist2d(rx,ry,bins, range=[histrange, histrange],normed=True) #plots a normalized 2D histogram of the smeared x and y positions
    plt.xlabel('position (m)') 
    plt.ylabel('position (m)')
    plt.gca().axis("on")
    plt.colorbar(h[3]) #plots the colorbar of probability
    plt.show()
    slicehistx=[] 
    slicehisty=[]
    halfbinwidth=(abs(histrange[1]-histrange[0])/bins)*0.5 
    for j in range(len(rx)): 
        if ry[j]<halfbinwidth and ry[j]>-halfbinwidth: #find the values along the y=0 axis
                slicehistx.append(rx[j]) #append the x values along the y=0 axis to an array 
    plt.title('No. of counts along y=0')
    plt.hist(slicehistx, bins=100, range=[-2, 2],edgecolor='black', linewidth=0.5,color='r') #plot a 1d histogram which is a slice of the 2d hist along the y=0 axis
    plt.ylabel('counts')
    plt.xlabel('position (m)')
    plt.show()
    print('Standard deviation is',np.std(slicehistx)) #prints the spread of the overall decay diffraction pattern in the x direction
    for k in range(len(rx)):
        if rx[k]<halfbinwidth and rx[k]>-halfbinwidth: #find the values along the x=0 axis
                slicehisty.append(ry[k]) #append the y values along the x=0 axis to an array
    plt.title('No. of counts along x=0')
    plt.hist(slicehisty, bins=100, range=[-2, 2], edgecolor='black',linewidth=0.5,color='r') #plot a 1d histogram which is a slice of the 2d hist along the x=0 axis
    plt.ylabel('counts')
    plt.xlabel('position (m)')
    plt.show()
    print('Standard Deviation is', np.std(slicehisty)) #prints the spread of the overall decay diffraction pattern in the y direction
    return 

def resolution(xvalue,yvalue): 
    x_smeared=xvalue+np.random.normal(loc=xvalue,scale=0.1) #Gaussian smear each x position with a resolution of 10cm
    y_smeared=yvalue+np.random.normal(loc=yvalue,scale=0.3) #Gaussian smear each y position with a resolution of 30cm
    return x_smeared, y_smeared #return the smeared values




#-------FIND X-SECTION LIMIT FOR COLLIDER EXPERIMENT-------

def xsection(uncertainty):  
    N=1000  #number of pseudoexperiments to be generated
    L=12    #Luminosity in /nb (10^-37 m^2)
    values=[] 
    bg=np.zeros(N) 
    signal=np.zeros(N)
    L=np.zeros(N) 
    for sigma in np.arange(0,1,0.001): #range of sigma to be tested
        values=[]
        for i in range(N): 
            L[i]=np.random.normal(loc=12, scale=uncertainty)
            bg[i]=np.random.normal(loc=5.7,scale=0.4) #generate gaussian distributed background values
            if L[i]>0:
                signal[i]=np.random.poisson(lam=L[i]*sigma+bg[i]) #generate signal value by adding the gaussian distributed bg value to L*sigma for lamda
                if signal[i]>5: #if the signal value > observed value, it satisfies the detection of the candidate particle X
                    values.append(signal[i])
                else: 
                    continue
        numvalues=len(values)/N #the fraction of values satisfying the detection criterion over the total number of pseudoexperiments, i.e.: the confidence level
        if numvalues==0.950: #if the confidence level =0.95, return the x-section value
            return sigma


def xsectionparams(uncertainty): #generate the average and standard deviation of several runs of the xsection function
    sigmas=[]
    for i in range(8): 
        sig=xsection(uncertainty)
        if sig!= None: #condition to remove all None values returned by xsection(), which may happen if the random variables generated in xsection() do not satisfy the detection criterion
            sigmas.append(sig)
    sigmas=np.array(sigmas,dtype='float64') #converts dtype to float64
    average= np.mean(sigmas,dtype='float64') #calculate the average of all sigma values
    sd=np.array(sigmas).std(dtype='float64') #calculate the standard deviation of all sigma values
    return average, sd 

def plotxsecthist(uncertainty): #plots a histogram using the returned values from xsectionparams()
    x_section=xsectionparams(uncertainty)
    xsects=np.random.normal(loc=x_section[0],scale=x_section[1],size=2000) #generates a number of random samples from the distribution with mean and scale defined in xsectionparams()
    bins=math.ceil(math.sqrt(len(xsects))) #calculates the no. of bins for the histogram
    plt.hist(xsects,bins,edgecolor='black') #plots histogram of the sigma values 
    plt.xlabel('Limit on x-section (nb)')
    plt.ylabel('No.of counts')
    plt.show()
    print('Average limit on the x-section is', round(x_section[0],2),'±', round(x_section[1],2),'nb')   
    return

MyInput = '0'

while MyInput != 'q':
    MyInput = input('Enter from the following choices:\n "1" to see Monte Carlo Methods \n "2" to see plot of radioactive decay\n "3" to find limit on x-section \n "q" to quit ')

    if MyInput == "1":
        print('You have chosen to see the PDF of sin(θ) using the Analytic Inversion method and the Reject-Accept method...')
        print(plotanalytic())
        print(plotrejacc())

    elif MyInput == "2":
        print('You have chosen to see the hitmap of radioactive decay for a beam of unstable nuclei...')
        print(decayhitmap())
     
    elif MyInput == "3":
        print('You have chosen to find the limit of the production x-section for the collider experiment...')  
        uncertainty=6
        while uncertainty>5:
            myinput=input('Enter a value for the uncertainty on luminosity between 0 and 5:')        
            uncertainty=float(myinput)
        print(plotxsecthist(uncertainty))

    elif MyInput != "q":
        print('This is not a valid choice')
print('You have chosen to finish - goodbye.')
    

