import math
import cmath
import numpy as np
import matplotlib.pyplot as plt

def simpson(x1,x2,f,N): #function which performs simple 1d simspon's rule integration taking as arguements lower limit, upper limit, function to be integrated and number of intervals respectively.
    
    h=(x2-x1)/N #h is the width of  a strip defined by the difference between the limits divided by number of intervals
    y=0.0
    x=x1 #set starting value of x to lower limit
    
    for i in np.arange(1, N/2 +1): #summing odd order y terms
    
        y+=4*f(x)
        x+=2*h #increment x by 2*strip width to ensure only the odd order terms are multiplied by 4 and summed.
    
    x=x1+2*h 
    for i in np.arange(0, N/2): #summing even order y terms
    
        y+=2*f(x)
        x+=2*h
    
    integral= (h/3)*(y+f(x1)+f(x2))    
    
    return integral



def expfuncXB(x,xp): #gives the function to be integrated for the 1d fresnel diffraction pattern in the x direction for part b
    
    j=cmath.sqrt(-1)
    g=(k/2*z)*((x-xp)**2)    
    
    return cmath.cos(g)+cmath.sin(g)*j #used the trig equivalent of the exponential function as the exponential did not seem to compute correctly
    


def simpsonX(xp1,xp2,x,f,N): #integrates 1d fresnel diffraction pattern in the x direction for part b
    
    h=(xp2-xp1)/N
    y=0
    xp=xp1
    for i in np.arange(1, N/2 +1): #summing odd order y terms
    
        y+=4*f(x,xp)
        xp+=2*h
    
    xp=xp1+2*h
    for i in np.arange(2,N/2): #summing even order y terms
    
        y+=2*f(x,xp)
        xp+=2*h
    
    integral= (h/3)*(y+f(x, xp1)+f(x, xp2))    
    
    return integral

def expfuncX(x,xp): #gives the x function to be integrated along the x direction for part c


    
    return cmath.exp(((1j*k)/(2*z))*((x-xp)**2))


def X(xp1,xp2,x,f,N): #integrates the x function expfuncX for part c
    
    h=(xp2-xp1)/N
    ff=0
    xp=xp1
    for i in np.arange(1, N/2 +1): #summing odd order func terms
    
        ff+=4*f(x,xp)
        xp+=2*h
    
    xp=xp1+2*h
    for i in np.arange(2,N/2): #summing even order func terms
    
        ff+=2*f(x,xp)
        xp+=2*h
    
    integral= (h/3)*(ff+f(x, xp1)+f(x, xp2))    
    
    return integral

 
    

def expfuncXY(x,y,yp):  #gives the 2d func to be integrated
    

    
    return X(xp1,xp2,x,expfuncX,N)*(cmath.exp(((1j*k)/(2*z))*((y-yp)**2)))    #multiplies the integration along x by the exp function for y to give the function to be integrated in 2d


def simpsonXY(x,yp1,yp2,y,f,N): #integrates 2d function
    
    h=(yp2-yp1)/N
    ff=0
    yp=yp1
    for i in np.arange(1, N/2 +1): #summing odd order func terms
    
        ff+=4*f(x,y,yp)
        yp+=2*h
    
    yp=yp1+2*h
    for i in np.arange(2,N/2): #summing even order func terms
    
        ff+=2*f(x,y,yp)
        yp+=2*h
    
    integral= ((E0*k)/(2*(math.pi)*z))*(h/3)*(ff+f(x,y, yp1)+f(x, y, yp2))    
    
    return integral

def expfunccircle(x,y,yp):  #gives the 2d func to be integrated with limits for a circle
    xp1=-math.sqrt((yp1)**2)-(yp**2)
    xp2=math.sqrt(((yp2)**2)-(yp**2))

    
    return X(xp1,xp2,x,expfuncX,N)*(cmath.exp(((1j*k)/(2*z))*((y-yp)**2)))

MyInput = '0'

while MyInput != 'q':
    MyInput = input('Enter a choice, "a", "b", "c", "d" or "q" to quit: ')


    if MyInput == 'a':
        print('You have chosen part (a)')

        N=1
        while N%2!=0 and N!=int():
            Ninput=input('Enter an even integer value of N:')
            N=int(Ninput)
        
        print('The integral of sin between 0 and pi is:', simpson(0,math.pi,math.sin,N))  


    elif MyInput == 'b':
        print('You have chosen part (b)')
        
        lamda=0
        while lamda<40 or lamda>80:
            lamdainput=input('Enter a wavelength between 40 and 80 nm:')
            
            lamda=(float(lamdainput))
            
            
        z=2
        while z<3 or z>8:
            
            zinput=input('Enter the aperture-screen distance between 3 and 8 cm:')
            z=(float(zinput))
            
        N=1
        while N%2!=0 and N!=int():
            Ninput=input('Enter an even integer value of N:')
            N=int(Ninput)
            
        k=(2*math.pi)/(lamda*1e-9) #k is the wave number
        z=z*0.01 #z is the distance between aperture and screen
        NumPoints = 200
        xmin = -1 
        xmax =1 #xmin and xmax define screen width
        dx = (xmax - xmin) / (NumPoints - 1)
        xvals = [0.0] * NumPoints
        yvals = np.zeros(NumPoints)
        for i in range(NumPoints):
            xvals[i] = xmin + i * dx
            yvals[i] = abs(simpsonX(-1e-6,1e-6,xvals[i],expfuncXB,N)**2) #intensity of the pattern is proportional to mod of integral squared
        plt.plot(xvals,yvals)
        plt.show()

    elif MyInput == 'c':
        print('You have chosen part (c)')
        lamda=0.0000006
        k=(2*math.pi)/lamda
        
        h=6.67e-34
        c=3e8
        e0=8.85e12
        E0=h*c/lamda
        xp1=-1e-3 
        xp2=1e-3
        N=50
        y=-1
        yp1=-1e-3
        yp2=1e-3
        x=-1
        zinput=input('Enter a value of z:') #allows user to enter a value for aperture-screen distance
        z=float(zinput)
      
        NumPoints = 200
        xmin = -0.0014
        xmax =0.0014
        ymin=-0.0014
        ymax=0.0014
        delta = (xmax - xmin) / (NumPoints - 1)
        delta1=(ymax-ymin)/(NumPoints-1)
        xvals = [0.0] * NumPoints
        yvals = [0.0] * NumPoints
        intensity = np.zeros( (NumPoints,NumPoints) )
        for j in range(NumPoints):
            yvals[j] =ymin+j * delta1
            for i in range(NumPoints):
                xvals[i] = xmin+i * delta
                intensity[i,j] =e0*c*((abs(simpsonXY(xvals[i],-1e-3,1e-3,yvals[j],expfuncXY,50))**2))
        plt.imshow(intensity)
        plt.show()


        
    elif MyInput =='d':
        print('You have chosen part (d)')
        lamda=0.0000006
        k=(2*math.pi)/lamda
        
        h=6.67e-34
        c=3e8
        e0=8.85e12
        E0=h*c/lamda
        N=50
        y=-1
        yp1=-1e-3
        yp2=1e-3
        x=-1
        zinput=input('Enter a value of z:') 
        z=float(zinput)
        yp1=-1e-3
        yp2=1e-3
        NumPoints = 200
        xmin = -0.0014
        xmax =0.0014
        ymin=-0.0014
        ymax=0.0014

        delta = (xmax - xmin) / (NumPoints - 1)
        delta1=(ymax-ymin)/(NumPoints-1)
        xvals = [0.0] * NumPoints
        yvals = [0.0] * NumPoints
        intensity = np.zeros( (NumPoints,NumPoints) )
        for j in range(NumPoints):
            yvals[j] =ymin+j * delta1
            for i in range(NumPoints):
                xvals[i] = xmin+i * delta
                intensity[i,j] =e0*c*((abs(simpsonXY(xvals[i],yp1,yp2,yvals[j],expfunccircle,50))**2))
        plt.imshow(intensity)
        plt.show()

    elif MyInput != 'q':
        print('This is not a valid choice')
print('You have chosen to finish - goodbye.')
    