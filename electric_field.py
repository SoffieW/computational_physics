# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 12:33:40 2019

@author: Soffie
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from math import ceil
from math import floor
from scipy import linalg as lin
import pylab

def grid(Nx):
    #grid dimensions


    xmax=15 #max. position in x dirn
    ymax=15#max. position in y dirn
    x=np.linspace(1,xmax+1,Nx+1) #define grid in xdirn with Nx+1 nodes
    y=np.linspace(1,ymax+1,Ny+1) #define grid in ydirn with Ny+1 nodes
    Vguess=0 #initialize the guess of the potential to 0
    X,Y=np.meshgrid(x,y) #create 2d grid in x and y space
    V=np.zeros_like(X) #create a zero matrix with dimensions of X
    V.fill(Vguess) #fill each element of V with Vguess
    return V,X,Y


def capacitor(v,MyInput):
    #capacitor plate dimensions and positioning
    
    #set the voltage across plates by filling the elements of V with the plate voltage value, and using the range from 
    #
    if MyInput=='1':
        
        v[plate1x,(midpointy-halfplatelength):(midpointy+halfplatelength)] = -10
        
        v[plate2x,(midpointy-halfplatelength):(midpointy+halfplatelength)] = 10
    elif MyInput=='2':
        
        v[plate1x,(midpointy-halfplatelength):(midpointy+halfplatelength)] = 0
    
        v[plate2x,(midpointy-halfplatelength):(midpointy+halfplatelength)] = 0

    return v 


#set B.C.s
def BCs(v): 
    V_t=0.0 #define what values to assign the boundaries of the grid. 
    V_b=0.0
    V_l=0.0 
    V_r=0.0 
    
    v[-1,:] = V_t  #top of voltage grid
    v[0,:] = V_b #bottom of voltage grid
    v[:,0] = V_l #left of voltage grid
    v[:,-1] = V_r #right of voltage grid
    return v

def solveV(v,rel_tol,tol): 
    iterations=0
    while (rel_tol > tol):  #convergence condition
        dVmax=0.0
        for i in range(1,n+1): 

            for j in range(1, n+1):
                if i == plate1x or i==plate2x: #do not calculate voltage values for nodes on which plate is located
                    continue
                if j <(midpointy-halfplatelength) and j>(midpointy+halfplatelength):
                    continue

                dV= 0.25*(v[i,j-1]+v[i-1,j]+v[i,j+1]+v[i+1,j]-4.0*v[i,j]) #Gauss Seidel finite difference formula
                v[i,j]+=dV #update the grid node values with dV
                dVmax=np.max([np.abs(dV),dVmax])  
                rel_tol=dVmax/np.max(np.abs(v)) 
                iterations+=1
    print('iteration no.',iterations)
    print('distance between plates is',d,'m')
    print('plate length is',platelength,'m')


    return v


def matrix(N,case): #create iteration matrix 
    h=xmax/(N) #grid spacing
    alpha=0.16 #thermal diffusivity
    dt=0.1 #timestep


    M=np.eye(N+1) 
    M=M*(1+((2*alpha*dt)/h**2)) #set values for the diagonal elements 



    for i in range(2,N+1): #set values for the nearest neighbours of the diagonals on the RHS

        M[i-1,i]=-(alpha*dt)/h**2

    for i in range(0,N-1): #set values for the nearest neighbours of the diagonals on the LHS
        M[i+1,i]=-(alpha*dt)/h**2
    
    if case=='1': #sets the matrix values for no heat loss at end of rod
        M[0,0]=1 #fixes furthest left node of rod to 1000degrees
        M[N,N-1]=(-(alpha*dt)/h**2)*2 #imposes 0 gradient at furthest right node of rod

    elif case=='2': #set matrix values for end submerged in ice bath
        M[0,0]=1
        M[N,N]=1  #fixes  furthest right node of rod to 0 degrees 

    return M



def soln(t): 
    b=np.zeros((N+1,1)) #create 1d grid
    #set b.c.s
    b[:,:]=20
    b[0,:]=1000



    soln = np.zeros((N+1, t)) #create empty matrix to append values 
    for j in range(0, t, 1):  
        
        LU, P = lin.lu_factor(matrix(N,case)) #solve using LU decomposition
        vector = lin.lu_solve((LU, P), b) #store solution in 'vector'
        b = vector  #update value of b
        flat_list = [] #empty array to append values
        for sublist in vector:
            for item in sublist: #index each element of the vector array and append to flat_list
                flat_list.append(item)
        soln[:, j] = flat_list #assign flat_list to jth element of soln

    return  soln[:,t-1] 

def testsoln(t): 
    b=np.zeros((N+1,1)) #create 1d grid
    #set b.c.s
    b[:,:]=0
    b[0,:]=0



    soln = np.zeros((N+1, t)) #create empty matrix to append values 
    for j in range(0, t, 1):  
        
        LU, P = lin.lu_factor(matrix(N,case)) #solve using LU decomposition
        vector = lin.lu_solve((LU, P), b) #store solution in 'vector'
        b = vector  #update value of b
        flat_list = [] #empty array to append values
        for sublist in vector:
            for item in sublist: #index each element of the vector array and append to flat_list
                flat_list.append(item)
        soln[:, j] = flat_list #assign flat_list to jth element of soln

    return  soln[:,t-1] 




MyInput = '0'

while MyInput != 'q':
    MyInput = input('Enter a choice:\n "1" to see capacitor problem \n "2" to test the capacitor problem \n "3" to see the temperature distribution along a rod \n "4" to test the temperature distribution along rod \n "q" to quit:')


    if MyInput == '1': 
        print("You have chosen to see the contour of electric potential and vector field, and contour of electric field for a parallel plate capacitor")

            #define global variables and call functions
        Nx=Ny=15
        n=Nx-1


        tol=1e-19
        rel_tol=1e-10
        d=5
        platelength=10  #length of the plates
        halfplatelength=int(floor(platelength/2)) #takes the lower bound of half the plate length so as not to get a decimal number
        halfd=d/2 #half the distance between plates
        midpointx=int(ceil(Nx/2)) #takes the upper bound of half the number of nodes in x dirn
        midpointy=int(ceil(Ny/2))
        plate1x=int(midpointx-halfd) #positions plate 1 in x dirn
        plate2x=int(midpointx+halfd) #positions plate 2 in x dirn
        grid1=grid(Nx) #assigns the function to a variable
        v=grid1[0] #assign the 1st element of grid1 to v
        X=grid1[1] #assign the 2nd element of grid1 to X
        Y=grid1[2] #assign the 3rd element of grid1 to Y
        v=capacitor(v,'1') #sets the voltage values for the capacitor plates
        v=BCs(v) #sets edge boundary conditions on v
        v=solveV(v,rel_tol,tol) #calculates the voltage values for v 
        colorinterpolation = 100 
        colourMap = plt.cm.jet
    
    #contour plot of potential
        plt.title("Electric potential and Electric field ")
        plt.contourf(X, Y, v, colorinterpolation, cmap=colourMap)
        plt.colorbar() 
        plt.ylabel('y (m)')
        plt.xlabel('x (m)')
    

        #quiver plot of electric field
        gradients = np.gradient(-v) #finds the (negative) gradient of v
            
            
        E_x = gradients[0] #finds electric field in x and y dirn
        E_y = gradients[1]
            
            
        quiv = plt.quiver(X, Y, E_y, E_x, color='Black',headlength=1,headwidth=3,scale=50,pivot='middle')
        plt.show()

    
    #  contour plot for electric field
    
        plt.title("Electric field,E (V/m)")
                
        plt.ylabel('y (m)')
        plt.xlabel('x (m)')
        colorinterpolation = 45
        
        plt.contourf(X,Y,E_x,100,cmap=plt.cm.coolwarm_r)  #plot only the E_x values as only care about the dirn of the field between the plates. (In fact describes field in y dirn as it is a matrix)
        plt.colorbar()
            
        plt.show() 
    elif MyInput=='2': #tests the coding routine against simple, known case (p.d. of 0v between capacitor plates)
        print("You have chosen to test the capacitor code by setting all grid nodes to 0")
        tol=1e-19
        rel_tol=1e-10
        d=5
        platelength=10 #length of the plates
        halfplatelength=int(floor(platelength/2)) #takes the lower bound of half the plate length so as not to get a decimal number
        halfd=d/2 #half the distance between plates
        midpointx=int(ceil(Nx/2)) #takes the upper bound of half the number of nodes in x dirn
        midpointy=int(ceil(Ny/2))
        plate1x=int(midpointx-halfd) #positions plate 1 in x dirn
        plate2x=int(midpointx+halfd) #positions plate 2 in x dirn
        grid1=grid(Nx)
        v=grid1[0] 
        X=grid1[1]
        Y=grid1[2]
        v=capacitor(v,'2')
        v=BCs(v)
        v=solveV(v,rel_tol,tol)
        colorinterpolation = 100
        colourMap = plt.cm.jet
    
    
        plt.title("Electric potential and Electric field ")
        
            #colourMap = plt.cm.hot
        plt.contourf(X, Y, v, colorinterpolation, cmap=colourMap)
            # Set Colorbar
            
        plt.colorbar()
            
        plt.ylabel('y (m)')
        plt.xlabel('x (m)')
            
        gradients = np.gradient(-v)
            
            
        E_x = gradients[0]
        E_y = gradients[1]
            
            
        quiv = plt.quiver(X, Y, E_y, E_x, color='Black',headlength=1,headwidth=3,scale=50,pivot='middle')
        plt.show()

    
    

        
            
    elif MyInput == '3': 
        print("You have chosen to see the temperature distribution along a rod")
        NewInput=input('Enter "1" to see case 1 or "2" to see case 2:')
        if NewInput== '1':
            print("You have chosen to see the temperature distribution for no heat loss at end of rod")
            N=30
            case='1' 
            t = 10000
            xmax=50
            x=np.linspace(0,xmax,N+1)
            for l in range(0,6):
                
                t=10**l
                plt.plot(x, soln(t), label= str(t) + 's')
            plt.xlabel("Distance along rod (cm)")
            plt.ylabel("Temperature (°C)")
            plt.legend()
            pylab.legend(loc='center right')
            plt.show() 
        elif NewInput== '2':
            print("You have chosen to see the temperature distribution for the far RHS of rod submerged in ice at 0degrees C")
            N=30
            case='2' 
            t = 10000
            xmax=50
            x=np.linspace(0,xmax,N+1)
            for l in range(0,6):
                
                t=10**l
                plt.plot(x, soln(t), label= str(t) + 's')
            plt.xlabel("Distance along rod (cm)")
            plt.ylabel("Temperature (°C)")
            plt.legend()
            plt.show() 
    elif MyInput=='4': #tests code against simple, known case (boundary nodes set to 0)
            print("You have chosen to test the code for temperature distribution along a rod for case 1 (no heat loss at end of rod) by setting both ends of the rod to 0")
            N=30
            case='1'  #uses case 1 as the test case (case defines the matrix elements)
            t = 10000
            xmax=50
            x=np.linspace(0,xmax,N+1)
            for l in range(0,6): #for time between 1 and 10^6 sec, find the temperature distribution along the rod
                
                t=10**l #increases t by factor of 10 each iteration
                plt.plot(x, testsoln(t), label= str(t) + 's') #plots the distribution for each t
            plt.xlabel("Distance along rod (cm)")
            plt.ylabel("Temperature (°C)")
            plt.legend()
            plt.show() 
            


    elif MyInput != 'q':
        
        print('You have chosen to finish - goodbye.')