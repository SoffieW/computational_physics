import numpy as np
from scipy import linalg as lin
import time
import matplotlib.pyplot as plt
import pylab
import random
from math import sqrt
from mpl_toolkits.mplot3d import Axes3D

def generateempty(n): 
    '''generates empty matrix for the co-factors'''
    C=[]
    for i in range(n):
        C.append([0.0]*n) 
    
    for j in range(n):
        for i in range(n):
            C[i][j] = 0 #each element of the matrix set to 0
            
            return C
        
def minormatrix(M,i,j): 
    '''deletes the ith row and jth column of a matrix to obtain its matrix of minors'''
    Mnew=np.delete(M,i,axis=0) #deletes ith row
    return np.delete(Mnew,j,axis=1) #deletes jth column



def detnxn(M,n):
    '''computes the determinant of a nxn matrix with n>2'''
    if len(M)==2: #if the length of the matrix is 2, returns the determinant of that 2x2 matrix.
        return M[0,0]*M[1,1]-M[1,0]*M[0,1]
    det=0
    for j in range(0,n): #for every column in the matrix to be inverted,
        J=minormatrix(M,0,j) #delete the jth column and the top row to get the matrix of minors.

        det+=M[0,j]*detnxn(J,n-1)*((-1)**j) #recursive formula to compute the determinant. Put the new matrix J back into the function to find a its matrix of minors until length of matrix is 2, then its determinant will be found.
            
    return det



def cofactorsnxn(M,n,C): 
    '''computes the cofactors of a nxn matrix with n>2'''
    
    for i in range(0,n):  # for every element in the matrix to be inverted,
        for j in range(0,n): # delete ith row and jth column 
            J=minormatrix(M,i,j) 
            minordet=detnxn(J,n-1) #use the determinant func to find the det of the shrinking matrix until 2x2 matrix reached.
            C[i,j]=minordet*((-1)**(i+j))  #recursive formula to get the matirx of minors 
    return C



def cofactors2x2(M,n):
    '''computes the cofactors for 2x2 matrix'''
    M_c=np.array(generateempty(n)) #generate empty array
    M_c[0,0]=M[1,1] 
    M_c[1,1]=M[0,0]
    M_c[1,0]=-M[1,0]
    M_c[0,1]=-M[0,1] #find the elements of the cofactor matrix
    return M_c

def generatematrix(N): 
    '''generates NXN matrix'''
    myarray=[] #empty array stored in myarray
    for i in range(N):
        myarray.append([0.0]*N) #declare array of 0s with N rows and N columns
    
    for j in range(N): #sets N columns in matrix
        for i in range(N):
            myarray[i][j] = random.randint(-100,100) #generates random integer between 0 and 100 and assigns them to each element of the matrix
            
    return myarray

def generatevector(N):
    '''generates b vector'''
    b=[0 for y in range(N)] #generates an empty array with dimensions 1xN
    for y in range(N):
        b[y]=random.randint(-100,100) #assign a random integer between -100 and 100 to each element of array

    b=np.array(np.transpose([b])) #swaps rows and columns so get a vector i.e.:Nx1 matrix
    return b

def solveAIM(n):
    ''' returns the unknown vector for a system of equations'''

    C=np.array(generateempty(n),dtype='int64') #inserts empty matrix into C
    M=np.array(generatematrix(n), dtype='int64') #inserts the matrix to be inverted into M
       
    C=np.transpose(cofactorsnxn(M,n,C)) #transposes the matrix of cofactors
    b=generatevector(n)     #generates a b vector
    if len(M)==2:           #if M is a 2x2 matrix, solve the inverse using special cofactor function defined for 2x2 matrix
        inverse=(1/detnxn(M,n))*cofactors2x2(M,n)

    else:
        inverse=(1/(detnxn(M,n)))*C     
        
    vector=np.dot(inverse,b) #solve by doing dot product of inverse with b vector
    print('Your matrix is:')
    print(M)
    print('Your inverse is:')
    print(inverse)
    print('Your solution to the system of linear equations is:')
    print(vector)
    print('Solving using the builtin method gives:')
    print(np.dot(lin.inv(M),b)) #check AIM value by computing the builtin inverse
    
    return 

def generateAIMruntime():
    '''plots the runtime of AIM against N'''
    print('Program runtime for Analytical Inversion Method against N:')
    X=[] 
    Y=[]
    
    for n in range(2,9): #for number of rows between 2 and 9, 
        C=np.array(generateempty(n),dtype='int64') #inserts empty matrix into C
        M=np.array(generatematrix(n), dtype='int64') #inserts the matrix to be inverted into M
       
        C=np.transpose(cofactorsnxn(M,n,C)) 
        b=generatevector(n)
        start = time.time() #starts monitoring the runtime
        if len(M)==2: 
            inverse=(1/detnxn(M,n))*cofactors2x2(M,n)

        else:
            inverse=(1/(detnxn(M,n)))*C 

        vector=np.dot(inverse,b)
        end = time.time() #finishes monitoring the runtime
        Y.append(end-start) #appends the runtime for each n to array named Y
        X.append(n) #appends each n to array named X
                

                
    plt.plot(X,Y) #plots runtime against n

    plt.xlabel("N")
    plt.ylabel("time taken (s)")
    plt.show()
    return

def generatematrix(N): 
    '''generates NXN matrix'''
    myarray=[] #empty array stored in myarray
    for i in range(N):
        myarray.append([0.0]*N) #declare array of 0s with N rows and N columns
    
    for j in range(N): #sets N columns in matrix
        for i in range(N):
            myarray[i][j] = random.randint(-100,100) #generates random integer between 0 and 100 and assigns them to each element of the matrix
            
    return myarray

def generatevector(N):
    '''generates b vector'''
    b=[0 for y in range(N)] #generates an empty array with dimensions 1xN
    for y in range(N):
        b[y]=random.randint(-100,100) #assign a random integer between -100 and 100 to each element of array

    b=np.array(np.transpose([b])) #swaps rows and columns so get a vector i.e.:Nx1 matrix
    return b



    
def LUD(M,N,b): 
    '''returns the time it takes to solve a system of equations using LUD'''
    start = time.time()
    LU, P = lin.lu_factor(M) #find the Lower and upper triangular matrices and Permutation matrix of M
    vector=lin.lu_solve((LU,P),b) #rearrranges equation such that the unknown (vector) is the subject
    end = time.time()
    return end-start

def SVD(M,N,b):
    '''returns the time it takes to solve a system of equations using SVD'''
    start = time.time()
    u,s,v_transpose= lin.svd(M) #finds the unitary arrays u, v whose columns form an orthonormal basis and finds s, the singular value vectors
    s=lin.diagsvd(s,N,N) #construct the diagonal matrix from singular values
    s_inv=lin.inv(s) #computes the inverse of s
    u_transpose=np.transpose(u) #computes transpose of u
    v=np.transpose(v_transpose) #computes transpose of v
    
    vector=np.dot(v,np.dot(s_inv,np.dot(u_transpose,b))) #solves to find the unknown vector v
    end = time.time()
    return end-start

def LUD_vectors(k): 
    '''returns the vectors for a system of linear equations with varying k using LUD'''
    M=np.array([[1,1,1],[1,2,-1],[2,3,k]])  #define M with varying k
    b=np.transpose(np.array([[5,10,15]])) 
    LU, P = lin.lu_factor(M) 
    vector=lin.lu_solve((LU,P),b)
    return vector[0], vector[1], vector[2] #return the 3 vector elements individually

def SVD_vectors(k):
    '''returns the vectors for a system of linear equations with varying k using SVD'''
    M=np.array([[1,1,1],[1,2,-1],[2,3,k]]) 
    b=np.transpose(np.array([[5,10,15]])) 
    u,s,v_transpose= lin.svd(M)
    s=lin.diagsvd(s,3,3)
    s_inv=lin.inv(s)
    u_transpose=np.transpose(u)
    v=np.transpose(v_transpose)
    
    vector=np.dot(v,np.dot(s_inv,np.dot(u_transpose,b)))

    return vector[0], vector[1], vector[2]

def generateLUSVDruntime():
    '''Generates plot of runtime for LUD and SVD against number of rows N'''
    print('Program runtime for LU and SV decomposition against N:')
    LU=[] #creat empty arrays for dependent and independent variables
    x=[]
    SV=[]


    for N in np.arange(2,500,2): #for number of rows between 2 and 500, perform these commands:
        M=generatematrix(N) #generates random matrix of dimensions NxN
        b=generatevector(N) #generates random b vector of dimensions Nx1
        l=LUD(M,N,b) #assign the runtimes for lu over each iteration to l
        x.append(N) #append the N to array x
        LU.append(l) #append l to array LU
        m=SVD(M,N,b) #assigns the runtimes for svd over each iteration to m
        SV.append(m) #appends m to array SV
    y=LU
    y1=SV
    
#plot the runtimes for LU and SVD against number of rows N:
    plt.plot(x,y, label='LU')
    plt.plot(x,y1,label='SVD')
    plt.xlabel("N")
    plt.ylabel("time taken (s)")
    pylab.legend(loc='upper left')
    plt.show()

    return 

def generateSVDsingplot():
    '''generates a plot of SVD solutions with varying k'''
    print('Vector elements x,y and z using SVD against k (last element in M):')
    SVvx=[]
    SVvy=[]
    SVvz=[]
    x_s=[]
  #  for k in np.arange(1.0e-17,1e-14,1.0e-16): #range in which SVD method displays strange behaviour
    for k in np.arange(1e-5,1e-2,0.001):
        x_s.append(k)
        SVDvects=SVD_vectors(k)

        SVvx.append(SVDvects[0])
        
        SVvy.append(SVDvects[1])
        
        SVvz.append(SVDvects[2])

    
    
    plt.plot(x_s,SVvx, label='x')
    plt.plot(x_s,SVvy,label='y')
    plt.plot(x_s, SVvz,label='z')
    plt.xlabel("k")
    plt.ylabel("vector element")
    pylab.legend(loc='upper right')
    plt.show()
    return

def AIM_vectors(k):
       M=np.array([[1,1,1],[1,2,-1],[2,3,k]]) 
       b=np.transpose(np.array([[5,10,15]])) 
       C=np.array(generateempty(3),dtype='int64') #inserts empty matrix into C
       Cnew=np.transpose(cofactorsnxn(M,3,C))

       inverse=(1/(detnxn(M,3)))*Cnew

       vector=np.dot(inverse,b)
       return vector[0], vector[1], vector[2]
    
    
    
    
def generateAIMsingplot():
    '''generate a plot of AIM solns with varying k'''
    print('Vector elements x,y and z using AIM against k (last element in M):')
    x_s=[]
    AIMvx=[]
    AIMvy=[]
    AIMvz=[]
    

    for k in np.arange(0.1,10,0.001):
     
        x_s.append(k)
        
        AIMvects=AIM_vectors(k)
        
        AIMvx.append(AIMvects[0])
        
        AIMvy.append(AIMvects[1])
        
        AIMvz.append(AIMvects[2])


    
    plt.plot(x_s,AIMvx, label='x')
    plt.plot(x_s,AIMvy,label='y')
    plt.plot(x_s, AIMvz,label='z')
    plt.xlabel("k")
    plt.ylabel("vector element")
    pylab.legend(loc='upper right')
    plt.show()
    return

def tension2d(x,z):
    '''function which returns an array containing T1 and T2'''
    m=70 #mass of acrobat
    g=9.81 #acceleration due to gravity
    w=m*g #weight 
    x=float(x)
    z=float(z)

    matrix = np.array([[x/sqrt(x**2+(8-z)**2),(x-15)/sqrt((x-15)**2+(8-z)**2)],
                  [(8-z)/sqrt(x**2+(8-z)**2),(8-z)/sqrt((x-15)**2+(8-z)**2)]]) #matrix of coefficients in terms of position 
    u,s,v_transpose= lin.svd(matrix) 
    s=lin.diagsvd(s,2,2)
    s_inv=lin.inv(s)
    u_transpose=np.transpose(u)
    v=np.transpose(v_transpose)
    tension=np.dot(v,np.dot(s_inv,np.dot(u_transpose,np.array([0,w])))) #solve tension using SVD, and the b vector is [0, w]

    return tension

def tension3d(x,y,z):
    '''function which returns an array containing T1, T2 and T3'''
    g= 9.81
    m = 70
    w = m*g

    matrix = np.array([[x/sqrt((x**2+(8-y)**2+(8-z)**2)),(15-x)/sqrt((15-x)**2+(8-y)**2+(8-z)**2),(7.5-x)/sqrt((7.5-x)**2+(y)**2+(8-z)**2)],
             [(8-y)/sqrt((x**2+(8-y)**2+(8-z)**2)),(8-y)/sqrt((15-x)**2+(8-y)**2+(8-z)**2),y/sqrt((7.5-x)**2+(y)**2+(8-z)**2)],
             [(8-z)/sqrt((x**2+(8-y)**2+(8-z)**2)),(8-z)/sqrt((15-x)**2+(8-y)**2+(8-z)**2),(8-z)/sqrt((7.5-x)**2+(y)**2+(8-z)**2)]])
    
    u,s,v= lin.svd(matrix)
    u_transpose=np.transpose(u)
    s=np.diag(1/s)
    v_transpose=np.transpose(v)
    Tension= np.dot(v_transpose,np.linalg.solve(s,np.dot(u_transpose,np.array([0,0,w]))))

    return  Tension[0], Tension[1], Tension[2]



def plot2dtension():
    '''generates surface plots of the tension in wires 1 and 2, neglecting the 3rd wire.'''
        
    T1=[]
    T2=[]

    X =[]
    Y= []
    T1_pos=[]
    T2_pos=[]

    X_pos =[]
    Y_pos= []
    
    
      
    for x in np.arange(1,15,1):
        for y in np.arange(1,8,1):
            X.append(x)
            Y.append(y)
            t= tension2d(x,y)
            T1.append(t[0])
            T2.append(t[1])
            
    for x_pos in np.arange(1,15,0.1):
        for y_pos in np.arange(1,8,0.1):
            X_pos.append(x_pos)
            Y_pos.append(y_pos)
            t_pos= tension2d(x_pos,y_pos)
            T1_pos.append(t_pos[0])
            T2_pos.append(t_pos[1])
            
        
                

    x = X
    y= Y
    z1 = T1
    z2=T2

    ax = plt.axes(projection='3d')
    ax.plot_trisurf(x, y, z1, edgecolor='none')
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    ax.set_zlabel("T1 👎")
    plt.show()
    ax = plt.axes(projection='3d')
    ax.plot_trisurf(x, y, z2, edgecolor='none')
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    ax.set_zlabel("T2 👎")
    plt.show()
#        ax.zaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        
    print('Max. tension in 1st wire:', max(T1),'N') 
    print('Max. tension in 2nd wire:', max(T2),'N') 
    position1=int(T1_pos.index(max(T1_pos)))
    position2=int(T2_pos.index(max(T2_pos)))
    print('Position of max tension in 1st wire:', 'x=', X_pos[position1],'m', 'z=',Y_pos[position1],'m')
    print('Position of max tension in 2nd wire:', 'x=', X_pos[position2],'m', 'z=',Y_pos[position2],'m')
    return



def plot3dtension():
    '''generates scatter plot of the tension in T1,T2,T3 with position on the axes, and colour map representing tension.'''
    X=[]
    Y=[]
    Z=[]
    T1=[]
    T2=[]
    T3=[]


    for x in np.arange(1,15,0.5):
        for y in np.arange(1,8,0.5):
            for z in np.arange(1,8,0.5):
                if y > (8/7.5)*x:
                    continue
                if y > -((8/7.5)*x) + 16:
                    continue
        
        
                t= tension3d(x,y,z)
                
                X.append(x)
                Y.append(y)
                Z.append(z)
                T1.append(t[0])
                T2.append(t[1])
                T3.append(t[2])     
        
    x=X
    y=Y
    z=Z


# 3D Plot
    fig = plt.figure()
    ax3D = fig.add_subplot(111, projection='3d')
    ax3D.scatter(x, y, z, s=10, c=T1, marker='o')  
    plt.xlabel("x (m)")
    plt.ylabel(" y (m)")
    plt.show()


# 3D Plot
    fig = plt.figure()
    ax3D = fig.add_subplot(111, projection='3d')
    ax3D.scatter(x, y, z, s=10, c=T2, marker='o')  
    plt.xlabel("x (m)")
    plt.ylabel(" y (m)")
    plt.show()
# 3D Plot
    fig = plt.figure()
    ax3D = fig.add_subplot(111, projection='3d')
    ax3D.scatter(x, y, z, s=10, c=T3, marker='o')  
    plt.xlabel("x (m)")
    plt.ylabel(" y (m)")  

    plt.show()
    return

def solve3dtension():
    '''finds the maximum values of T1,T2,T3 and the position at which this occurs.''' 
    
    X=[]
    Y=[]
    Z=[]
    T1=[]
    T2=[]
    T3=[]
    for x in np.arange(1,15,1):  #increment the position by 1 over the ranges for x,y,z and find the tensions of the 3 wires at each point.
        for y in np.arange(1,8,1):
            for z in np.arange(1,8,1):
                if y > (8/7.5)*x: #sets boundary conditions to restrict the acrobat's position- he cannot go outside triangle formed by the 3 wire drums, as then tension would become negative.
                    continue
                if y > -((8/7.5)*x) + 16:
                    continue #for these regions defined by the if statement, do not append values 
        
        
                t= tension3d(x,y,z) 
                
                X.append(x)
                Y.append(y)
                Z.append(z)
                T1.append(t[0])
                T2.append(t[1])
                T3.append(t[2])     
        
    x=X
    y=Y
    z=Z    
    X_pos=[]
    Y_pos=[]
    Z_pos=[]
    T1_pos=[]
    T2_pos=[]
    T3_pos=[]
    for x_pos in np.arange(1,15,0.5): #iterations to find position- require smaller increments of 0.5 
        for y_pos in np.arange(1,8,0.5):
            for z_pos in np.arange(1,8,0.5):
                if y_pos > (8/7.5)*x_pos:
                    continue
                if y_pos > -((8/7.5)*x_pos) + 16:
                    continue
        
        
                t_pos= tension3d(x_pos,y_pos,z_pos)
                
                X_pos.append(x_pos)
                Y_pos.append(y_pos)
                Z_pos.append(z_pos)
                T1_pos.append(t_pos[0])
                T2_pos.append(t_pos[1])
                T3_pos.append(t_pos[2])     
        
    x_pos=X_pos
    y_pos=Y_pos
    z_pos=Z_pos 
    print('Max. tension in 1st wire:', max(T1),'N') 
    print('Max. tension in 2nd wire:', max(T2),'N') 
    print('Max.tension in 3rd wire:',max(T3),'N')
    position1=int(T1_pos.index(max(T1_pos))) #assigns the maximum value in T1_pos array to position1
    position2=int(T2_pos.index(max(T2_pos)))
    position3=int(T3_pos.index(max(T3_pos)))
    print('Position of max tension in 1st wire:', 'x=', X_pos[position1],'m', 'y=',Y_pos[position1],'m','z=',8-Z_pos[position1],'m')
    print('Position of max tension in 2nd wire:', 'x=', X_pos[position2],'m', 'y=',Y_pos[position2],'m','z=',8-Z_pos[position2],'m')
    print('Position of max tension in 3rd wire:', 'x=', X_pos[position3],'m', 'y=',Y_pos[position3],'m','z=',8-Z_pos[position3],'m')
    return


MyInput = '0'

while MyInput != 'q':
    MyInput = input('Enter a choice, "a", "b", "c", "d", "e" or "q" to quit: ')


    if MyInput == 'a': 
        myinput=input('Enter a value for N between 2 and 9:')
        n=int(myinput)
        print(solveAIM(n))
        
    elif MyInput == 'b': 
        print(generateAIMruntime())
        print(generateLUSVDruntime())

    elif MyInput== 'c':
        print(generateAIMsingplot())
        print(generateSVDsingplot())
        
 
    elif MyInput== 'd':
        print(plot2dtension())
        
    elif MyInput=='e':
        print(solve3dtension())
        print(plot3dtension())



    elif MyInput != 'q':
        
        print('You have chosen to finish - goodbye.')