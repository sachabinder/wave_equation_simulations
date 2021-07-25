"""
This file was built to solve numerically a classical PDE, 2D wave equation. The equation corresponds to 

$\dfrac{\partial}{\partial x} \left( \dfrac{\partial c^2 U}{\partial x} \right) + \dfrac{\partial}{\partia
l y} \left( \dfrac{\partial c^2 U}{\partial y} \right) = \dfrac{\partial^2 U}{\partial t^2}$

where
 - U represent the signal
 - x represent the position
 - t represent the time
 - c represent the velocity of the wave (depends on space parameters)

The numerical scheme is based on finite difference method. This program is also providing several boundary conditions. More particularly the Neumann, Dirichlet and Mur boundary conditions.
Copyright - © SACHA BINDER - 2021
"""

############## MODULES IMPORTATION ###############
import numpy as np
import matplotlib.pyplot as plt
import viz_tools    #self-developed module that groups animation functions



#Def of the initial condition   
def I(x,y):
    """
    two space variables depending function 
    that represent the wave form at t = 0
    """
    return 0.2*np.exp(-((x-1)**2/0.1 + (y-1)**2/0.1))

def V(x,y):
    """
    initial vertical speed of the wave
    """
    return 0
    
    
    
############## SET-UP THE PROBLEM ###############

#Def of velocity (spatial scalar field)
def celer(x,y):
    """
    constant velocity
    """
    return 1
    
loop_exec = 1 # Processing loop execution flag

bound_cond = 2  #Boundary cond 1 : Dirichlet, 2 : Neumann, 3 Mur

if bound_cond not in [1,2,3]:
    loop_exec = 0
    print("Please choose a correct boundary condition")



#Spatial mesh - i indices
L_x = 5 #Range of the domain according to x [m]
dx = 0.05 #Infinitesimal distance in the x direction
N_x = int(L_x/dx) #Points number of the spatial mesh in the x direction
X = np.linspace(0,L_x,N_x+1) #Spatial array in the x direction

#Spatial mesh - j indices
L_y = 5 #Range of the domain according to y [m]
dy = 0.05 #Infinitesimal distance in the x direction
N_y = int(L_y/dy) #Points number of the spatial mesh in the y direction
Y = np.linspace(0,L_y,N_y+1) #Spatial array in the y direction



#Temporal mesh with CFL < 1 - n indices
L_t = 4 #Duration of simulation [s]
dt = dt = 0.1*min(dx, dy)   #Infinitesimal time with CFL (Courant–Friedrichs–Lewy condition)
N_t = int(L_t/dt) #Points number of the temporal mesh
T = np.linspace(0,L_t,N_t+1) #Temporal array

#Velocity array for calculation (finite elements)
c = np.zeros((N_x+1,N_y+1), float)
for i in range(0,N_x+1):
    for j in range(0,N_y+1):
        c[i,j] = celer(X[i],Y[j])




############## CALCULATION CONSTANTS ###############
Cx2 = (dt/dx)**2
Cy2 = (dt/dy)**2 
CFL_1 = dt/dy*c[:,0]
CFL_2 = dt/dy*c[:,N_y]
CFL_3 = dt/dx*c[0,:]
CFL_4 =dt/dx*c[N_x,:]




############## PROCESSING LOOP ###############

if loop_exec:
    # $\forall i \in {0,...,N_x}$
    U = np.zeros((N_x+1,N_x+1,N_t+1),float) #Tableau de stockage de la solution

    u_nm1 = np.zeros((N_x+1,N_y+1),float)   #Vector array u_{i,j}^{n-1}
    u_n = np.zeros((N_x+1,N_y+1),float)     #Vector array u_{i,j}^{n}
    u_np1 = np.zeros((N_x+1,N_y+1),float)  #Vector array u_{i,j}^{n+1}
    V_init = np.zeros((N_x+1,N_y+1),float)
    q = np.zeros((N_x+1, N_y+1), float)
    
    #init cond - at t = 0
    for i in range(0, N_x+1):
        for j in range(0, N_y+1):
            q[i,j] = c[i,j]**2
    
    for i in range(0, N_x+1):
        for j in range(0, N_y+1):
            u_n[i,j] = I(X[i],Y[j])
            
    for i in range(0, N_x+1):
        for j in range(0, N_y+1):
            V_init[i,j] = V(X[i],Y[j])
    
    U[:,:,0] = u_n.copy()

    

    #init cond - at t = 1
    #without boundary cond
    u_np1[1:N_x,1:N_y] = 2*u_n[1:N_x,1:N_y] - (u_n[1:N_x,1:N_y] - 2*dt*V_init[1:N_x,1:N_y]) + Cx2*(  0.5*(q[1:N_x,1:N_y] + q[2:N_x+1,1:N_y ])*(u_n[2:N_x+1,1:N_y] - u_n[1:N_x,1:N_y])  - 0.5*(q[0:N_x -1,1:N_y] + q[1:N_x,1:N_y ])*(u_n[1:N_x,1:N_y] - u_n[0:N_x -1,1:N_y]) ) + Cy2*(  0.5*(q[1:N_x,1:N_y] + q[1:N_x ,2:N_y+1])*(u_n[1:N_x,2:N_y+1] - u_n[1:N_x,1:N_y])  - 0.5*(q[1:N_x,0:N_y -1] + q[1:N_x ,1:N_y])*(u_n[1:N_x,1:N_y] - u_n[1:N_x,0:N_y -1]) )


    #boundary conditions
    if bound_cond == 1:
        #Dirichlet bound cond
        u_np1[0,:] = 0
        u_np1[-1,:] = 0
        u_np1[:,0] = 0
        u_np1[:,-1] = 0



    elif bound_cond == 2:
        #Nuemann bound cond
        i,j = 0,0
        u_np1[i,j] = 2*u_n[i,j] - (u_n[i,j] - 2*dt*V_init[i,j]) + Cx2*(q[i,j] + q[i+1,j])*(u_n[i+1,j] - u_n[i,j]) + Cy2*(q[i,j] + q[i,j+1])*(u_n[i,j+1] - u_n[i,j])
        
        i,j = 0,N_y
        u_np1[i,j] = 2*u_n[i,j] - (u_n[i,j] - 2*dt*V_init[i,j]) + Cx2*(q[i,j] + q[i+1,j])*(u_n[i+1,j] - u_n[i,j]) + Cy2*(q[i,j] + q[i,j-1])*(u_n[i,j-1] - u_n[i,j])
                    
        i,j = N_x,0
        u_np1[i,j] = 2*u_n[i,j] - (u_n[i,j] - 2*dt*V_init[i,j]) + Cx2*(q[i,j] + q[i-1,j])*(u_n[i-1,j] - u_n[i,j]) + Cy2*(q[i,j] + q[i,j+1])*(u_n[i,j+1] - u_n[i,j])
                
        i,j = N_x,N_y
        u_np1[i,j] = 2*u_n[i,j] - (u_n[i,j] - 2*dt*V_init[i,j]) + Cx2*(q[i,j] + q[i-1,j])*(u_n[i-1,j] - u_n[i,j]) + Cy2*(q[i,j] + q[i,j-1])*(u_n[i,j-1] - u_n[i,j])      
        
        i = 0
        u_np1[i,1:N_y -1] = 2*u_n[i,1:N_y -1] - (u_n[i,1:N_y -1] - 2*dt*V_init[i,1:N_y -1]) + Cx2*(q[i,1:N_y -1] + q[i+1,1:N_y -1])*(u_n[i+1,1:N_y -1] - u_n[i,1:N_y -1]) + Cy2*(  0.5*(q[i,1:N_y -1] + q[i,2:N_y])*(u_n[i,2:N_y] - u_n[i,1:N_y -1])  - 0.5*(q[i,0:N_y -2] + q[i,1:N_y -1])*(u_n[i,1:N_y -1] - u_n[i,0:N_y -2]) )
              
        j = 0
        u_np1[1:N_x -1,j] = 2*u_n[1:N_x -1,j] - (u_n[1:N_x -1,j] - 2*dt*V_init[1:N_x -1,j]) + Cx2*(  0.5*(q[1:N_x -1,j] + q[2:N_x,j])*(u_n[2:N_x,j] - u_n[1:N_x -1,j])  - 0.5*(q[0:N_x -2,j] + q[1:N_x -1,j])*(u_n[1:N_x -1,j] - u_n[0:N_x -2,j]) ) + Cy2*(q[1:N_x -1,j] + q[1:N_x -1,j+1])*(u_n[1:N_x -1,j+1] - u_n[1:N_x -1,j])
       
        i = N_x
        u_np1[i,1:N_y -1] = 2*u_n[i,1:N_y -1] - (u_n[i,1:N_y -1] - 2*dt*V_init[i,1:N_y -1]) + Cx2*(q[i,1:N_y -1] + q[i-1,1:N_y -1])*(u_n[i-1,1:N_y -1] - u_n[i,1:N_y -1]) + Cy2*(  0.5*(q[i,1:N_y -1] + q[i,2:N_y])*(u_n[i,2:N_y] - u_n[i,1:N_y -1])  - 0.5*(q[i,0:N_y -2] + q[i,1:N_y -1])*(u_n[i,1:N_y -1] - u_n[i,0:N_y -2]) )
              
        j = N_y
        u_np1[1:N_x -1,j] = 2*u_n[1:N_x -1,j] - (u_n[1:N_x -1,j] - 2*dt*V_init[1:N_x -1,j]) + Cx2*(  0.5*(q[1:N_x -1,j] + q[2:N_x,j])*(u_n[2:N_x,j] - u_n[1:N_x -1,j])  - 0.5*(q[0:N_x -2,j] + q[1:N_x -1,j])*(u_n[1:N_x -1,j] - u_n[0:N_x -2,j]) ) + Cy2*(q[1:N_x -1,j] + q[1:N_x -1,j-1])*(u_n[1:N_x -1,j-1] - u_n[1:N_x -1,j])
               
               
    elif bound_cond == 3:
        #Nuemann bound cond
        i = 0
        u_np1[i,:] = u_n[i+1,:] + (CFL_3 - 1)/(CFL_3 + 1)*(u_np1[i+1,:] - u_n[i,:])
        
        j = 0
        u_np1[:,j] = u_n[:,j+1] + (CFL_1 - 1)/(CFL_1 + 1)*(u_np1[:,j+1] - u_n[:,j])
        
        i = N_x
        u_np1[i,:] = u_n[i-1,:] + (CFL_4 - 1)/(CFL_4 + 1)*(u_np1[i-1,:] - u_n[i,:])
        
        j = N_y
        u_np1[:,j] = u_n[:,j-1] + (CFL_2 - 1)/(CFL_2 + 1)*(u_np1[:,j-1] - u_n[:,j])
        
        
    
    
    u_nm1 = u_n.copy()
    u_n = u_np1.copy()
    U[:,:,1] = u_n.copy()

    
    #Process loop (on time mesh)
    for n in range(2, N_t):
        
        #calculation at step j+1  
        #without boundary cond           
        u_np1[1:N_x,1:N_y] = 2*u_n[1:N_x,1:N_y] - u_nm1[1:N_x,1:N_y] + Cx2*(  0.5*(q[1:N_x,1:N_y] + q[2:N_x+1,1:N_y])*(u_n[2:N_x+1,1:N_y] - u_n[1:N_x,1:N_y])  - 0.5*(q[0:N_x - 1,1:N_y] + q[1:N_x,1:N_y])*(u_n[1:N_x,1:N_y] - u_n[0:N_x - 1,1:N_y]) ) + Cy2*(  0.5*(q[1:N_x ,1:N_y] + q[1:N_x,2:N_y+1])*(u_n[1:N_x,2:N_y+1] - u_n[1:N_x,1:N_y])  - 0.5*(q[1:N_x,0:N_y - 1] + q[1:N_x,1:N_y])*(u_n[1:N_x,1:N_y] - u_n[1:N_x,0:N_y - 1]) )
            
            
            
        #bound conditions
        if bound_cond == 1:
            #Dirichlet bound cond
            u_np1[0,:] = 0
            u_np1[-1,:] = 0
            u_np1[:,0] = 0
            u_np1[:,-1] = 0
            
        
        elif bound_cond == 2:
            #Nuemann bound cond
            i,j = 0,0
            u_np1[i,j] = 2*u_n[i,j] - u_nm1[i,j] + Cx2*(q[i,j] + q[i+1,j])*(u_n[i+1,j] - u_n[i,j]) + Cy2*(q[i,j] + q[i,j+1])*(u_n[i,j+1] - u_n[i,j])
                        
            i,j = 0,N_y
            u_np1[i,j] = 2*u_n[i,j] - u_nm1[i,j] + Cx2*(q[i,j] + q[i+1,j])*(u_n[i+1,j] - u_n[i,j]) + Cy2*(q[i,j] + q[i,j-1])*(u_n[i,j-1] - u_n[i,j])
                            
            i,j = N_x,0
            u_np1[i,j] = 2*u_n[i,j] - u_nm1[i,j] + Cx2*(q[i,j] + q[i-1,j])*(u_n[i-1,j] - u_n[i,j]) + Cy2*(q[i,j] + q[i,j-1])*(u_n[i,j-1] - u_n[i,j])
                    
            i,j = N_x,N_y
            u_np1[i,j] = 2*u_n[i,j] - u_nm1[i,j] + Cx2*(q[i,j] + q[i-1,j])*(u_n[i-1,j] - u_n[i,j]) + Cy2*(q[i,j] + q[i,j-1])*(u_n[i,j-1] - u_n[i,j])
                    
                    
            i = 0
            u_np1[i,1:N_y -1] = 2*u_n[i,1:N_y -1] - u_nm1[i,1:N_y -1] + Cx2*(q[i,1:N_y -1] + q[i+1,1:N_y -1])*(u_n[i+1,1:N_y -1] - u_n[i,1:N_y -1]) + Cy2*(  0.5*(q[i,1:N_y -1] + q[i,2:N_y])*(u_n[i,2:N_y] - u_n[i,1:N_y -1])  - 0.5*(q[i,0:N_y -2] + q[i,j])*(u_n[i,1:N_y -1] - u_n[i,0:N_y -2]) )
                        
            j = 0
            u_np1[1:N_x - 1,j] = 2*u_n[1:N_x - 1,j] - u_nm1[1:N_x - 1,j] + Cx2*(  0.5*(q[1:N_x - 1,j] + q[2:N_x,j])*(u_n[2:N_x,j] - u_n[1:N_x - 1,j])  - 0.5*(q[0:N_x - 2,j] + q[1:N_x - 1,j])*(u_n[1:N_x - 1,j] - u_n[0:N_x - 2,j]) ) + Cy2*(q[1:N_x - 1,j] + q[1:N_x - 1,j+1])*(u_n[1:N_x - 1,j+1] - u_n[1:N_x - 1,j])
                    
            i = N_x
            u_np1[i,1:N_y -1] = 2*u_n[i,1:N_y -1] - u_nm1[i,1:N_y -1] + Cx2*(q[i,1:N_y -1] + q[i-1,1:N_y -1])*(u_n[i-1,1:N_y -1] - u_n[i,1:N_y -1]) + Cy2*(  0.5*(q[i,1:N_y -1] + q[i,2:N_y])*(u_n[i,2:N_y] - u_n[i,1:N_y -1])  - 0.5*(q[i,0:N_y -2] + q[i,1:N_y -1])*(u_n[i,1:N_y -1] - u_n[i,0:N_y -2]) )
                    
            j = N_y
            u_np1[1:N_x - 1,j] = 2*u_n[1:N_x - 1,j] - u_nm1[1:N_x - 1,j] + Cx2*(  0.5*(q[1:N_x - 1,j] + q[2:N_x,j])*(u_n[2:N_x,j] - u_n[1:N_x - 1,j])  - 0.5*(q[0:N_x - 2,j] + q[1:N_x - 1,j])*(u_n[1:N_x - 1,j] - u_n[0:N_x - 2,j]) ) + Cy2*(q[1:N_x - 1,j] + q[1:N_x - 1,j-1])*(u_n[1:N_x - 1,j-1] - u_n[1:N_x - 1,j])
                

        elif bound_cond == 3:
            #Mur bound cond
            i = 0
            u_np1[i,:] = u_n[i+1,:] + (CFL_3 - 1)/(CFL_3 + 1)*(u_np1[i+1,:] - u_n[i,:])
            
            j = 0
            u_np1[:,j] = u_n[:,j+1] + (CFL_1 - 1)/(CFL_1 + 1)*(u_np1[:,j+1] - u_n[:,j])
            
            i = N_x
            u_np1[i,:] = u_n[i-1,:] + (CFL_4 - 1)/(CFL_4 + 1)*(u_np1[i-1,:] - u_n[i,:])
            
            j = N_y
            u_np1[:,j] = u_n[:,j-1] + (CFL_2 - 1)/(CFL_2 + 1)*(u_np1[:,j-1] - u_n[:,j])

        
        
        u_nm1 = u_n.copy()      
        u_n = u_np1.copy() 
        U[:,:,n] = u_n.copy()
        
######################### PLOT #############################

anim = viz_tools.anim_2D(X,Y,U,dt,5)
anim2 = viz_tools.anim_2D_flat(X,Y,U,dt,2)
plt.show()



