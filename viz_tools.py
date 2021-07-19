"""
File with several visualization functions intended to be used 
with results from 1D/2D wave equations simulation
"""

############## MODULES IMPORTATION ###############
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from scipy import meshgrid
import os

##################################################

def plot_a_frame_1D(x, y, xmin, xmax, ymin, ymax, titre = "My title", type = "-"):
    """
    Plot a random 1D solution with its plot window
    
    (x:np.ndarray (format 1D), y:np.ndarray (format 1D), xmin:float, xmax:float, ymin:float, ymax:float, titre:str, type:str) -> plot
    """
    plt.axis([xmin,xmax,ymin,ymax])
    plt.plot(x,y, type , color = "black")
    #plt.ion()
    plt.title(titre)
    plt.xlabel("x-axis [m]")
    plt.ylabel("y-axis [m]")
    plt.show()

##################################################

def plot_spatio_temp_3D(x,y,z):
    """
    Plot a 2 two parameters function z = f(x,t) where x-axis is spatial and y-axis is time.
    
    (x:np.ndarray (format 1D), y:np.ndarray (format 1D), z:np.ndarray (format 1D)) -> plot
    """
    fig = plt.figure(figsize=(14,8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel(r'$x \ [ $', fontsize = 16)
    ax.set_ylabel('temps', fontsize = 16)
    ax.set_zlabel('amplitude', fontsize = 16)
    ax.view_init(elev = 15, azim = 120)
    
    ST,SX = meshgrid(y,x)
    p = ax.plot_surface(SX,ST,z,color = 'white')       
    plt.show()

##################################################

def anim_1D(x,y, pas_de_temps, pas_d_images, save = False, myxlim = (0, 4) , myylim = (-4,4)):
    """
    Function allowing to display an annimation based on calculation result with a given time step. This function can be used to save the images sequence in the current directory.
    
    The y parameter is a list containing several functions to display : y = [ [f_1(x)], ... , [f_n(x)] ].
    
    (x:np.ndarray (format 1D), y:np.ndarray (format 2D), pas_de_temps:float , pas_d_images:int, save:bool , myxlim:tuple , myylim:tuple) -> plot (+ .mp4)
    """
    
    
    fig = plt.figure()
    ax = plt.axes(xlim= myxlim , ylim= myylim)
    line, = ax.plot([], [])
    ax.set_title("t = 0 s", fontname = "serif", fontsize = 16)
    ax.set_xlabel("x [m]", fontname = "serif", fontsize = 14)
    ax.set_ylabel("$u$ [m]", fontname = "serif", fontsize = 14)
    def init():
        line.set_data([],[])
        return line,
    
    # animation function.  This is called sequentially
    def animate(i):
        line.set_data(x, y[:,pas_d_images*i])
        ax.set_title("$u(x)$ à t = {} s".format(np.round(i*pas_d_images*pas_de_temps, 4)), fontname = "serif", fontsize = 16)
        return line,
        
    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=y.shape[1]//pas_d_images, interval=10, blit=True)

    if save:
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=25, metadata=dict(artist='Me'), bitrate=1800)

        anim.save('lines.mp4', writer=writer)

    return anim
    
##################################################

def anim_2D(X, Y, L, pas_de_temps, pas_d_images, save = False, myzlim = (-0.15, 0.15)):
    """
    Fonction qui peremt d'annimer un représentation graphique en 3D où z(x,y). Pour cela, on stock dans la liste L l'ensembles des images à affiicher. L = [ [z_1(x,y)], ... , [z_n(x,y)] ].
    On peut éventuellement enregistrer l'annimation.
    
    (X:np.ndarray (format 1D), Y:np.ndarray (format 1D), L:np.ndarray (format 3D), pas_de_temps:float, pas_d_images:int, save:bool, myzlim:tuple) -> plot (+ .mp4)
    """
    
    fig = plt.figure(figsize = (8, 8), facecolor = "white")
    ax = fig.add_subplot(111, projection='3d')
    SX,SY = meshgrid(X,Y)
    surf = ax.plot_surface(SX, SY, L[:,:,0],cmap = plt.cm.RdBu_r)
    ax.set_zlim(myzlim[0], myzlim[1])
    ax.set_title("t = 0 s", fontname = "serif", fontsize = 16)
    
    # animation function.  This is called sequentially
    def update_surf(num):
        ax.clear()
        surf = ax.plot_surface(SX, SY, L[:,:,pas_d_images*num],cmap = plt.cm.viridis)
        ax.set_xlabel("x [m]", fontname = "serif", fontsize = 14)
        ax.set_ylabel("y [m]", fontname = "serif", fontsize = 14)
        ax.set_zlabel("$u$ [m]", fontname = "serif", fontsize = 16)
        ax.set_title("$u(x,y)$ à t = {} s".format(np.round(pas_d_images*num*pas_de_temps, 4)), fontname = "serif", fontsize = 16)
        ax.set_zlim(myzlim[0], myzlim[1])
        plt.tight_layout()
        return surf,
        
    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, update_surf, frames = L.shape[2]//pas_d_images, interval = 50, blit = False)
    
    # Save the result
    if save:
        writer = animation.FFMpegWriter(fps = 24, bitrate = 10000, codec = "libx264", extra_args = ["-pix_fmt", "yuv420p"])
        anim.save('file.mp4',writer=writer)
    
    return anim
    
##################################################

def anim_2D_flat(X, Y, L, pas_de_temps, pas_d_images, save = False,  myzlim = (-0.15, 0.15), vmin, vmax):
    """
    Fonction qui peremt d'annimer un représentation graphique en 2D z(x,y). Pour cela, on stock dans la liste L l'ensembles des images à affiicher. L = [ [z_1(x,y)], ... , [z_n(x,y)] ].
    On peut éventuellement enregistrer l'annimation.
    
    (X:np.ndarray (format 1D), Y:np.ndarray (format 1D), L:np.ndarray (format 3D), pas_de_temps:float, pas_d_images:int, save:bool) -> plot (+ .mp4)
    """
    
    
    fig, ax = plt.subplots(1, 1)
    plt.title("0.0 s", fontname = "serif", fontsize = 17)
    plt.xlabel("x [m]", fontname = "serif", fontsize = 12)
    plt.ylabel("y [m]", fontname = "serif", fontsize = 12)
    #Color map possible : RdBu_r
    mesh = plt.pcolormesh(X, Y,  L[:,:,0], vmin = myzlim[0] , vmax = myzlim[1], cmap = plt.cm.viridis)
    plt.colorbar(mesh, orientation = "vertical")

    # Update function
    def update_surf(num):
        ax.set_title("$u(x,y)$ à t = {} s".format(np.round(num*pas_d_images*pas_de_temps, 4)), fontname = "serif", fontsize = 16)
        mesh.set_array(L[:,:,pas_d_images*num][:-1, :-1].flatten())
        return mesh,

    anim = animation.FuncAnimation(fig, update_surf, frames = L.shape[2]//pas_d_images, interval = 10, blit = False)
    # Save the result
    if save:
        writer = animation.FFMpegWriter(fps = 24, bitrate = 10000, codec = "libx264", extra_args = ["-pix_fmt", "yuv420p"])
        anim.save('file.mp4',writer=writer)
        
    return anim
    
    
    
    
    
    
    
    
    
    







