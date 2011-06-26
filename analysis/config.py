#=============================================================
#	SIMULATIONS PARAMETERS
#=============================================================

#Files location
files = "../datas/Galaxy-"

#Number of Disk particles
N_D = 5000
#Number of Halo particles
N_H = 1000
#Number of Bulge particles
N_B = 1000
#Number of Total particles
N_T = N_D + N_H + N_B
#Time step
t_step = 1.
#T_max
t_max = 1000.
#Time Array
T = np.arange( 0, t_max, t_step )
#Number of time steps
NT = len(T)


#=============================================================
#	ANALYSIS PARAMETERS
#=============================================================

#Number of particles to analyse
N_a = 2000

#Reference Frecuency
Omega = 0.2

#Propierty to descompose
#N-Transient
N_Transient = 50.
#Time array
Tk = T[NT/N_Transient:NT]