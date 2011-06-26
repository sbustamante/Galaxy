import pylab as plt
import numpy as np
import scipy.integrate as integ
import scipy as sp

execfile("config.py")

#Load Datas Function
def LoadDatas(  ):
    j = 0
    Datos = []
    for t in np.arange(0, t_max, t_step):
	filename = "%s%05d.0"%( files, j )
	print filename
	Datos.append( np.loadtxt( filename ) )
	j += 1
    return Datos;
    
#Frecuencies Interval
def fk(k):
    return k/(NN*t_max*1.0/NT)
    
#DFT
def Sk(k,F):
    Skh = 0
    for n in xrange(0,NN):
	Skh += F[n]*np.exp( 2*np.pi*1j*k*(n+1)/NN )
    return abs(Skh*1/(NN)**(0.5))**2

#========================================================================================
# MODEL
#========================================================================================
#Loading datas
Datos = LoadDatas()

pr = []
Frec = np.zeros( (N_a,2) )
for i in xrange(0, N_a):
    #Random Particle
    pr.append( np.random.randint( N_D ) )
    
    #Load Positions of each random particle
    X = [1]; Y = [0]; Z = [0];
    VX = [1]; VY = [0]; VZ = [0];
    j = 0
    th = 0
    for t in np.arange(0, t_max, t_step):
	X.append( Datos[j][pr[i]][1] )
	Y.append( Datos[j][pr[i]][2] )
	Z.append( Datos[j][pr[i]][3] )
	VX.append( Datos[j][pr[i]][4] )
	VY.append( Datos[j][pr[i]][5] )
	VZ.append( Datos[j][pr[i]][6] )
	j += 1
	
    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)
    R = np.sqrt(X**2 + Y**2)
    THp = 1/R**2*( X*VY - Y*VX )
    
    #POWER SPECTRUM---------------------------------------------
    F_R = R[NT/N_Transient:NT] - sum(R)/len(R)
    NN = len(F_R)
    #K array
    k = np.linspace(0,NN/2.,NN/2.)
    #Frecuencies array
    f_k = fk(k)
    S_kR = Sk(k,F_R)[1:-2]
    Frec[i][0] = f_k[ np.argsort( S_kR )[-1] ]
    
    F_TH = THp[NT/N_Transient:NT] - min(THp)
    NN = len(F_TH)
    #K array
    k = np.linspace(0,NN/2.,NN/2.)
    #Frecuencies array
    f_k = fk(k)
    S_kT = Sk(k,F_TH)[1:-2]
    Frec[i][1] = f_k[ np.argsort( S_kT )[-1] ]


plt.title("Frequencies for Bar Galaxy")
plt.xlabel( "Angular Frequency  $\Omega-\Omega_B$ $[TS^{-1}]$" )
plt.ylabel( "Radial Frequency $\kappa$ $[TS^{-1}]$" )
plt.plot( f_k-Omega, (f_k) , color='red');
plt.plot( f_k-Omega, (f_k)/2, color='red' );
plt.plot( f_k-Omega, (f_k)/3, color='red' );
plt.plot(Frec[:,1]-Omega, Frec[:,0], '.', color='black', label='Orbital Frequencies' )
plt.ylim( (0.0,0.05) )
plt.xlim( (0-Omega,0.15-Omega) )
plt.legend()
plt.show()

