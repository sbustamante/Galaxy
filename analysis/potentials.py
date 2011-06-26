import pylab as plt
import numpy as np
import scipy.integrate as integ
import scipy as sp
import scipy.special as spc

#GENERAL PARAMETERS----------------------------------------------------------------------
G = 1

def sech(x):
    return 2*np.exp(x)/( np.exp(2*x) + 1 )

#========================================================================================
# MODEL
#========================================================================================
#BULGE:	PLUMMER SPHERE-------------------------------------------------------------------
  #PARAMETERS
    #scale-lenght of Bulge
eps_s = 0.4 	#[Manos, Athanassoula]
    #Bulge Mass
M_S = 0.01 	#[Manos, Athanassoula]

  #POTENTIAL
def bulge_plummer_p( x, y, z ):
      V = -G*M_S/np.sqrt(  x**2 + y**2 + z**2 + eps_s**2  )
      return V
      
  #DENSITY
def bulge_plummer_d( x, y, z ):
      rho = 3*M_S/(4*np.pi)*eps_s**2/(  x**2 + y**2 + z**2 + eps_s**2  )**(2.5)
      return rho


#DISK:	MIYAMOTO DISK-----------------------------------------------------------------
  #PARAMETERS
    #scale height of disk
z_0 = 0.1 	#[Ceverino, Klypin]
    #exponential lentgh
R_D = 1.5 	#[Ceverino, Klypin]
    #Disk Mass
M_D = 0.05 	#[Ceverino, Klypin]

  #POTENTIAL
def disk_miyamoto_p( x, y, z ):
      rho = -G*M_D/( x**2 + y**2 + (R_D + ( z**2 + z_0**2 )**0.5 )**2 )**0.5
      return rho*(1+0.2*(x**2-y**2)/(x**2+y**2))
      
  #DENSITY
def disk_miyamoto_d( x, y, z ):
      A = R_D
      B = z_0
      M = M_D
      rho = (3*M*(2*A*B**2*z**2 + 2*A*z**4 + B**2*x**2*np.sqrt(B**2 + z**2) + \
      B**2*y**2*np.sqrt(B**2 + z**2) + A**2*z**2*np.sqrt(B**2 + z**2) + \
      B**2*z**2*np.sqrt(B**2 + z**2) + x**2*z**2*np.sqrt(B**2 + z**2) + \
      y**2*z**2*np.sqrt(B**2 + z**2) + z**4*np.sqrt(B**2 + z**2))/(4*np.pi*(B**2 + \
      z**2)**(3/2.)*(x**2 + y**2 + (A + np.sqrt(B**2 + z**2))**2)**(5/2.)) + \
      (M*(A*B**2 + 3*(B**2 + z**2)*(3/2.)))/(4 *np.pi*(B**2 + z**2)**(3/2.)*(x**2 + y**2 + \
      (A + np.sqrt(B**2 + z**2))**2)**(3/2.)))
      return rho*(1+0.2*(x**2-y**2)/(x**2+y**2))
      
      
#HALO:	NFW PROFILE----------------------------------------------------------------------
  #PARAMETERS
    #Concentration parameter
C = 14. 	#[Ceverino, Klypin]
    #scale radius
R_H = 1.5 	#[Ceverino, Klypin]
    #Halo Mass
M_H = 0.94	#[Ceverino, Klypin]
    #Virial radius
R_Vr = C*R_D
    #caracteristic density
rho_0 = M_H/( 4*np.pi*R_H**3*( np.log(1+C) - C/(1+C) ) )

  #DENSITY
def halo_NFW_d( x, y, z ):
      x = np.sqrt( x**2+ y**2 + z**2 )/R_H
      rho = rho_0/( x*(1+x)**2 )
      return rho
      
  #POTENTIAL
def halo_NFW_p( x, y, z ):
      x = np.sqrt( x**2+ y**2 + z**2 )/R_H
      rho = -4*np.pi*rho_0*G*R_H**2*np.log(1+x)/x
      return rho


#========================================================================================
# CONTOURS PLOTS
#========================================================================================
#Contornos
Nres = 200.
#Tamano Grid
Xg = 3.

V_dat = np.zeros( (Nres, Nres) )
rho_dat = np.zeros( (Nres, Nres) )
X = np.linspace( -Xg, Xg, Nres )
Z = np.linspace( -Xg, Xg, Nres )
i = 0
for x in X:
    j = 0
    for z in Z:
	#Bulge
	#V_dat[i,j] = bulge_plummer_p( x, 0, z )
	#rho_dat[i,j] = bulge_plummer_d( x, 0, z )
	
	#Disk
	V_dat[j,i] = disk_miyamoto_p( x, z, 0 )
	rho_dat[j,i] = disk_miyamoto_d( x, z, 0 )
	
	#Halo
	#V_dat[i,j] = halo_NFW_p( x, 0, z )
	#rho_dat[j,i] = np.log10(halo_NFW_d( x, 0, z ))
	
	#Total
	#V_dat[j,i] =  disk_miyamoto_p( x, 0, z ) + bulge_plummer_p( x, 0, z ) + halo_NFW_p( x, 0, z ) 
	#rho_dat[j,i] = np.log10( disk_miyamoto_d( x, 0, z ) + bulge_plummer_d( x, 0, z ) + halo_NFW_d( x, 0, z ) )
	j+=1
    i+=1

plt.subplot(121)
plt.contourf( X, Z, V_dat, 20 )
plt.legend()
plt.xlabel('x [SU]')
plt.ylabel('y [SU]')

#Bulge =====
#plt.title('Bulge Model: Plummer sphere potential')
#Disk ======
#plt.title('Disk Model: Miyamoto disk potential')
#Halo ======
#plt.title('Halo Model: NFW potential')
#Total =====
#plt.title('Total Model potential')
#Bar =====
plt.title('Bar Model potential')

plt.subplot(122)
plt.contourf( X, Z, rho_dat, 20 )
#plt.matshow( rho_dat )
plt.legend()
plt.xlabel('x [SU]')
plt.ylabel('y [SU]')

#Bulge =====
#plt.title('Bulge Model: Plummer sphere density')
#Disk ======
#plt.title('Disk Model: Miyamoto disk density')
#Halo ======
#plt.title('Halo Model: NFW density  [$\log \\rho_H$]')
#Total =====
#plt.title('Total Model density'
#Bar =====
plt.title('Bar Model density')

#x = np.linspace(0,1.0, 100)
#plt.plot( x, bulge_plummer_d( x, 0, 0 ) )

plt.show()