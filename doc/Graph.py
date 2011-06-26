import pylab as plt
import numpy as np
import scipy.integrate as integ
import scipy as sp

Graphic = 1
#========================================================================================
# GRAPH 1 ( Kuzmin Potential )
#========================================================================================
if Graphic==1:
    
    def kuzmin_phi(x,y,z,a,Vo):
	xi = ( x**2 + y**2 + (a + abs(z))**2 )**(0.5)
	phi = Vo*np.log(xi)
	return phi
	
    def kuzmin_rho(x,y,z,a,Vo):
	xi = ( x**2 + y**2 + (a + abs(z))**2 )**(0.5)
	rho = Vo/(4*np.pi)*(  -1/(xi**2) + (2/xi**2)  )
	return rho
	
    #Parametros
    Nres = 200.
      #numero de contornos
    con = 15
    a = 10.
    Vo = 1.
    y = 0.0
    maxZ = 10
    maxX = 10
    
    #CONTOURS
    phi_dat = np.zeros( (Nres, Nres) )
    rho_dat = np.zeros( (Nres, Nres) )
    
    X = np.linspace( -maxX, maxX, Nres )
    Z = np.linspace( -maxZ, maxZ, Nres )
    i = 0
    for x in X:
	j = 0
	for z in Z:
	    phi_dat[j,i] = kuzmin_phi( x, y, z, a, Vo )
	    rho_dat[j,i] = kuzmin_rho( x, y, z, a, Vo )
	    j+=1
	i+=1
    
    plt.subplot(121)
    plt.contourf(X, Z, phi_dat, con)
    #plt.colorbar()
    plt.legend()
    plt.xlabel('x [Kpc]')
    plt.ylabel('z [Kpc]')
    plt.title('Contours of Kuzmin Potential $\\phi$')
    
    plt.subplot(122)
    plt.contourf(X, Z, rho_dat, con)
    #plt.colorbar()
    plt.legend()
    plt.xlabel('x [Kpc]')
    plt.ylabel('z [Kpc]')
    plt.title('Contours of Kuzmin Density $\\rho$')
    
    plt.show()