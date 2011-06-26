#include <allvars.h>

/*****************************************************************
 NAME:       distance
 FUNCTION:   calculate distance between two vectors
 INPUTS:     two 3D arrays
 RETURN:     0
 LOG:        27/01/2011    Creation
*****************************************************************/
double distance( double r1[3],
		 double r2[3] )
{    
    return pow( (r1[X]-r2[X])*(r1[X]-r2[X]) +
		(r1[Y]-r2[Y])*(r1[Y]-r2[Y]) +
		(r1[Z]-r2[Z])*(r1[Z]-r2[Z]), 0.5 );
}

/*****************************************************************
 NAME:       factorial
 FUNCTION:   calculate factorial of number
 INPUTS:     two 3D arrays
 RETURN:     0
 LOG:        27/01/2011    Creation
*****************************************************************/
int factorial( int n )
{    
    int N=1, i;
    for( i=1; i<=n; i++ )
	N = N*i;
    return N;
}

/*****************************************************************
 NAME:       integration_step
 FUNCTION:   integrate a step of system of N particles
 INPUTS:     'cluster' structure with particles states, step of 
	     integration and integer with method to use.
	     0 --- Four Order Runge Kutta
	     1 --- First Order Euler
	     2 --- Simplectic First Order Euler
	     3 --- Simplectic Second LeapFrog
	     4 --- Simplectic Four Order
 RETORNA:    0
 LOG:        27/01/2011    Creation
*****************************************************************/
int integration_step( struct cluster galaxy[1], 
		      double t,
		      double h, 
		      int method )
{    
  
  //RK4 de primer orden
    if( method==0 )
	rk4_step( galaxy, t, h );
    
    //Euler de primer orden
    if( method==1 )
 	euler_step( galaxy, t, h );

    return 0;
}

/*****************************************************************
 NAME:       euler_step
 FUNCTION:   euler method integrator
 INPUTS:     'cluster' structure with particles state, step of 
	     integration.
 RETURN:     0
 LOG:        27/01/2011    Creation
*****************************************************************/
int euler_step( struct cluster galaxy[1],
		double t,
		double h )
{
    int i;
    struct intg_step f[MAXNPARTICLES];

    //funcion dinamica
    dynamic_function( galaxy, t, f );

    for( i=0; i<galaxy[0].N; i++ ){
	galaxy[0].prts[i].r[X] = galaxy[0].prts[i].r[X] + h*f[i].r[X];
	galaxy[0].prts[i].r[Y] = galaxy[0].prts[i].r[Y] + h*f[i].r[Y];
	galaxy[0].prts[i].r[Z] = galaxy[0].prts[i].r[Z] + h*f[i].r[Z];
	
	galaxy[0].prts[i].v[X] = galaxy[0].prts[i].v[X] + h*f[i].v[X];
	galaxy[0].prts[i].v[Y] = galaxy[0].prts[i].v[Y] + h*f[i].v[Y];
	galaxy[0].prts[i].v[Z] = galaxy[0].prts[i].v[Z] + h*f[i].v[Z];}
    
    return 0;
}


/*****************************************************************
 NAME:       rk4_step
 FUNCTION:   rk4 method integrator
 INPUTS:     'cluster' structure with particles state, step of 
	     integration.
 RETURN:     0
 LOG:        27/01/2011    Creation
*****************************************************************/
int rk4_step( struct cluster galaxy[1],
	      double t,
	      double h )
{
    int i, k;
    struct cluster glx_h[1];
    struct intg_step k1[MAXNPARTICLES], k2[MAXNPARTICLES], k3[MAXNPARTICLES], k4[MAXNPARTICLES];

    //Inicializacion
    for( i=0; i<galaxy[0].N; i++ ){
	for( k=0; k<3; k++ ){
	    glx_h[0].prts[i].r[k] = galaxy[0].prts[i].r[k];
	    glx_h[0].prts[i].v[k] = galaxy[0].prts[i].v[k];}}
    
    //Calculo k1
    dynamic_function( galaxy, t, k1 );
    
    for( i=0; i<galaxy[0].N; i++ ){
	for( k=0; k<3; k++ ){
	    galaxy[0].prts[i].r[k] = glx_h[0].prts[i].r[k] + h*k1[i].r[k]/2.0;
	    galaxy[0].prts[i].v[k] = glx_h[0].prts[i].v[k] + h*k1[i].v[k]/2.0;}}
	    
    //Calculo k2
    dynamic_function( galaxy, t + h/2.0, k2 );
    
    for( i=0; i<galaxy[0].N; i++ ){
      	for( k=0; k<3; k++ ){
	    galaxy[0].prts[i].r[k] = glx_h[0].prts[i].r[k] + h*k2[i].r[k]/2.0;
	    galaxy[0].prts[i].v[k] = glx_h[0].prts[i].v[k] + h*k2[i].v[k]/2.0;}}

    //Calculo k3
    dynamic_function( galaxy, t + h/2.0, k3 );
    
    for( i=0; i<galaxy[0].N; i++ ){
	for( k=0; k<3; k++ ){
	    galaxy[0].prts[i].r[k] = glx_h[0].prts[i].r[k] + h*k3[i].r[k];
	    galaxy[0].prts[i].v[k] = glx_h[0].prts[i].v[k] + h*k3[i].v[k];}}
    
    //Calculo k4
    dynamic_function( galaxy, t + h, k4 );

    //Calculo del siguiente paso
    for( i=0; i<galaxy[0].N; i++ ){
	for( k=0; k<3; k++ ){
	    galaxy[0].prts[i].r[k] = glx_h[0].prts[i].r[k] + h*( k1[i].r[k] + 2*k2[i].r[k] + 2*k3[i].r[k] + k4[i].r[k] )/6.;
	    galaxy[0].prts[i].v[k] = glx_h[0].prts[i].v[k] + h*( k1[i].v[k] + 2*k2[i].v[k] + 2*k3[i].v[k] + k4[i].v[k] )/6.;}
	    //Total energy
	    galaxy[0].prts[i].E = 0.5*( pow(galaxy[0].prts[i].v[X],2)+pow(galaxy[0].prts[i].v[Y],2)+\
	    pow(galaxy[0].prts[i].v[Z],2) ) + total_phi( galaxy[0].prts[i].r, galaxy[0].prts[i].v, galaxy[0].pot_conf, galaxy[0].distros );
	    //Total Angular momentum z
	    galaxy[0].prts[i].Lz = (galaxy[0].prts[i].r[X]*galaxy[0].prts[i].v[Y] - galaxy[0].prts[i].r[Y]*galaxy[0].prts[i].v[X] );}

    return 0;
}