#include <allvars.h>

/*****************************************************************
 NOMBRE:     dynamic_function
 FUNCION:    integra un paso del sistema de N particulas
 ARGUMENTOS: arreglo 'cluster' con estado de particulas, paso de
	     de integracion,y entero con especificacion de metodo
	     de integracion a utilizar.
	     0 --- Runge Kutta 4
	     1 --- Euler primer orden
 RETORNA:    0
 LOG:        27/01/2011    Creacion
*****************************************************************/
int dynamic_function( struct cluster galaxy[1],
		      double t,
		      struct intg_step f[MAXNPARTICLES] )
{
    int i, j;
    
    double phi[4];
    double phi_H[4];
    double phi_tmp[4];

    for( i=0; i<galaxy[0].N; i++ ){
	//Positions derivatives
	f[i].r[X] = galaxy[0].prts[i].v[X];
	f[i].r[Y] = galaxy[0].prts[i].v[Y];
	f[i].r[Z] = galaxy[0].prts[i].v[Z];
	
	//Velocities derivatives
	for( j=0; j<4; j++ ){
	    phi_H[j] = 0;; }
	    
	//POTENTIAL: Homogeneos Disk Component
	switch( galaxy[0].distros[HOM_DIS] ){
	    case 0:{
		none_disk( galaxy[0].prts[i].r, galaxy[0].prts[i].v, t,
			   galaxy[0].pot_conf, phi_tmp );
		break;}
	    case 1:{
		log_kuzmin_disk_phi( galaxy[0].prts[i].r, galaxy[0].prts[i].v, t,
				     galaxy[0].pot_conf, phi_tmp );
		break;}
	    case 2:{
		miyamoto_disk_phi( galaxy[0].prts[i].r, galaxy[0].prts[i].v, t,
				     galaxy[0].pot_conf, phi_tmp );
		break;}}
	for( j=0; j<4; j++ ){
	    phi_H[j] += phi_tmp[j]; }
	    
	//POTENTIAL: Bulge Component
	switch( galaxy[0].distros[BULGE] ){
	    case 0:{
		none_bulge( galaxy[0].prts[i].r, galaxy[0].prts[i].v, t,
			   galaxy[0].pot_conf, phi_tmp );
		break;}
	    case 1:{
		plummer_bulge_phi( galaxy[0].prts[i].r, galaxy[0].prts[i].v, t,
				   galaxy[0].pot_conf, phi_tmp );
		break;}}
	for( j=0; j<4; j++ ){
	    phi_H[j] += phi_tmp[j]; }
	    
	//POTENTIAL: Halo Component
	switch( galaxy[0].distros[HALO_DK] ){
	    case 0:{
		none_halo( galaxy[0].prts[i].r, galaxy[0].prts[i].v, t,
			   galaxy[0].pot_conf, phi_tmp );
		break;}
	    case 1:{
		NFW_phi( galaxy[0].prts[i].r, galaxy[0].prts[i].v, t,
			 galaxy[0].pot_conf, phi_tmp );}
	    case 2:{
		plummer_halo_phi( galaxy[0].prts[i].r, galaxy[0].prts[i].v, t,
				  galaxy[0].pot_conf, phi_tmp );
		break;}}
	for( j=0; j<4; j++ ){
	    phi_H[j] += phi_tmp[j]; }
	    
	//POTENTIAL: bar Component
	switch( galaxy[0].distros[BAR_DIS] ){
	    case 1:{
		miyamoto_bar_phi( galaxy[0].prts[i].r, galaxy[0].prts[i].v, t,
				  galaxy[0].pot_conf, phi_tmp );}}
	for( j=0; j<4; j++ ){
	    phi_H[j] += phi_tmp[j]; }
	    
	
	f[i].v[X] = -phi_H[X];
	f[i].v[Y] = -phi_H[Y];
	f[i].v[Z] = -phi_H[Z];
    }
	
    return 0;
}


/*****************************************************************
 NAME:       Total Potential
 FUNCTION:   calculate total potential
 INPUTS:     particle position and velocity, time, parameters
 RETURN:     0
 LOG:        28/05/2011    Creation
*****************************************************************/
double total_phi( double r[3],
		  double v[3],	  
		  double prmt[],
		  int distros[])
{
    double phi = 0;
    
    	//POTENTIAL: Homogeneos Disk Component
	switch( distros[HOM_DIS] ){
	    case 1:{
		phi += log_kuzmin_disk_phi(  r, v, prmt );
		break;}
	    case 2:{
		phi += miyamoto_disk_phi0( r, v, prmt );
		break;}}
	    
	//POTENTIAL: Bulge Component
	switch( distros[BULGE] ){
	    case 1:{
		phi += plummer_bulge_phi0( r, v, prmt );
		break;}}
	    
	//POTENTIAL: Halo Component
	switch( distros[HALO_DK] ){
	    case 1:{
		phi += NFW_phi0( r, v, prmt );
		break;}
	    case 2:{
		phi += plummer_halo_phi0( r, v, prmt );
		break;}}
    return phi;
}


//================================================================
//			HOMOGENEUS DISK COMPONENTS
//================================================================

/*****************************************************************
 NAME:       none_disk
 FUNCTION:   return None
 INPUTS:     particle position and velocity, time, parameters
 RETURN:     0
 LOG:        28/05/2011    Creation
*****************************************************************/
int none_disk( double r[3],
	       double v[3],
	       double t,
	       double prmt[],
	       double phi_H[4])
{
    phi_H[0] = 0;
    phi_H[1] = 0;
    phi_H[2] = 0;
    phi_H[3] = 0;
    
    return 0;
}

/*****************************************************************
 NAME:       log_kuzmin_disk_phi
 FUNCTION:   calculate potnetial of logarithmic Kuzmin Disk
 INPUTS:     particle position and velocity, time, parameters
 RETURN:     0
 LOG:        28/05/2011    Creation
*****************************************************************/
int log_kuzmin_disk_phi( double r[3],
			 double v[3],
			 double t,
			 double prmt[],
			 double phi_H[4] )
{
    double V0, a;
    double xi;
    
    V0 = prmt[HD_PAR+0];
    a = prmt[HD_PAR+1];
    
    //Xi coordinate
    xi = sqrt( r[X]*r[X] + r[Y]*r[Y] +
    ( a + fabs(r[Z])*( a + fabs(r[Z]) ) ) );
    
    //potential
    phi_H[3] = V0*log( xi );
    //X derivative of potential
    phi_H[X] = V0*r[X]/(xi*xi);
    //Y derivative of potential
    phi_H[Y] = V0*r[Y]/(xi*xi);
    //Z derivative of potential
    phi_H[Z] = V0*( a + fabs(r[Z]) )/(xi*xi)*fabs(r[Z])/r[Z];
    
    return 0;
}

/*****************************************************************
 NAME:       miyamoto_disk_phi
 FUNCTION:   calculate potential of Miyamoto Disk
 INPUTS:     particle position and velocity, time, parameters
 RETURN:     0
 LOG:        28/05/2011    Creation
*****************************************************************/
int miyamoto_disk_phi( double r[3],
		       double v[3],
		       double t,
		       double prmt[],
		       double phi_H[4] )
{
    double M_D, R_D, z0;
    double cc;
    
    //Disk Mass
    M_D = prmt[HD_PAR+0];
    //Caracteristic lenght disk
    R_D = prmt[HD_PAR+1];
    //Caracteristic height disk
    z0 = prmt[HD_PAR+2];
    
    cc = r[X]*r[X] + r[Y]*r[Y] + pow( R_D + sqrt( r[Z]*r[Z] + z0*z0 ), 2);
    
    //potential
    phi_H[3] = -GC*M_D/sqrt( cc );
    //X derivative of potential
    phi_H[X] = GC*M_D*r[X]/pow( cc, 1.5 );
    //Y derivative of potential
    phi_H[Y] = GC*M_D*r[Y]/pow( cc, 1.5 );
    //Z derivative of potential
    phi_H[Z] = ( R_D*GC*M_D*r[Z]/sqrt( r[Z]*r[Z] + z0*z0 ) + GC*M_D*r[Z] )/pow( cc, 1.5 );
    
    return 0;
}


/*****************************************************************
 NAME:       miyamoto_disk_phi0
 FUNCTION:   calculate potential of Miyamoto Disk
 INPUTS:     particle position and velocity, time, parameters
 RETURN:     0
 LOG:        28/05/2011    Creation
*****************************************************************/
double miyamoto_disk_phi0( double r[3],
			   double v[3],	  
			   double prmt[] )
{
    double M_D, R_D, z0;
    double cc;
    double phi;
    
    //Disk Mass
    M_D = prmt[HD_PAR+0];
    //Caracteristic lenght disk
    R_D = prmt[HD_PAR+1];
    //Caracteristic height disk
    z0 = prmt[HD_PAR+2];
    
    cc = r[X]*r[X] + r[Y]*r[Y] + pow( R_D + sqrt( r[Z]*r[Z] + z0*z0 ), 2);
    
    //potential
    phi = -GC*M_D/sqrt( cc );
    
    return phi;
}

/*****************************************************************
 NAME:       miyamoto_disk_rho
 FUNCTION:   calculate density of Miyamoto Disk
 INPUTS:     particle position and velocity, time, parameters
 RETURN:     0
 LOG:        28/05/2011    Creation
*****************************************************************/
double miyamoto_disk_rho( double r[3],
			  double v[3],	  
			  double prmt[] )
{
    double M_D, R_D, z0;
    double cc;
    double rho;
    double x, y, z;
    
    //Disk Mass
    M_D = prmt[HD_PAR+0];
    //Caracteristic lenght disk
    R_D = prmt[HD_PAR+1];
    //Caracteristic height disk
    z0 = prmt[HD_PAR+2];
    
    x = r[X]; y = r[Y]; z = r[Z];
    
    //density
    rho = -((3*GC*x*x*M_D)/pow(x*x+y*y+pow(sqrt(z*z+z0*z0)+R_D,2),2.5))-\
    (3*GC*y*y*M_D)/pow(x*x+y*y+pow(sqrt(z*z+z0*z0)+R_D,2),2.5)-\
    (3*pow(z,4)*M_D)/(4*M_PI*(z*z+z0*z0)*pow(x*x+y*y+pow(sqrt(z*z+z0*z0)+R_D,2),2.5))-\
    (3*z*z*z0*z0*M_D)/(4*M_PI*(z*z+z0*z0)*pow(x*x+y*y+pow(sqrt(z*z+z0*z0)+R_D,2),2.5))-\
    (3*z*z*M_D*R_D)/(2*M_PI*sqrt(z*z+z0*z0)*pow(x*x+y*y+pow(sqrt(z*z+z0*z0)+R_D,2),2.5))-\
    (3*z*z*M_D*R_D*R_D)/(4*M_PI*(z*z+z0*z0)*pow(x*x+y*y+pow(sqrt(z*z+z0*z0)+R_D,2),2.5))+\
    (2*GC*M_D)/pow(x*x+y*y+pow(sqrt(z*z+z0*z0)+R_D,2),1.5)+\
    M_D/(4*M_PI*pow(x*x+y*y+pow(sqrt(z*z+z0*z0)+R_D,2),1.5))-\
    (z*z*M_D*R_D)/(4*M_PI*pow(z*z+z0*z0,1.5)*pow(x*x+y*y+pow(sqrt(z*z+z0*z0)+R_D,2),1.5))+\
    (M_D*R_D)/(4*M_PI*sqrt(z*z+z0*z0)*pow(x*x+y*y+pow(sqrt(z*z+z0*z0)+R_D,2),1.5));
    
    
    return rho;
}

/*****************************************************************
 NAME:       miyamoto_disk_rho2
 FUNCTION:   calculate approximatte density of Miyamoto Disk
 INPUTS:     particle position and velocity, time, parameters
 RETURN:     0
 LOG:        28/05/2011    Creation
*****************************************************************/
double miyamoto_disk_rho2( double r[3],
			  double v[3],	  
			  double prmt[] )
{
    double M_D, R_D, z0, rho0;
    double R[3];
    double rho;
    
    //Disk Mass
    M_D = prmt[HD_PAR+0];
    //Caracteristic lenght disk
    R_D = prmt[HD_PAR+1];
    //Caracteristic height disk
    z0 = prmt[HD_PAR+2];
    
    //density
    R[X] = 0; R[Y] = 0; R[Z] = 0;
    rho0 = miyamoto_disk_rho( R, R, prmt );
    
    //Disk Mass
    M_D = prmt[HD_PAR+0];
    //Caracteristic lenght disk
    R_D = prmt[HD_PAR+1];
    //Caracteristic height disk
    z0 = prmt[HD_PAR+2];
    
    //Inverse
    rho = rho0/( 1 + pow( 1.5*( sqrt(r[X]*r[X]+r[Y]*r[Y]) )/R_D, 3. ) );

    return rho;
}

/*****************************************************************
 NAME:       miyamoto_disk_inv
 FUNCTION:   calculate position for a density given
 INPUTS:     density, parameters
 RETURN:     density
 LOG:        15/06/2011    Creation
*****************************************************************/
double miyamoto_disk_inv( double rho,
			  double prmt[])
{
    double M_D, R_D, z0, rho0;
    double cc;
    double R[3];
    double r;
    
    R[X] = 0; R[Y] = 0; R[Z] = 0;
    rho0 = miyamoto_disk_rho( R, R, prmt );
    
    //Disk Mass
    M_D = prmt[HD_PAR+0];
    //Caracteristic lenght disk
    R_D = prmt[HD_PAR+1];
    //Caracteristic height disk
    z0 = prmt[HD_PAR+2];
    
    //Inverse
    r = 1.5*pow( rho0/rho - 1, 1/3. );
   
    return r;
}

/*****************************************************************
 NAME:       miyamoto_disk_DF
 FUNCTION:   calculate DF
 INPUTS:     density, parameters
 RETURN:     density
 LOG:        15/06/2011    Creation
*****************************************************************/
double miyamoto_disk_DF( double E,
			 double Lz,
			 double prmt[])
{
    double M_D, R_D, z0, sigma, R_off;
    double q, vc, Sigma0, F;
    double df = 0;

    //Disk Mass
    M_D = prmt[HD_PAR+0];
    //Caracteristic lenght disk
    R_D = prmt[HD_PAR+1];
    //Caracteristic height disk
    z0 = prmt[HD_PAR+2];
    //Off radius
    sigma = prmt[HD_PAR+3];
    //Off radius
    R_off = prmt[HD_PAR+4];
    
    //Circular velocity
    vc = sqrt( GC*M_D/R_off );
    //Dispersion index
    q = vc*vc/(sigma*sigma) - 1;
    //Disk superficial density
    Sigma0 = M_D/(2*M_PI*R_D*R_off);
    //F index
    F = Sigma0*pow(vc,q)/(pow(2,q/2.)*sqrt(M_PI)*factorial( 0.5*q - 0.5 )*pow(sigma,q+2));
    
    //DF
    if( E<0 && Lz>0 ){
	df = F*pow(Lz/(R_D*vc), q)*exp(E/(sigma*sigma));}
    
    return df;
}

//================================================================
//			BULGE COMPONENTS
//================================================================

/*****************************************************************
 NAME:       none_bulge
 FUNCTION:   return None
 INPUTS:     particle position and velocity, time, parameters
 RETURN:     0
 LOG:        28/05/2011    Creation
*****************************************************************/
int none_bulge( double r[3],
	       double v[3],
	       double t,
	       double prmt[],
	       double phi_H[4])
{
    phi_H[0] = 0;
    phi_H[1] = 0;
    phi_H[2] = 0;
    phi_H[3] = 0;
    
    return 0;
}

/*****************************************************************
 NAME:       plummer_bulge_phi
 FUNCTION:   calculate potential of plummer sphere
 INPUTS:     particle position and velocity, time, parameters
 RETURN:     0
 LOG:        28/05/2011    Creation
*****************************************************************/
int plummer_bulge_phi( double r[3],
		       double v[3],
		       double t,
		       double prmt[],
		       double phi_H[4] )
{
    double M_B, R_B;
    double R;
    
    //Bulge Mass
    M_B = prmt[BL_PAR+0];
    //Bulge caracteristic lenght
    R_B = prmt[BL_PAR+1];
    
    //r coordinate
    R = sqrt( r[X]*r[X] + r[Y]*r[Y] + r[Z]*r[Z] );
    
    //potential
    phi_H[3] = -GC*M_B/pow( R*R + R_B*R_B, 0.5 );
    //X derivative of potential
    phi_H[X] = GC*M_B*r[X]/pow( R*R + R_B*R_B, 1.5 );
    //Y derivative of potential
    phi_H[Y] = GC*M_B*r[Y]/pow( R*R + R_B*R_B, 1.5 );
    //Z derivative of potential
    phi_H[Z] = GC*M_B*r[Z]/pow( R*R + R_B*R_B, 1.5 );
    
    return 0;
}

/*****************************************************************
 NAME:       plummer_bulge_phi0
 FUNCTION:   calculate potential of plummer sphere
 INPUTS:     particle position and velocity, time, parameters
 RETURN:     potential
 LOG:        28/05/2011    Creation
*****************************************************************/
double plummer_bulge_phi0( double r[3],
			   double v[3],
		           double prmt[] )
{
    double M_B, R_B;
    double R;
    
    //Bulge Mass
    M_B = prmt[BL_PAR+0];
    //Bulge caracteristic lenght
    R_B = prmt[BL_PAR+1];
    
    //r coordinate
    R = sqrt( r[X]*r[X] + r[Y]*r[Y] + r[Z]*r[Z] );
    
    //potential
    return  -GC*M_B/pow( R*R + R_B*R_B, 0.5 );
}

/*****************************************************************
 NAME:       plummer_bulge_rho
 FUNCTION:   calculate density of plummer sphere
 INPUTS:     positions, velocity, parameters
 RETURN:     density
 LOG:        28/05/2011    Creation
*****************************************************************/
double plummer_bulge_rho( double r[3],
			  double v[3],
			  double prmt[])
{
    double M_B, R_B;
    double rho;
    double R;
    
    //Bulge Mass
    M_B = prmt[BL_PAR+0];
    //Bulge caracteristic lenght
    R_B = prmt[BL_PAR+1];
    
    //r coordinate
    R = sqrt( r[X]*r[X] + r[Y]*r[Y] + r[Z]*r[Z] );
    
    //potential
    rho = 3*M_B/(4*M_PI)*R_B*R_B/pow( R*R + R_B*R_B, 2.5 );

    return rho;
}

/*****************************************************************
 NAME:       plummer_bulge_inv
 FUNCTION:   calculate position for a density
 INPUTS:     density, parameters
 RETURN:     density
 LOG:        28/05/2011    Creation
*****************************************************************/
double plummer_bulge_inv( double rho,
			  double prmt[])
{
    double M_B, R_B;
    double R;
    
    //Bulge Mass
    M_B = prmt[BL_PAR+0];
    //Bulge caracteristic lenght
    R_B = prmt[BL_PAR+1];
        
    //potential
    R = pow( pow( 4*M_PI*rho/( 3*M_B*R_B*R_B ), -2/5. ) - R_B*R_B ,0.5 );

    return R;
}

/*****************************************************************
 NAME:       plummer_bulge_DF
 FUNCTION:   calculate Distribution Function
 INPUTS:     Energy, parameters
 RETURN:     density
 LOG:        28/05/2011    Creation
*****************************************************************/
double plummer_bulge_DF( double E,
			 double prmt[])
{
    double M_B, R_B;
    double DF = 0;
    
    //Bulge Mass
    M_B = prmt[BL_PAR+0];
    //Bulge caracteristic lenght
    R_B = prmt[BL_PAR+1];
        
    //DF
    if( E<0 )
	DF = 96/( 7*sqrt(8)*M_PI*M_PI*M_PI )*R_B*R_B/( pow(GC,5)*pow(M_B,4) )*pow( -E, 3.5);

    return DF;
}

//================================================================
//			DARK HALO COMPONENTS
//================================================================

/*****************************************************************
 NAME:       none_halo
 FUNCTION:   return None
 INPUTS:     particle position and velocity, time, parameters
 RETURN:     0
 LOG:        28/05/2011    Creation
*****************************************************************/
int none_halo( double r[3],
	       double v[3],
	       double t,
	       double prmt[],
	       double phi_H[4])
{
    phi_H[0] = 0;
    phi_H[1] = 0;
    phi_H[2] = 0;
    phi_H[3] = 0;
    
    return 0;
}

/*****************************************************************
 NAME:       NFW_phi
 FUNCTION:   calculate potnetial of NFW
 INPUTS:     particle position and velocity, time, parameters
 RETURN:     0
 LOG:        28/05/2011    Creation
*****************************************************************/
int NFW_phi( double r[3],
	     double v[3],
	     double t,
	     double prmt[],
	     double phi_H[4] )
{
    double M_H, R_H, R_HV;
    double rho_s, C;
    double R, x;
    
    //Halo Mass
    M_H = prmt[DH_PAR+0];
    //Halo caracteristic lenght
    R_H = prmt[DH_PAR+1];
    //Virial radius
    R_HV = prmt[DH_PAR+4];
    
    //r coordinate
    R = sqrt( r[X]*r[X] + r[Y]*r[Y] + r[Z]*r[Z] );
    x = R/R_H;
    
    //Virial Parameter
    C = R_HV/R_H;
    // Caracteristic density
    rho_s = M_H/( 4*M_PI*pow(R_H,3)*( log(1+C) - C/(1+C) ) );
    
        
    //potential
    phi_H[3] = -rho_s/( x*pow(1+x,2) );
    //X derivative of potential
    phi_H[X] = (rho_s/R_H)*(r[X]*(3*x+1))/( pow(x,3)*pow(x+1,3) );
    //Y derivative of potential
    phi_H[Y] = (rho_s/R_H)*(r[Y]*(3*x+1))/( pow(x,3)*pow(x+1,3) );
    //Z derivative of potential
    phi_H[Z] = (rho_s/R_H)*(r[Z]*(3*x+1))/( pow(x,3)*pow(x+1,3) );
    
    return 0;
}

/*****************************************************************
 NAME:       NFW_phi0
 FUNCTION:   calculate potnetial of NFW
 INPUTS:     particle position and velocity, time, parameters
 RETURN:     0
 LOG:        16/06/2011    Creation
*****************************************************************/
double NFW_phi0( double r[3],
		 double v[3],
		 double prmt[] )
{
    double M_H, R_H, R_HV;
    double rho_s, C;
    double R, x;
    double phi;
    
    //Halo Mass
    M_H = prmt[DH_PAR+0];
    //Halo caracteristic lenght
    R_H = prmt[DH_PAR+1];
    //Virial radius
    R_HV = prmt[DH_PAR+4];
    
    //r coordinate
    R = sqrt( r[X]*r[X] + r[Y]*r[Y] + r[Z]*r[Z] );
    x = R/R_H;
    
    //Virial Parameter
    C = R_HV/R_H;
    // Caracteristic density
    rho_s = M_H/( 4*M_PI*pow(R_H,3)*( log(1+C) - C/(1+C) ) );
    
        
    //potential
    phi = rho_s/( x*pow(1+x,2) );
    
    return phi;
}


/*****************************************************************
 NAME:       plummer_halo_phi
 FUNCTION:   calculate potential of plummer sphere
 INPUTS:     particle position and velocity, time, parameters
 RETURN:     0
 LOG:        28/05/2011    Creation
*****************************************************************/
int plummer_halo_phi( double r[3],
		       double v[3],
		       double t,
		       double prmt[],
		       double phi_H[4] )
{
    double M_B, R_B;
    double R;
    
    //Bulge Mass
    M_B = prmt[DH_PAR+0];
    //Bulge caracteristic lenght
    R_B = prmt[DH_PAR+1];
    
    //r coordinate
    R = sqrt( r[X]*r[X] + r[Y]*r[Y] + r[Z]*r[Z] );
    
    //potential
    phi_H[3] = -GC*M_B/pow( R*R + R_B*R_B, 0.5 );
    //X derivative of potential
    phi_H[X] = GC*M_B*r[X]/pow( R*R + R_B*R_B, 1.5 );
    //Y derivative of potential
    phi_H[Y] = GC*M_B*r[Y]/pow( R*R + R_B*R_B, 1.5 );
    //Z derivative of potential
    phi_H[Z] = GC*M_B*r[Z]/pow( R*R + R_B*R_B, 1.5 );
    
    return 0;
}

/*****************************************************************
 NAME:       plummer_halo_phi0
 FUNCTION:   calculate potential of plummer sphere
 INPUTS:     particle position and velocity, time, parameters
 RETURN:     potential
 LOG:        28/05/2011    Creation
*****************************************************************/
double plummer_halo_phi0( double r[3],
			   double v[3],
		           double prmt[] )
{
    double M_B, R_B;
    double R;
    
    //Bulge Mass
    M_B = prmt[DH_PAR+0];
    //Bulge caracteristic lenght
    R_B = prmt[DH_PAR+1];
    
    //r coordinate
    R = sqrt( r[X]*r[X] + r[Y]*r[Y] + r[Z]*r[Z] );
    
    //potential
    return  -GC*M_B/pow( R*R + R_B*R_B, 0.5 );
}

/*****************************************************************
 NAME:       plummer_halo_rho
 FUNCTION:   calculate density of plummer sphere
 INPUTS:     positions, velocity, parameters
 RETURN:     density
 LOG:        28/05/2011    Creation
*****************************************************************/
double plummer_halo_rho( double r[3],
			  double v[3],
			  double prmt[])
{
    double M_B, R_B;
    double rho;
    double R;
    
    //Bulge Mass
    M_B = prmt[DH_PAR+0];
    //Bulge caracteristic lenght
    R_B = prmt[DH_PAR+1];
    
    //r coordinate
    R = sqrt( r[X]*r[X] + r[Y]*r[Y] + r[Z]*r[Z] );
    
    //potential
    rho = 3*M_B/(4*M_PI)*R_B*R_B/pow( R*R + R_B*R_B, 2.5 );

    return rho;
}

/*****************************************************************
 NAME:       plummer_halo_inv
 FUNCTION:   calculate position for a density
 INPUTS:     density, parameters
 RETURN:     density
 LOG:        28/05/2011    Creation
*****************************************************************/
double plummer_halo_inv( double rho,
			 double prmt[])
{
    double M_B, R_B;
    double R;
    
    //Bulge Mass
    M_B = prmt[DH_PAR+0];
    //Bulge caracteristic lenght
    R_B = prmt[DH_PAR+1];
        
    //potential
    R = pow( pow( 4*M_PI*rho/( 3*M_B*R_B*R_B ), -2/5. ) - R_B*R_B ,0.5 );

    return R;
}

/*****************************************************************
 NAME:       plummer_halo_DF
 FUNCTION:   calculate Distribution Function
 INPUTS:     Energy, parameters
 RETURN:     density
 LOG:        28/05/2011    Creation
*****************************************************************/
double plummer_halo_DF( double E,
			 double prmt[])
{
    double M_B, R_B;
    double DF = 0;
    
    //Bulge Mass
    M_B = prmt[DH_PAR+0];
    //Bulge caracteristic lenght
    R_B = prmt[DH_PAR+1];
        
    //DF
    if( E<0 )
	DF = 96/( 7*sqrt(8)*M_PI*M_PI*M_PI )*R_B*R_B/( pow(GC,5)*pow(M_B,4) )*pow( -E, 3.5);

    return DF;
}

//================================================================
//			BARS COMPONENTS
//================================================================

/*****************************************************************
 NAME:       miyamoto_bar_phi
 FUNCTION:   calculate potential of Miyamoto Bar Disk
 INPUTS:     particle position and velocity, time, parameters
 RETURN:     0
 LOG:        17/05/2011    Creation
*****************************************************************/
int miyamoto_bar_phi( double r[3],
		      double v[3],
		      double t,
		      double prmt[],
		      double phi_H[4] )
{
    double M_D, R_D, z0, OmegaB;
    double phi_part, pphi_part;
    double cc;
    
    //Bar Mass
    M_D = prmt[BM_PAR+0];
    //Caracteristic lenght bar disk
    R_D = prmt[BM_PAR+1];
    //Caracteristic height bar disk
    z0 = prmt[BM_PAR+2];
    //Rotation rate
    OmegaB = -prmt[BM_PAR+3];
    
    cc = r[X]*r[X] + r[Y]*r[Y] + pow( R_D + sqrt( r[Z]*r[Z] + z0*z0 ), 2);
    
    phi_part = ( cos( OmegaB*t )*(r[X]*r[X] - r[Y]*r[Y]) - 2*sin( OmegaB*t )*r[X]*r[Y] )/( r[X]*r[X] + r[Y]*r[Y] );
    
    //potential
    phi_H[3] = -GC*M_D/sqrt( cc );
    //X derivative of potential
    phi_H[X] = GC*M_D*r[X]/pow( cc, 1.5 )*phi_part - GC*M_D/sqrt( cc )*\
    ( 2*r[Y]*( 2*cos(OmegaB*t)*r[X]*r[Y] + sin(OmegaB*t)*(r[X]*r[X] - r[Y]*r[Y]))/pow( r[X]*r[X] + r[Y]*r[Y], 2 ) );
    //Y derivative of potential
    phi_H[Y] = GC*M_D*r[Y]/pow( cc, 1.5 )*phi_part + GC*M_D/sqrt( cc )*\
    ( 2*r[X]*( 2*cos(OmegaB*t)*r[X]*r[Y] + sin(OmegaB*t)*(r[X]*r[X] - r[Y]*r[Y]))/pow( r[X]*r[X] + r[Y]*r[Y], 2 ) );
    //Z derivative of potential
    phi_H[Z] = ( R_D*GC*M_D*r[Z]/sqrt( r[Z]*r[Z] + z0*z0 ) + GC*M_D*r[Z] )/pow( cc, 1.5 )*phi_part;
    
    return 0;
}