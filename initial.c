#include <allvars.h>

/********************************************************************
 NOMBRE:     initial_conditions
 FUNCION:    genera condiciones iniciales para distribucion en un 
	     archivo de salida.
 ARGUMENTOS: estructura 'cluster' para inicializar las particulas
	     en la distribucion, archivo de parametros.Nombre de 
	     archivo de salida
 RETORNA:    0							    
 LOG:        21/01/2011    Creacion                                
********************************************************************/
int initial_conditions( struct cluster galaxy[1], 
			double prmt[],
			char filename[] )
{  
    int i;
    char distribution[10];
    FILE *file;
    
    //Uniform
//     uniform_cube( galaxy, prmt );
//     sprintf( distribution, "uniform_cube" );
    
    //Acording potential model
    initial_DF( galaxy, prmt );
    sprintf( distribution, "initial_DF" );

    //File printing
    file = fopen( filename, "w" );
    fprintf( file, "#Label\tX\t\tY\t\tZ\t\tVx\t\tVy\t\tVz\t\tEnergy\t\tLz\t\tMass\t\ttype\n" );
    for( i=0; i<galaxy[0].N; i++ )
	fprintf( file,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n", 
		 galaxy[0].prts[i].label, 
		 galaxy[0].prts[i].r[0], galaxy[0].prts[i].r[1], galaxy[0].prts[i].r[2],
		 galaxy[0].prts[i].v[0], galaxy[0].prts[i].v[1], galaxy[0].prts[i].v[2], 
		 galaxy[0].prts[i].E, galaxy[0].prts[i].Lz, galaxy[0].prts[i].m, 
		 galaxy[0].prts[i].type);
	  
    fclose(file);
    
    printf( "  * Condiciones iniciales generadas!\n\tDistribución: %s\n\tArchivo de salida: '%s'\n", 
	    distribution, filename );
  
    return 0;
}


/********************************************************************
 NOMBRE:     uniform_cube                                           
 FUNCION:    genera una distribucion uniforme sobre un cubo de     
             lados dados                                           
 ARGUMENTOS: estructura 'cluster' para retorno de posicion y      
	     velocidad de cada particula y arreglo con parametros 
 COORDENADAS:							    
	     *Coord 1 -- Lado x
	     *Coord 2 -- Lado y
	     *Coord 3 -- Lado z
	     *Veloc 1 -- Magnitud maxima velocidad
	     *Veloc 2 -- NULL
	     *Veloc 3 -- NULL                 			    
 RETORNA:    0							    
 LOG:        15/09/2009    Creacion           
	     26/01/2011	   Implementacion estructura cluster
********************************************************************/
int uniform_cube( struct cluster galaxy[1], 
		  double prmt[NMAX2] )
{
    int i,j;

    srand48(time(NULL));
    //Positions an velocities for disk stars
    for( i=0; i<galaxy[0].N1; i++ ){
	//Mass
	galaxy[0].prts[i].m = prmt[MASS];
	
        //Position
	galaxy[0].prts[i].r[X] = prmt[COOR_1]*(1-2*drand48());
	galaxy[0].prts[i].r[Y] = prmt[COOR_2]*(1-2*drand48());
	galaxy[0].prts[i].r[Z] = prmt[COOR_3]*(1-2*drand48())/10.0;

	//Velocity
	galaxy[0].prts[i].v[X] = (1-2*drand48())*prmt[VEL_1];
	galaxy[0].prts[i].v[Y] = (1-2*drand48())*prmt[VEL_2];
	galaxy[0].prts[i].v[Z] = 0;
	
	//Label
	galaxy[0].prts[i].label = i+1;
	galaxy[0].prts[i].type = 0;}

    //Positions an velocities for halo and bulge stars
    for( i=galaxy[0].N1; i<galaxy[0].N; i++ ){
	//Mass
	galaxy[0].prts[i].m = prmt[MASS];
	
        //Position
	galaxy[0].prts[i].r[X] = prmt[COOR_1]*(1-2*drand48());
	galaxy[0].prts[i].r[Y] = prmt[COOR_2]*(1-2*drand48());
	galaxy[0].prts[i].r[Z] = prmt[COOR_3]*(1-2*drand48());

	//Velocity
	galaxy[0].prts[i].v[X] = (1-2*drand48())*prmt[VEL_1];
	galaxy[0].prts[i].v[Y] = (1-2*drand48())*prmt[VEL_2];
	galaxy[0].prts[i].v[Z] = (1-2*drand48())*prmt[VEL_3];
	
	//Label
	galaxy[0].prts[i].label = i;
	galaxy[0].prts[i].type = 1;}
	    
    return 0;
}


/********************************************************************
 NAME:       initial_DF
 FUNCTION:   Generate initial conditions from distribution function
	     of pair density potential choosen.
 INPUTS:     structure 'cluster' and parameters arrays 
 PARAMETERS:							    
 RETURN:     0							    
 LOG:        01/06/2011    Creation
********************************************************************/
int initial_DF( struct cluster galaxy[1], 
		  double prmt[NMAX2] )
{
    int i,j,k,l,m,n,eng;
    int ri, zi;
    
    //Spherical coordinates
    double rr, th, ph;
    //Cilyndrical coordinates
    double zz, vr;
    
    double rho1, rho2;
    double R_grid[NMAX2];
    double Z_grid[NMAX2];
    int RES;
    int ZRES;
    double R[3];
    
    double V_grid[NMAX2][2];
    double EL_grid[NMAX2][NMAX2][4];
    int VRES;
    int ERES;
    int LRES;
    double DFnorm;
    double phi_i, phi0, phi1;
    double dN, dNv;
    double sg;
    
    double Mtot;
    
    Mtot = prmt[MASS]*( galaxy[0].N );
    
    srand48(time(NULL));
    //BULGE PARTICLES====================================================================
    //Bulge parameters
    double M_B, R_B, R_B_off;
    //Radial resolution of grid
    RES = 50;
    //Velocity resolution of grid
    VRES = 10;
    
    //Bulge Mass
    M_B = galaxy[0].pot_conf[BL_PAR+0];
    //Bulge caracteristic lenght
    R_B = galaxy[0].pot_conf[BL_PAR+1];
    //Cut off length
    R_B_off = galaxy[0].pot_conf[BL_PAR+4];
    
    switch( galaxy[0].distros[BULGE] ){
	case 0:{
	    break;}
	case 1:{
	    //Center density
	    R[X] = R_B_off/(RES*RES); R[Y] = 0; R[Z] = 0;
	    rho1 = plummer_bulge_rho( R, R, galaxy[0].pot_conf );
	    phi0 = total_phi( R, R, galaxy[0].pot_conf, galaxy[0].distros );
	    
	    //Cut off radius density
	    R[X] = R_B_off; R[Y] = 0; R[Z] = 0;
	    rho2 = plummer_bulge_rho( R, R, galaxy[0].pot_conf );
	    
	    for( i=0; i<=RES; i++ ){
		R_grid[i] = plummer_bulge_inv( rho1+i*(rho2-rho1)/RES, galaxy[0].pot_conf );}
	    
	    break;}}
    
    i = 0;
    while( i<galaxy[0].N3 ){
	for( ri=0; ri<RES; ri++ ){
	    switch( galaxy[0].distros[BULGE] ){
		case 0:{
		    dN = 0;
		    break;}
		case 1:{
		    R[X] = R_grid[ri]; R[Y] = 0; R[Z] = 0;
		    rr = R_grid[ri+1] - R_grid[ri];
		    dN = (plummer_bulge_rho( R, R, galaxy[0].pot_conf )*(4*M_PI*R_grid[ri]*R_grid[ri]*rr)/M_B)*galaxy[0].N3;
		    
		    //Velocities
		    phi_i = total_phi( R, R, galaxy[0].pot_conf, galaxy[0].distros );
		    
		    V_grid[VRES][0] = -2*phi0;
		    DFnorm = 0;
		    for( l=0; l<VRES; l++ ){
			V_grid[l][0] = V_grid[VRES][0]*(l+1)/VRES;
			V_grid[l][1] = plummer_bulge_DF( 0.5*pow(V_grid[l][0],2) + phi_i, galaxy[0].pot_conf )*\
			4*M_PI*V_grid[l][0]*V_grid[l][0]*V_grid[VRES][0]/VRES;
			DFnorm += V_grid[l][1];}
		    for( l=0; l<VRES; l++ ){
			V_grid[l][1] = V_grid[l][1]*dN/DFnorm;
			if( V_grid[l][1]>0.5 && V_grid[l][1]<1 ) V_grid[l][1] = 1;}
			
		    break;}}
				
	    //Positions
	    l = 0;
	    m = 0;
	    for( j=0; j<dN; j++ ){
		k = galaxy[0].N1+galaxy[0].N2+i+j;
		if( k==galaxy[0].N )
		    break;
		//Positions
		rr = R_grid[ri] + drand48()*(R_grid[ri+1]-R_grid[ri]);
		th = M_PI*drand48();
		ph = 2*M_PI*drand48();
		galaxy[0].prts[k].r[X] = rr*cos(ph)*sin(th);
		galaxy[0].prts[k].r[Y] = rr*sin(ph)*sin(th);
		galaxy[0].prts[k].r[Z] = rr*cos(th);
		galaxy[0].prts[k].m = prmt[MASS];
		galaxy[0].prts[k].label = k+1;
		galaxy[0].prts[k].type = 3;

 		if( m <= (int)V_grid[l][1] ){
 		    m++;}
		
		while( (int)V_grid[l][1]==0 || isnan(V_grid[l][1])==1 ){
		    l++;
		    m = 0;
		    if( l == VRES )	l==0;}
		
		do{	
		    galaxy[0].prts[k].v[X] = (1-2*drand48())*V_grid[l][0];
		    galaxy[0].prts[k].v[Y] = (1-2*drand48())*V_grid[l][0];
		    galaxy[0].prts[k].v[Z] = (1-2*drand48())*V_grid[l][0];
		}while( (pow(galaxy[0].prts[k].v[X],2)+pow(galaxy[0].prts[k].v[Y],2)+\
		pow(galaxy[0].prts[k].v[Z],2) - pow(V_grid[l][0],2))/pow(V_grid[l][0],2)<=0.0001);
		
		
		R[X] = rr; R[Y] = 0; R[Z] = 0;
		phi_i = total_phi( R, R, galaxy[0].pot_conf, galaxy[0].distros );
		
		//Total energy
		galaxy[0].prts[k].E = 0.5*( pow(galaxy[0].prts[k].v[X],2)+pow(galaxy[0].prts[k].v[Y],2)+\
		pow(galaxy[0].prts[k].v[Z],2) ) + phi_i;
		//Total Angular momentum z
		galaxy[0].prts[k].Lz = (galaxy[0].prts[k].r[X]*galaxy[0].prts[k].v[Y] - galaxy[0].prts[k].r[Y]*galaxy[0].prts[k].v[X] );

	    }
		
	    i += dN;}}
		    
		    
	    
    //DISK PARTICLES====================================================================
    double M_D, R_D, z0, R_D_off;
    //Radial resolution of grid
    RES = 100;
    //Vertical resolution of grid
    ZRES = 10;
    //Velocity resolution of grid
    VRES = 10;
    //Energy resolution of grid
    ERES = 10;
    //Energy angular of grid
    LRES = 10;

    //Disk Mass
    M_D = galaxy[0].pot_conf[HD_PAR+0];
    //Disk caracteristic lenght
    R_D = galaxy[0].pot_conf[HD_PAR+1];
    //Disk caracteristic height
    z0 = galaxy[0].pot_conf[HD_PAR+2];
    //Cut off length
    R_D_off = galaxy[0].pot_conf[HD_PAR+4];

    
    switch( galaxy[0].distros[HOM_DIS] ){
	case 0:{
	    break;}
	case 2:{
	    //Center density
	    R[X] = 0; R[Y] = 0; R[Z] = 0;
	    rho1 = miyamoto_disk_rho2( R, R, galaxy[0].pot_conf );
	    phi0 = total_phi( R, R, galaxy[0].pot_conf, galaxy[0].distros );
	    
	    //Cut off radius density
	    R[X] = R_D_off; R[Y] = 0; R[Z] = 0;
	    rho2 = miyamoto_disk_rho2( R, R, galaxy[0].pot_conf );
	    phi1 = total_phi( R, R, galaxy[0].pot_conf, galaxy[0].distros );

	    for( i=0; i<=RES; i++ ){
		R_grid[i] = miyamoto_disk_inv( rho1+i*(rho2-rho1)/RES, galaxy[0].pot_conf );}
		
	    for( i=0; i<=ZRES; i++ ){
		Z_grid[i] = -3*z0 + 6*i*z0/ZRES;}

	    break;}}

    i = 0;
    while( i<galaxy[0].N1 ){
	for( ri=0; ri<RES; ri++ ){
	for( zi=0; zi<ZRES; zi++ ){
	    switch( galaxy[0].distros[HOM_DIS] ){
		case 0:{
		    dN = 0;
		    break;}
		case 2:{
  		    R[X] = R_grid[ri]; R[Y] = 0; R[Z] = Z_grid[zi];
		    rr = R_grid[ri+1] - R_grid[ri];
		    zz = fabs(Z_grid[1] - Z_grid[0]);
 		    dN = (miyamoto_disk_rho( R, R, galaxy[0].pot_conf )*(M_PI*(R_grid[ri]+R_grid[ri+1])*rr*zz)/M_D)*galaxy[0].N1;
		    //Currente potential
// 		    phi_i = miyamoto_disk_phi0( R, R, galaxy[0].pot_conf )/3;
		    phi_i = total_phi( R, R, galaxy[0].pot_conf, galaxy[0].distros )/3;
		    
		    //Energy - Angular Momentum distribution
		    DFnorm = 0;
		    for( l=1; l<=ERES; l++ ){
		    for( j=0; j<LRES; j++ ){
			EL_grid[l][j][0] = phi_i*(l)/ERES;
			EL_grid[l][j][1] = R_grid[ri]*sqrt(2*(EL_grid[l][j][0]-phi_i))*(j+1)/LRES;
			EL_grid[l][j][2] = miyamoto_disk_DF( EL_grid[l][j][0] , EL_grid[l][l][1], galaxy[0].pot_conf  );
			DFnorm += EL_grid[l][j][2];
			EL_grid[l][j][3] = DFnorm;
		    }}
		    
		    break;}}
		    
	    //Normalization Factor
	    DFnorm = dN/EL_grid[ERES][LRES-1][3];

	    //Positions
	    m = 0;
	    l = 0;
	    n = 1;
	    for( j=0; j<dN; j++ ){
		k = i+j;
		if( k<galaxy[0].N1 ){
		//Positions
		rr = R_grid[ri] + drand48()*(R_grid[ri+1]-R_grid[ri]);
		th = 2*M_PI*drand48();
		zz = Z_grid[zi] + drand48()*(Z_grid[zi+1]-Z_grid[zi]);
		//Current Potential
		R[X] = rr; R[Y] = 0; R[Z] = zz;
// 		phi_i = miyamoto_disk_phi0( R, R, galaxy[0].pot_conf );
		phi_i = total_phi( R, R, galaxy[0].pot_conf, galaxy[0].distros )/3;
		
		galaxy[0].prts[k].r[X] = rr*cos(th);
		galaxy[0].prts[k].r[Y] = rr*sin(th);
		galaxy[0].prts[k].r[Z] = zz;
		galaxy[0].prts[k].m = prmt[MASS];
		galaxy[0].prts[k].label = k+1;
		galaxy[0].prts[k].type = 1;
		
		if( m >= EL_grid[n][l][2]*DFnorm ){
 		    l++;
		    m = 0; }
  		if( l >= LRES ){
 		    l = 0;
		    n++; }
		if( n >= ERES ){
		    n = 1; }
		
 		while( (int)(EL_grid[n][l][3]*DFnorm)==0.0 || isnan(EL_grid[n][l][2]*DFnorm)==1 ){
   		    l++;
		    if( l >= LRES ){
			l = 0;
			n++; }
  		    m = 0;}

  		do{
		EL_grid[n][l][1] = R_grid[ri]*sqrt(2*(EL_grid[n][l][0]-phi_i))*(l+1)/LRES;
		
		galaxy[0].prts[k].v[Z] = (1-2*drand48())*sqrt( 2*(EL_grid[n][l][0] - phi_i) - \
		pow(EL_grid[n][l][1]/rr,2) )*z0/R_D;
		
		sg = (1-2*drand48()); sg = sg/fabs(sg);
		vr = sg*sqrt( 2*(EL_grid[n][l][0] - phi_i) - \
		pow(EL_grid[n][l][1]/rr,2) - pow(galaxy[0].prts[k].v[Z],2) );
		
		galaxy[0].prts[k].v[X] = vr*cos(th) - sin(th)*EL_grid[n][l][1]/rr;
		galaxy[0].prts[k].v[Y] = vr*sin(th) + cos(th)*EL_grid[n][l][1]/rr;
				
 		}while( fabs(EL_grid[n][l][1]-(galaxy[0].prts[k].r[X]*galaxy[0].prts[k].v[Y] - galaxy[0].prts[k].r[Y]*galaxy[0].prts[k].v[X] ))>=0.01 || 
  		(galaxy[0].prts[k].r[X]*galaxy[0].prts[k].v[Y] - galaxy[0].prts[k].r[Y]*galaxy[0].prts[k].v[X] )<0); 
				
		//Total energy
		galaxy[0].prts[k].E = 0.5*( pow(galaxy[0].prts[k].v[X],2)+pow(galaxy[0].prts[k].v[Y],2)+\
		pow(galaxy[0].prts[k].v[Z],2) ) + 3*phi_i;
		//Total Angular momentum z
		galaxy[0].prts[k].Lz = (galaxy[0].prts[k].r[X]*galaxy[0].prts[k].v[Y] - galaxy[0].prts[k].r[Y]*galaxy[0].prts[k].v[X] );
		
		m++;
	    }}
	    i += dN;}}}
	    
	    
    //HALO PARTICLES====================================================================
    //Radial resolution of grid
    RES = 50;
    //Velocity resolution of grid
    VRES = 10;
    
    //Halo Mass
    M_B = galaxy[0].pot_conf[DH_PAR+0];
    //Halo caracteristic lenght
    R_B = galaxy[0].pot_conf[DH_PAR+1];
    //Cut off length
    R_B_off = galaxy[0].pot_conf[DH_PAR+4];
    
    switch( galaxy[0].distros[HALO_DK] ){
	case 0:{
	    break;}
	case 2:{
	    //Center density
	    R[X] = R_B_off/(RES*RES); R[Y] = 0; R[Z] = 0;
	    rho1 = plummer_halo_rho( R, R, galaxy[0].pot_conf );
	    phi0 = total_phi( R, R, galaxy[0].pot_conf, galaxy[0].distros );
	    
	    //Cut off radius density
	    R[X] = R_B_off; R[Y] = 0; R[Z] = 0;
	    rho2 = plummer_halo_rho( R, R, galaxy[0].pot_conf );
	    
	    for( i=0; i<=RES; i++ ){
		R_grid[i] = plummer_halo_inv( rho1+i*(rho2-rho1)/RES, galaxy[0].pot_conf );}
	    
	    break;}}
    
    i = 0;
    while( i<galaxy[0].N2 ){
	for( ri=0; ri<RES; ri++ ){
	    switch( galaxy[0].distros[HALO_DK] ){
		case 0:{
		    dN = 0;
		    break;}
		case 2:{
		    R[X] = R_grid[ri]; R[Y] = 0; R[Z] = 0;
		    rr = R_grid[ri+1] - R_grid[ri];
		    dN = (plummer_halo_rho( R, R, galaxy[0].pot_conf )*(4*M_PI*R_grid[ri]*R_grid[ri]*rr)/M_B)*galaxy[0].N2;
		    
		    //Velocities
		    phi_i = total_phi( R, R, galaxy[0].pot_conf, galaxy[0].distros );
		    
		    V_grid[VRES][0] = -2*phi0;
		    DFnorm = 0;
		    for( l=0; l<VRES; l++ ){
			V_grid[l][0] = V_grid[VRES][0]*(l+1)/VRES;
			V_grid[l][1] = plummer_halo_DF( 0.5*pow(V_grid[l][0],2) + phi_i, galaxy[0].pot_conf )*\
			4*M_PI*V_grid[l][0]*V_grid[l][0]*V_grid[VRES][0]/VRES;
			DFnorm += V_grid[l][1];}
		    for( l=0; l<VRES; l++ ){
			V_grid[l][1] = V_grid[l][1]*dN/DFnorm;
			if( V_grid[l][1]>0.5 && V_grid[l][1]<1 ) V_grid[l][1] = 1;}
			
		    break;}}
				
	    //Positions
	    l = 0;
	    m = 0;
	    for( j=0; j<dN; j++ ){
		k = galaxy[0].N1+i+j;
		if( k<=galaxy[0].N1 + galaxy[0].N2 ){
		//Positions
		rr = R_grid[ri] + drand48()*(R_grid[ri+1]-R_grid[ri]);
		th = M_PI*drand48();
		ph = 2*M_PI*drand48();
		galaxy[0].prts[k].r[X] = rr*cos(ph)*sin(th);
		galaxy[0].prts[k].r[Y] = rr*sin(ph)*sin(th);
		galaxy[0].prts[k].r[Z] = rr*cos(th);
		galaxy[0].prts[k].m = prmt[MASS];
		galaxy[0].prts[k].label = k+1;
		galaxy[0].prts[k].type = 2;

 		if( m <= (int)V_grid[l][1] ){
 		    m++;}
		
		while( (int)V_grid[l][1]==0 || isnan(V_grid[l][1])==1 ){
		    l++;
		    m = 0;
		    if( l == VRES )	l==0;}
		
		do{	
		    galaxy[0].prts[k].v[X] = (1-2*drand48())*V_grid[l][0];
		    galaxy[0].prts[k].v[Y] = (1-2*drand48())*V_grid[l][0];
		    galaxy[0].prts[k].v[Z] = (1-2*drand48())*V_grid[l][0];
		}while( (pow(galaxy[0].prts[k].v[X],2)+pow(galaxy[0].prts[k].v[Y],2)+\
		pow(galaxy[0].prts[k].v[Z],2) - pow(V_grid[l][0],2))/pow(V_grid[l][0],2)<=0.0001);
		
		
		R[X] = rr; R[Y] = 0; R[Z] = 0;
		phi_i = total_phi( R, R, galaxy[0].pot_conf, galaxy[0].distros );
		
		//Total energy
		galaxy[0].prts[k].E = 0.5*( pow(galaxy[0].prts[k].v[X],2)+pow(galaxy[0].prts[k].v[Y],2)+\
		pow(galaxy[0].prts[k].v[Z],2) ) + phi_i;
		//Total Angular momentum z
		galaxy[0].prts[k].Lz = (galaxy[0].prts[k].r[X]*galaxy[0].prts[k].v[Y] - galaxy[0].prts[k].r[Y]*galaxy[0].prts[k].v[X] );

	    }}
		
	    i += dN;}}
	    
		    
    return 0;
}