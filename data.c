#include <allvars.h>

/*****************************************************************
 NOMBRE:     conf2dump
 FUNCTION:   convert a data file text in plain text 
 INPUTS:     configuration file name
 RETORNA:    0
 LOG:        18/01/2011    Creation
*****************************************************************/
int conf2dump( char filename[] )
{
    char cmd[100];
    sprintf( cmd, "grep -v \"#\" %s | grep -v \"^$\" | gawk -F\"=\" '{print $2}' > %s.dump", 
	     filename, filename );
    system( cmd );

    return 0;
}

/*****************************************************************
 NAME:       in2dump
 FUNCTION:   convert a data file text in plain text 
 INPUTS:     configuration file name
 RETORNA:    0
 LOG:        18/01/2011    Creation
*****************************************************************/
int in2dump( char filename[] )
{
    char cmd[100];
    sprintf( cmd, "grep -v \"#\" %s > %s.dump", 
	     filename, filename );
    system( cmd );

    return 0;
}

/*****************************************************************
 NAME:       read_parameters
 FUNCTION:   read the file with given name  and type 
	     information of array given
 INPUTS:     array where it returns reading datas and file name 
	     with congigiration file
 RETURN:     0 if file read ok
	     1 if file dont exist
 LOG:        18/01/2011    Creation
*****************************************************************/
int read_parameters( double parameters[],
		     char filename[] )
{
    char cmd[100], filenamedump[100];
    int i=0;
    FILE *file;

    //File Detection
    file = fopen( filename, "r" );
    if( file==NULL ){
	printf( "  * El archivo '%s' no existe!\n\tFIN DE LA SIMULACION.\n", filename );
	return 1;}
    fclose(file);
    
    //Convertion to plain text
    conf2dump( filename );
    sprintf( filenamedump, "%s.dump", filename );
    file = fopen( filenamedump, "r" );
    
    //Reading
    while( getc( file ) != EOF ){
	fscanf( file, "%lf", &parameters[i] );
	i++;}
    
    fclose( file );
    
    printf( "  * El archivo '%s' ha sido leido!\n", filename );

    sprintf( cmd, "rm -rf %s.dump", filename );
    system( cmd );
    
    return 0;
}

/********************************************************************
 NAME:       data_in
 FUNCTION:   read input file with masses, positions, velocities and
	     softening parameters of all particles
 INPUTS:     'cluster' structure and input file name
 RETURN:     0 if file read ok
	     1 if file dont exist
 LOG:        25/05/2008    Creation
             30/06/2008    Structres implementation
             21/01/2011    Modification
********************************************************************/
int data_in( struct cluster galaxy[1],
	     char filename[] )
{
    int i=0;
    char cmd[100], filenamedump[100];
    FILE *file;
    

    //File Detection
    file = fopen( filename, "r" );
    if( file==NULL ){
	printf( "  * El archivo '%s' no existe!\n\tFIN DE LA SIMULACION.\n", filename );
	return 1;}
    fclose(file);

    //Conversed to plain text
    in2dump( filename );
    sprintf( filenamedump, "%s.dump", filename );
    file = fopen( filenamedump, "r" );
    
    //Read datas
    while( getc( file ) != EOF ){
	fscanf( file,"%d %lf %lf %lf %lf %lf %lf %lf %d", 
		&galaxy[0].prts[i].label, 
		&galaxy[0].prts[i].r[0], &galaxy[0].prts[i].r[1], &galaxy[0].prts[i].r[2],
		&galaxy[0].prts[i].v[0], &galaxy[0].prts[i].v[1], &galaxy[0].prts[i].v[2], 
		&galaxy[0].prts[i].m, &galaxy[0].prts[i].type );
	i++;}
	  
    fclose( file );

    printf( "  * El archivo '%s' ha sido leido!\n", filename );
    
    sprintf( cmd, "rm -rf %s.dump", filename );
    system( cmd );

    return 0;
}

/********************************************************************
 NAME:       data_out
 FUNCTION:   write a file with masses, positions and velocities of
	     all particles
 INPUTS:    'cluster' structure with datas, string with out file name 
 RETURN:     0
 LOG:        25/05/2008    Creation
********************************************************************/
int data_out( struct cluster galaxy[1],
	     char filename[] )
{
    int i=0;
    char filename1[NMAX2];
    FILE *file;
    
    //GENERAL FILE
    file = fopen( filename, "w" );
//     fprintf( file, "#Galaxy Simulation\t\tSnapshot Time = %lf\n", galaxy[0].t_snap ); 
//     fprintf( file, "#Masa\t\tX\t\tY\t\tZ\t\tVx\t\tVy\t\tVz\n" ); 
    //Writing
    for( i=0; i<galaxy[0].N; i++ )
	fprintf( file, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n", 
		 galaxy[0].prts[i].label, 
		 galaxy[0].prts[i].r[0], galaxy[0].prts[i].r[1], galaxy[0].prts[i].r[2],
		 galaxy[0].prts[i].v[0], galaxy[0].prts[i].v[1], galaxy[0].prts[i].v[2],
		 galaxy[0].prts[i].m, galaxy[0].prts[i].E, galaxy[0].prts[i].Lz,
		 galaxy[0].prts[i].type );
    fclose( file );
    
    //DISK FILE
    sprintf( filename1, "%s.0" , filename );
    file = fopen( filename1, "w" );
    
//     fprintf( file, "#Galaxy Simulation\t\tSnapshot Time = %lf\n", galaxy[0].t_snap ); 
//     fprintf( file, "#Masa\t\tX\t\tY\t\tZ\t\tVx\t\tVy\t\tVz\n" ); 
    //Writing
    for( i=0; i<galaxy[0].N1; i++ )
	fprintf( file, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n", 
		 galaxy[0].prts[i].label, 
		 galaxy[0].prts[i].r[0], galaxy[0].prts[i].r[1], galaxy[0].prts[i].r[2],
		 galaxy[0].prts[i].v[0], galaxy[0].prts[i].v[1], galaxy[0].prts[i].v[2],
		 galaxy[0].prts[i].m, galaxy[0].prts[i].E, galaxy[0].prts[i].Lz,
		 galaxy[0].prts[i].type );
    fclose( file );
    
    //HALO FILE	
    sprintf( filename1, "%s.1" , filename );
    file = fopen( filename1, "w" );
    
//     fprintf( file, "#Galaxy Simulation\t\tSnapshot Time = %lf\n", galaxy[0].t_snap ); 
//     fprintf( file, "#label\t\tX\t\tY\t\tZ\t\tVx\t\tVy\t\tVz\t\tMass\n" ); 
    //Writing
    for( i=galaxy[0].N1; i<galaxy[0].N1+galaxy[0].N2; i++ )
	fprintf( file, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n", 
		 galaxy[0].prts[i].label, 
		 galaxy[0].prts[i].r[0], galaxy[0].prts[i].r[1], galaxy[0].prts[i].r[2],
		 galaxy[0].prts[i].v[0], galaxy[0].prts[i].v[1], galaxy[0].prts[i].v[2],
		 galaxy[0].prts[i].m, galaxy[0].prts[i].E, galaxy[0].prts[i].Lz,
		 galaxy[0].prts[i].type );
    fclose( file );
    
    //BULGE FILE	
    sprintf( filename1, "%s.2" , filename );
    file = fopen( filename1, "w" );
    
//     fprintf( file, "#Galaxy Simulation\t\tSnapshot Time = %lf\n", galaxy[0].t_snap ); 
//     fprintf( file, "#label\t\tX\t\tY\t\tZ\t\tVx\t\tVy\t\tVz\t\tMass\n" ); 
    //Writing
    for( i=galaxy[0].N1+galaxy[0].N2; i<galaxy[0].N; i++ )
	fprintf( file, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n", 
		 galaxy[0].prts[i].label, 
		 galaxy[0].prts[i].r[0], galaxy[0].prts[i].r[1], galaxy[0].prts[i].r[2],
		 galaxy[0].prts[i].v[0], galaxy[0].prts[i].v[1], galaxy[0].prts[i].v[2],
		 galaxy[0].prts[i].m, galaxy[0].prts[i].E, galaxy[0].prts[i].Lz,
		 galaxy[0].prts[i].type );
    fclose( file );
    
    return 0;
}