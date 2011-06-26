// Data Module
int conf2dump( char * );
int in2dump( char * );
int read_parameters( double *, char * );
int data_in( struct cluster *, char * );
int data_out( struct cluster *, char * );

// Numeric Module
double distance( double *, double * );
int factorial( int  );
int integration_step( struct cluster *, double , double , int );
int euler_step( struct cluster *, double , double );
int rk4_step( struct cluster *, double , double );

// Potentilas Module
int dynamic_function( struct cluster *, double , struct intg_step * );
double total_phi( double *, double *, double *, int * );

int miyamoto_disk_phi( double *, double *, double , double *, double * );
double miyamoto_disk_phi0( double *, double *, double * );
double miyamoto_disk_rho( double *, double *, double * );
double miyamoto_disk_rho2( double *, double *, double * );
double miyamoto_disk_inv( double , double * );
double miyamoto_disk_DF( double , double,  double * );

int plummer_bulge_phi( double *, double *, double , double *, double * );
double plummer_bulge_phi0( double *, double *, double * );
double plummer_bulge_rho( double *, double *, double * );
double plummer_bulge_inv( double , double * );
double plummer_bulge_DF( double , double * );

int plummer_halo_phi( double *, double *, double , double *, double * );
double plummer_halo_phi0( double *, double *, double * );
double plummer_halo_rho( double *, double *, double * );
double plummer_halo_inv( double , double * );
double plummer_halo_DF( double , double * );

int miyamoto_bar_phi( double *, double *, double , double *, double * );

int NFW_phi( double *, double *, double , double *, double * );
double NFW_phi0( double *, double *, double * );

// Initial Distribution Module
int initial_conditions( struct cluster *, double *, char * );
int uniform_cube( struct cluster *, double * );