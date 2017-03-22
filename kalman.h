#ifndef _KALMANH_
#define _KALMANH_

#include "MKL25Z4.h"                    // Device header
#include "arm_math.h"                   // ARM::CMSIS:DSP


#define _DIM_	8  											// State vector length 
#define _NME_	6												// Number of measurements 

#define EXIT_SUCCES 0

typedef float32_t data_t;
typedef arm_matrix_instance_f32  array_t ;

#define array_init(a,b,c,d) arm_mat_init_f32(a,b,c,d)
	

	
/********************************************************
****************MODEL PARAMTERES*************************
*********************************************************/

/**
*@breif Global sample time
*/
const data_t sampleTime = 0.1f;	

/**
*@breif Moment of inertia  on Xaxis 
*/
const data_t I1 = 0.1f;
/**
*@breif Moment of inertia  on Yaxis 
*/
const data_t I2 = 0.1f;
/**
*@breif Moment of inertia  on Zaxis 
*/
const data_t I3 = 0.1f;						



/********************************************************
*****DEFINITON AND INITIALIZATION OF DATA STUCTURES******
*********************************************************/
/**
* Every cell of array that is not modified in  dedicated funciton 
*	needs be initialised here with its proper value
*/


//===============================================================/
/**
*@brief system state vector (matrix)
**/

#define X_ROWS _DIM_
#define X_COLS 1

array_t X_k;
/**
*@brief system state vector data array
**/
data_t X_k_data[X_ROWS*X_COLS] = {0};

	
//===============================================================/
/**
*@brief measurements vector (matrix)
**/

#define Z_ROWS _NME_
#define Z_COLS 1

array_t Z_k;
/**
*@brief measurements vector data array
**/
data_t Z_k_data[Z_ROWS*Z_COLS] = {0};



//===============================================================/
/**
*@breif System noise covariance matrix of phase 1
**/
#define QPH1_ROWS _DIM_
#define QPH1_COLS _DIM_

array_t Qph1_k;
/**
*@breif System noise covariance matrix phase 1 data array
**/
data_t	Qph1_k_data[QPH1_ROWS*QPH1_COLS] =
{1, 0, 0, 0, 0, 0, 0, 0,
 0, 1, 0, 0, 0, 0, 0, 0,
 0, 0, 1, 0, 0, 0, 0, 0,
 0, 0, 0, 1, 0, 0, 0, 0,
 0, 0, 0, 0, 1, 0, 0, 0,
 0, 0, 0, 0, 0, 1, 0, 0, 
 0, 0, 0, 0, 0, 0, 1, 0,
 0, 0, 0, 0, 0, 0, 0, 1};

 
 //===============================================================/
 /**
*@breif Measurement noise covariance matrix of phase 1
**/
#define RPH1_ROWS _NME_
#define RPH1_COLS _NME_
 
array_t Rph1_k;
 /**
*@breif SMeasurement noise covariance matrix phase 1 data array
**/ 
data_t	Rph1_k_data[RPH1_COLS*RPH1_ROWS] =
{1, 0, 0, 0, 0, 0, 
 0, 1, 0, 0, 0, 0, 
 0, 0, 1, 0, 0, 0, 
 0, 0, 0, 1, 0, 0, 
 0, 0, 0, 0, 1, 0, 
 0, 0, 0, 0, 0, 1};

 
//===============================================================/
/**
*@breif System noise covariance matrix of phase 2
**/
#define QPH2_ROWS _DIM_
#define QPH2_COLS _DIM_ 

array_t Qph2_k;
 /**
*@breif System noise covariance matrix phase 2 data array
**/
data_t	Qph2_k_data[QPH2_ROWS*QPH2_COLS] =
{1, 0, 0, 0, 0, 0, 0, 0,
 0, 1, 0, 0, 0, 0, 0, 0,
 0, 0, 1, 0, 0, 0, 0, 0,
 0, 0, 0, 1, 0, 0, 0, 0,
 0, 0, 0, 0, 1, 0, 0, 0,
 0, 0, 0, 0, 0, 1, 0, 0, 
 0, 0, 0, 0, 0, 0, 1, 0,
 0, 0, 0, 0, 0, 0, 0, 1};

 //===============================================================/
 /**
*@breif Measurement noise covariance matrix of phase 2
**/
#define RPH2_ROWS _NME_
#define RPH2_COLS _NME_
 
array_t Rph2_k;
 /**
 *@breif Measurement noise covariance matrix phase 2 data array
*/
data_t	Rph2_k_data[QPH2_ROWS*QPH2_COLS] =
{1, 0, 0, 0, 0, 0, 
 0, 1, 0, 0, 0, 0, 
 0, 0, 1, 0, 0, 0, 
 0, 0, 0, 1, 0, 0, 
 0, 0, 0, 0, 1, 0, 
 0, 0, 0, 0, 0, 1};

//===============================================================/
/**
*@breif Global jacobian of f funciton (A^J) matirx
**/
#define A_ROWS _DIM_
#define A_COLS _DIM_
 
array_t A_k;
/**
*@breif Global jacobian of f funciton (A^J) matrix data array
**/
 
data_t	A_k_data[A_ROWS*A_COLS] =
{1, 0, 0, 0, 0, 0, 0, 0,
 0, 1, 0, 0, 0, 0, 0, 0,
 0, 0, 1, 0, 0, 0, 0, 0,
 0, 0, 0, 1, 0, 0, 0, 0,
 0, 0, 0, 0, 1, 0, 0, 0,
 0, 0, 0, 0, 0, 1, 0, 0, 
 0, 0, 0, 0, 0, 0, 1, 0,
 0, 0, 0, 0, 0, 0, 0, 1};

//===============================================================/
/**
*@breif Global jacobian of h funciton (H^J) matrix
**/
#define H_ROWS _NME_
#define H_COLS _DIM_
 
array_t H_k;
/**
 *@brief Global jacoban of h function (H^J) matrix data array
*/

data_t	H_k_data[_NME_*_DIM_] =
{0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 1, 0, 0,
 0, 0, 0, 0, 0, 0, 1, 0,
 0, 0, 0, 0, 0, 0, 0, 1};

//===============================================================/
/**
*@breif Global state prediction covariance matrix
**/
#define P_ROWS _DIM_
#define P_COLS _DIM_
 
array_t P_k;
 
/**
*@breif Global state pretiction covariance matrix data array
**/
 
data_t	P_k_data[P_ROWS*P_COLS] =
{1, 0, 0, 0, 0, 0, 0, 0,
 0, 1, 0, 0, 0, 0, 0, 0,
 0, 0, 1, 0, 0, 0, 0, 0,
 0, 0, 0, 1, 0, 0, 0, 0,
 0, 0, 0, 0, 1, 0, 0, 0,
 0, 0, 0, 0, 0, 1, 0, 0, 
 0, 0, 0, 0, 0, 0, 1, 0,
 0, 0, 0, 0, 0, 0, 0, 1};

 
/******************************************************************
***************KALMAN FILTER SPECIFIC FUNCTIONS********************
*******************************************************************
 /*
*Extended Kalman filter 
* x = f(x_k-1, u_k-1) + w_k-1;
* z = h(x_k) + v_k
*/
 
 
/**
*@biref computes f function
*/
void kalmanF			(array_t* dest,
									data_t q0,
									data_t q1,
									data_t q2,
									data_t q3,
									data_t w1,
									data_t w2,
									data_t w3, 
									data_t B,
									data_t T1,
									data_t T2,
									data_t T3);
/**
*@biref computes h function
*/
void kalmanH 			(array_t* dest,
									data_t q0,
									data_t q1,
									data_t q2,
									data_t q3,
									data_t w1,
									data_t w2,
									data_t w3, 
									data_t B);
/**
*@biref computes Jacobian of h function
*/
void jacobianHph1 (array_t* dest,
									data_t q0,
									data_t q1,
									data_t q2,
									data_t q3,
									data_t w1,
									data_t w2,
									data_t w3, 
									data_t B);


/******************************************************************
***************KALMAN FILTER API FUNCTIONS*************************
*******************************************************************
/**
*@breif initialization of Kalman filter in terms of data structures
*Initializes arm_math arrays and ...
*/
void kallmanInit();

/**
*@breif reads data from sensors
*/
void getSensors();




#endif

