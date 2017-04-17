#ifndef _KALMANH_
#define _KALMANH_

#include "MKL25Z4.h"                    // Device header
#include "arm_math.h"                   // ARM::CMSIS:DSP
#include "math.h"


#define _DIM_	8  											// System dimension 
#define _NOI_	6												// Number of inputs

#define EXIT_SUCCES 0

typedef float32_t data_t;
typedef arm_matrix_instance_f32  array_t ;


//Data type specific macros
#define array_init(a,b,c,d) 	(arm_mat_init_f32(a,b,c,d))
#define sqrtf(a,b) 						(arm_sqrt_f32(a, b))
#define array_transpose(a,b)	(arm_mat_trans_f32(a,b))
#define array_mult(a,b,c)     (arm_mat_mult_f32(a,b,c))
#define array_add(a,b,c)      (arm_mat_add_f32(a,b,c))
#define array_inv(a,b)				(arm_mat_inverse_f32(a,b))
#define array_sub(a,b,c)      (arm_mat_sub_f32(a,b,c))
	

	
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

#define X_Q0 0
#define X_Q1 1
#define X_Q2 2
#define X_Q3 3
#define X_W1 4
#define X_W2 5
#define X_W3 6
#define X_B 7

//===============================================================/
/**
*@brief system  state temp vector 1(matrix)
**/


#define X_T1_ROWS _DIM_
#define X_T1_COLS 1

array_t X_t1_k;
/**
*@brief system state temp vector 1 data array
**/
data_t X_t1_k_data[X_T1_ROWS*X_T1_COLS] = {0};

#define X_T1_Q0 0
#define X_T1_Q1 1
#define X_T1_Q2 2
#define X_T1_Q3 3
#define X_T1_W1 4
#define X_T1_W2 5
#define X_T1_W3 6
#define X_T1_B  7

//===============================================================/
/**
*@brief system  state temp vector 2 (matrix)
**/


#define X_T2_ROWS _DIM_
#define X_T2_COLS 1

array_t X_t2_k;
/**
*@brief system state temp vector 2 data array
**/
data_t X_t2_k_data[X_T2_ROWS*X_T2_COLS] = {0};

#define X_T2_Q0 0
#define X_T2_Q1 1
#define X_T2_Q2 2
#define X_T2_Q3 3
#define X_T2_W1 4
#define X_T2_W2 5
#define X_T2_W3 6
#define X_T2_B  7
//===============================================================/
/**
*@brief system state vector estimate apriori 
**/
#define X_AP_ROWS X_ROWS
#define X_AP_COLS X_COLS

array_t X_ap_k;
/**
*@brief system state vector estimate apriori data array
**/
data_t X_ap_k_data[X_AP_ROWS*X_AP_COLS] = {0};


//===============================================================/
/**
*@brief measurements vector (matrix)
**/

#define Z_ROWS _NOI_
#define Z_COLS 1

array_t Z_k;
/**
*@brief measurements vector data array
**/
data_t Z_k_data[Z_ROWS*Z_COLS] = {0};

#define Z_B1 0
#define Z_B2 1
#define Z_B3 2
#define Z_W1 3
#define Z_W2 4
#define Z_W3 5

//===============================================================/

/**
*@brief innovation vector (matrix)
**/
#define Y_ROWS Z_ROWS
#define Y_COLS Z_COLS

array_t Y_k;
/**
*@brief innovation vector data array
**/
data_t Y_k_data[Y_ROWS*Y_COLS] = {0};

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
{0.00762129, 0.00436500, 0.00436500, 0.00436500, 0.00436500, 8.73e-4, 8.73e-4, 8.73e-7,
 0.00436500, 0.00250000, 0.00250000, 0.00250000, 0.00250000, 5e-4, 		5e-4, 	 5e-7,
 0.00436500, 0.00250000, 0.00250000, 0.00250000, 0.00250000, 5e-4, 		5e-4,		 5e-7,
 0.00436500, 0.00250000, 0.00250000, 0.00250000, 0.00250000, 5e-4, 		5e-4,    5e-7,
 8.73e-4,		 5e-4, 			 5e-4,			 5e-4, 			 5e-4, 			 1e-4, 	  1e-4,    1e-7,
 8.73e-4,		 5e-4, 			 5e-4, 			 5e-4,       5e-4,			 1e-4, 		1e-4,    1e-7,
 8.73e-4, 	 5e-4, 			 5e-4, 			 5e-4,       5e-4, 			 1e-4, 		1e-4,		 1e-7,
 8.73e-7, 	 5e-7, 			 5e-7, 			 5e-7,       5e-7, 			 1e-7, 		1e-7,		 1e-10};

 
 //===============================================================/
 /**
*@breif Measurement noise covariance matrix of phase 1
**/
#define RPH1_ROWS _NOI_
#define RPH1_COLS _NOI_
 
array_t Rph1_k;
 /**
*@breif SMeasurement noise covariance matrix phase 1 data array
**/ 
data_t	Rph1_k_data[RPH1_COLS*RPH1_ROWS] =
{3.844e-13,				  0, 				 0, 			 0, 			 0, 
 0, 		 3.844e-13, 0,				 0, 			 0, 			 0, 
 0, 		 0,         3.844e-13, 0, 			 0,     	 0, 
 0, 		 0, 				0, 				 1.369e-5, 0, 			 0, 
 0, 		 0, 				0, 				 0,        1.369e-5, 0, 
 0, 		 0, 				0, 				 0, 			 0, 			 1.369e-5};

 
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
#define RPH2_ROWS _NOI_
#define RPH2_COLS _NOI_
 
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

//---------------------------------------------------------------/
 
array_t A_k_t;
/**
 *@brief Global jacoban of f function transposition 
*/
data_t A_k_t_data[A_COLS*A_ROWS] = {0}; 

 
//===============================================================/
/**
*@breif Global jacobian of h funciton (H^J) matrix
**/
#define H_ROWS _NOI_
#define H_COLS _DIM_
 
array_t H_k;
/**
 *@brief Global jacoban of h function (H^J) matrix data array
*/

data_t	H_k_data[H_ROWS*H_COLS] =
{0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 1, 0, 0,
 0, 0, 0, 0, 0, 0, 1, 0,
 0, 0, 0, 0, 0, 0, 0, 1};

//---------------------------------------------------------------/
 
array_t H_k_t;
/**
 *@brief Global jacoban of h function transposition 
*/
data_t H_k_t_data[H_COLS*H_ROWS] = {0}; 

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

//===============================================================/
/**
*@breif Global state prediction covariance matrix estimate apriori
**/
#define P_AP_ROWS P_ROWS
#define P_AP_COLS P_COLS
 
array_t P_ap_k;
 
/**
*@breif Global state pretiction covariance matrix estimate apriori data array
**/
 
data_t	P_ap_k_data[P_AP_ROWS*P_AP_COLS] =
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
*@breif Global kalman gain matrix
**/
#define KG_ROWS _DIM_
#define KG_COLS _NOI_
 
array_t Kg_k;
 
/**
*@breif Global kalman gain matrix data array
**/
 
data_t	Kg_k_data[KG_ROWS*KG_COLS] = {0};

//===============================================================/
/**
*@brief Scaling matrix for better conditioning
**/
#define S_ROWS _NOI_
#define S_COLS _NOI_

array_t S_k;

/**
*@breif Scaling matrix data array
**/
 
data_t S_k_data[S_ROWS*S_COLS] =
{10e6, 0, 	 0, 	 0,				 0, 			   0, 
 0, 	 10e6, 0, 	 0,				 0,     		 0, 
 0, 	 0, 	 10e6, 0, 			 0,          0, 
 0, 	 0, 	 0, 	 31.6228,	 0,   			 0,         //10*sqrt(10)
 0, 	 0, 	 0, 	 0, 			 31.6228,    0,
 0,		 0,		 0, 	 0,				 0, 			   31.6228};	

//===============================================================/
/**
*@brief V matrix 
**/
#define V_ROWS _NOI_
#define V_COLS _NOI_

array_t V_k;

/**
*@breif Vk data array
**/
 
data_t V_k_data[S_ROWS*S_COLS] = 
{1, 0, 0, 0, 0, 0,
 0, 1, 0, 0, 0, 0,
 0, 0, 1, 0, 0, 0,
 0, 0, 0, 1, 0, 0,
 0, 0, 0, 0, 1, 0, 
 0, 0, 0, 0, 0, 1};

//---------------------------------------------------------------/
 
array_t V_k_t;
/**
 *@brief Vk transposition 
*/
data_t V_k_t_data[V_COLS*V_ROWS] = {0}; 

//===============================================================/
/**
*@breif W matrix
**/
#define W_ROWS _DIM_
#define W_COLS _DIM_

array_t W_k;
/**
*@breif W matrix data array
**/
data_t	W_k_data[W_ROWS*W_COLS] =
{1, 0, 0, 0, 0, 0, 0, 0,
 0, 1, 0, 0, 0, 0, 0, 0,
 0, 0, 1, 0, 0, 0, 0, 0,
 0, 0, 0, 1, 0, 0, 0, 0,
 0, 0, 0, 0, 1, 0, 0, 0,
 0, 0, 0, 0, 0, 1, 0, 0, 
 0, 0, 0, 0, 0, 0, 1, 0,
 0, 0, 0, 0, 0, 0, 0, 1};

//---------------------------------------------------------------/
 
array_t W_k_t;
/**
 *@brief Global jacoban of h function transposition 
*/
data_t W_k_t_data[W_COLS*W_ROWS] = {0}; 

 
 //===============================================================/
/**
*@breif Temporary matrix for storing results
**/
#define TEMP_1_ROWS _DIM_
#define TEMP_1_COLS _DIM_

array_t Temp_1;
/**
*@breif Temporary matrix data array
**/
data_t	Temp_1_data[TEMP_1_ROWS*TEMP_1_COLS] =
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
*@breif Temporary matrix for storing results
**/
#define TEMP_2_ROWS _DIM_
#define TEMP_2_COLS _DIM_

array_t Temp_2;
/**
*@breif Temporary matrix data array
**/
data_t	Temp_2_data[TEMP_2_ROWS*TEMP_2_COLS] =
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
*@breif Temporary matrix for storing results
**/
#define TEMP_3_ROWS _DIM_
#define TEMP_3_COLS _DIM_

array_t Temp_3;
/**
*@breif Temporary matrix data array
**/
data_t	Temp_3_data[TEMP_3_ROWS*TEMP_3_COLS] =
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
*@breif Temporary matrix for storing results - Num of inputs x Num of states
**/
#define TEMP_4_ROWS _NOI_
#define TEMP_4_COLS _DIM_

array_t Temp_4;
/**
*@breif Temporary matrix data array
**/
data_t	Temp_4_data[TEMP_4_ROWS*TEMP_4_COLS] =
{1, 0, 0, 0, 0, 0,
 0, 1, 0, 0, 0, 0,
 0, 0, 1, 0, 0, 0,
 0, 0, 0, 1, 0, 0,
 0, 0, 0, 0, 1, 0,
 0, 0, 0, 0, 0, 1};

  //===============================================================/
 /**
*@breif Temporary matrix for storing results - Num of inputs x Num of inputs
**/
#define TEMP_5_ROWS _NOI_
#define TEMP_5_COLS _NOI_

array_t Temp_5;
/**
*@breif Temporary matrix data array
**/
data_t	Temp_5_data[TEMP_5_ROWS*TEMP_5_COLS] =
{1, 0, 0, 0, 0, 0,
 0, 1, 0, 0, 0, 0,
 0, 0, 1, 0, 0, 0,
 0, 0, 0, 1, 0, 0,
 0, 0, 0, 0, 1, 0,
 0, 0, 0, 0, 0, 1};

 //===============================================================/
 /**
*@breif Temporary matrix for storing results - Num of inputs x Num of inputs
**/
#define TEMP_6_ROWS _NOI_
#define TEMP_6_COLS _NOI_

array_t Temp_6;
/**
*@breif Temporary matrix data array
**/
data_t	Temp_6_data[TEMP_6_ROWS*TEMP_6_COLS] =
{1, 0, 0, 0, 0, 0,
 0, 1, 0, 0, 0, 0,
 0, 0, 1, 0, 0, 0,
 0, 0, 0, 1, 0, 0,
 0, 0, 0, 0, 1, 0,
 0, 0, 0, 0, 0, 1};

 //===============================================================/
 /**
*@breif Temporary matrix for storing results - Num of inputs x Num of inputs
**/
#define TEMP_7_ROWS _NOI_
#define TEMP_7_COLS _NOI_

array_t Temp_7;
/**
*@breif Temporary matrix data array
**/
data_t	Temp_7_data[TEMP_7_ROWS*TEMP_7_COLS] =
{1, 0, 0, 0, 0, 0,
 0, 1, 0, 0, 0, 0,
 0, 0, 1, 0, 0, 0,
 0, 0, 0, 1, 0, 0,
 0, 0, 0, 0, 1, 0,
 0, 0, 0, 0, 0, 1};

/******************************************************************
***************KALMAN FILTER SPECIFIC FUNCTIONS********************
*******************************************************************
*Extended Kalman filter 
* x = f(x_k-1, u_k-1) + w_k-1;
* z = h(x_k) + v_k
*/
 
 
/**
*@biref computes f function
*/
void kalmanF		(	data_t q0,
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
*@biref computes output innovation
*  y = z_k - h(x_ap)
*/
void kalmanInv  (	data_t q0,
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
void jacobianHph1(data_t q0,
									data_t q1,
									data_t q2,
									data_t q3,
									data_t w1,
									data_t w2,
									data_t w3, 
									data_t B);



/******************************************************************
***************KALMAN FILTER API FUNCTIONS*************************
*******************************************************************/
/**
*@breif initialization of Kalman filter in terms of data structures
*Initializes arm_math arrays and ...
*/
void kallmanInit();

/**
*@breif reads data from sensors
*/
void kalmanGetSensors();

/**
*@breif state vector initalizatnion 
*calls kalmanGetSenors()
*/
void kalmanInitializeStateVec();

/**
*@breif Kalman filter iteration
*/
void kalmanStep();
#endif

