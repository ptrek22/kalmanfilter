#define ARM_MATH_CM0PLUS
#include "kalman.h"



/********************************************************
****************MODEL PARAMTERES*************************
*********************************************************/

/**
*@breif Global sample time
*/
 data_t sampleTime = 0.1f;	

/**
*@breif Moment of inertia  on Xaxis 
*/
 data_t I1 = 0.1f;
/**
*@breif Moment of inertia  on Yaxis 
*/
 data_t I2 = 0.1f;
/**
*@breif Moment of inertia  on Zaxis 
*/
 data_t I3 = 0.1f;						



/********************************************************
*************INITIALIZATION OF DATA STUCTURES************
*********************************************************/

//===============================================================/
/**
*@brief system state vector (matrix)
**/
array_t X_k;
/**
*@brief system state vector data array
**/
data_t X_k_data[X_ROWS*X_COLS] = {0};

//===============================================================/
/**
*@brief system  state temp vector 1(matrix)
**/
array_t X_t1_k;
/**
*@brief system state temp vector 1 data array
**/
data_t X_t1_k_data[X_T1_ROWS*X_T1_COLS] = {0};

//===============================================================/
/**
*@brief system  state temp vector 2 (matrix)
**/
array_t X_t2_k;
/**
*@brief system state temp vector 2 data array
**/
data_t X_t2_k_data[X_T2_ROWS*X_T2_COLS] = {0};

//===============================================================/
/**
*@brief system state vector estimate apriori 
**/
array_t X_ap_k;
/**
*@brief system state vector estimate apriori data array
**/
data_t X_ap_k_data[X_AP_ROWS*X_AP_COLS] = {0};


//===============================================================/
/**
*@brief measurements vector (matrix)
**/
array_t Z_k;
/**
*@brief measurements vector data array
**/
data_t Z_k_data[Z_ROWS*Z_COLS] = {0};

//===============================================================/

/**
*@brief innovation vector (matrix)
**/
array_t Y_k;
/**
*@brief innovation vector data array
**/
data_t Y_k_data[Y_ROWS*Y_COLS] = {0};

//===============================================================/
/**
*@breif System noise covariance matrix of phase 1
**/
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
array_t Kg_k;
/**
*@breif Global kalman gain matrix data array
**/
 
data_t	Kg_k_data[KG_ROWS*KG_COLS] = {0};

//===============================================================/
/**
*@brief Scaling matrix for better conditioning
**/
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

void jacobianHph1(data_t q0,
									data_t q1,
									data_t q2,
									data_t q3,
									data_t w1,
									data_t w2,
									data_t w3, 
									data_t B)

{	//Row 1 d B(q1q3-q0q2)
	H_k.pData[Z_B1*H_COLS + X_Q0] = -2*q2*B;
	H_k.pData[Z_B1*H_COLS + X_Q1] =  2*q3*B;
	H_k.pData[Z_B1*H_COLS + X_Q2] = -2*q0*B;
	H_k.pData[Z_B1*H_COLS + X_Q3] =  2*q1*B;
//	H_k.pData[Z_B1*H_COLS + X_W1] =  0;
//	H_k.pData[Z_B1*H_COLS + X_W2] =  0;
//	H_k.pData[Z_B1*H_COLS + X_W3] =  0;
	H_k.pData[Z_B1*H_COLS + X_B] =  2*(q1*q3-q0*q2);
	
	//Row 2 d B(q2q3+q0q1)
	H_k.pData[Z_B2*H_COLS + X_Q0] =  2*q1*B;
	H_k.pData[Z_B2*H_COLS + X_Q1] =  2*q0*B;
	H_k.pData[Z_B2*H_COLS + X_Q2] =  2*q3*B;
	H_k.pData[Z_B2*H_COLS + X_Q3] =  2*q2*B;
//	H_k.pData[Z_B2*H_COLS + X_W1] =  0;
//	H_k.pData[Z_B2*H_COLS + X_W2] =  0;
//	H_k.pData[Z_B2*H_COLS + X_W3] =  0;
	H_k.pData[Z_B2*H_COLS + X_B] =  2*(q2*q3+q0*q1);
	
	//Row 3 d B(q0^2 - q1^2 + q2^2 + q3^2)/2
	H_k.pData[Z_B3*H_COLS + X_Q0] =  2*q0*B;
	H_k.pData[Z_B3*H_COLS + X_Q1] = -2*q1*B;
	H_k.pData[Z_B3*H_COLS + X_Q2] = -2*q2*B;
	H_k.pData[Z_B3*H_COLS + X_Q3] =  2*q3*B;
//  H_k.pData[Z_B3*H_COLS + X_W1] =  0;
//  H_k.pData[Z_B3*H_COLS + X_W2] =  0;
//  H_k.pData[Z_B3*H_COLS + X_W3] =  0;
	H_k.pData[Z_B3*H_COLS + X_B] =  (q0*q0 - q1*q1 - q2*q2 + q3*q3);
	
	/*Constant Part*/
	
	//Row 4 d Omega1
//	H_k.pData[Z_W1*H_COLS + X_Q0] =  0;
//	H_k.pData[Z_W1*H_COLS + X_Q1] =  0;
//	H_k.pData[Z_W1*H_COLS + X_Q2] =  0;
//	H_k.pData[Z_W1*H_COLS + X_Q3] =  0;
//	H_k.pData[Z_W1*H_COLS + X_W1] =  0;
//	H_k.pData[Z_W1*H_COLS + X_W2] =  1;
//	H_k.pData[Z_W1*H_COLS + X_W3] =  0;
//	H_k.pData[Z_W1*H_COLS + X_B]  =  0;
	
	//Row 5 d Omega2
//	H_k.pData[Z_W2*H_COLS + X_Q0] =  0;
//	H_k.pData[Z_W2*H_COLS + X_Q1] =  0;
//	H_k.pData[Z_W2*H_COLS + X_Q2] =  0;
//	H_k.pData[Z_W2*H_COLS + X_Q3] =  0;
//	H_k.pData[Z_W2*H_COLS + X_W1] =  0;
//	H_k.pData[Z_W2*H_COLS + X_W2] =  0;
//	H_k.pData[Z_W2*H_COLS + X_W3] =  1;
//	H_k.pData[Z_W2*H_COLS + X_B]  =  0;
	
	//Row 6 d Omega3 
//	H_k.pData[Z_W3*H_COLS + X_Q0] =  0;
//	H_k.pData[Z_W3*H_COLS + X_Q1] =  0;
//	H_k.pData[Z_W3*H_COLS + X_Q2] =  0;
//	H_k.pData[Z_W3*H_COLS + X_Q3] =  0;
//	H_k.pData[Z_W3*H_COLS + X_W1] =  0;
//	H_k.pData[Z_W3*H_COLS + X_W2] =  0;
//	H_k.pData[Z_W3*H_COLS + X_W3] =  0;
//	H_k.pData[Z_W3*H_COLS + X_B]  =  1;


//Compute transpose
array_transpose(&H_k, &H_k_t);
}

/*****************************************************/

void jacobianAph1(data_t q0,
									data_t q1,
									data_t q2,
									data_t q3,
									data_t w1,
									data_t w2,
									data_t w3, 
									data_t B)

{	
	//Row 1 d  1/2*(-Omega1*q1 - Omega2*q2 - Omega3*q3)
//	A_k.pData[X_Q0*A_COLS + X_Q0] = 	1;
	A_k.pData[X_Q0*A_COLS + X_Q1] =  -w1/2*sampleTime;
	A_k.pData[X_Q0*A_COLS + X_Q2] =  -w2/2*sampleTime;
	A_k.pData[X_Q0*A_COLS + X_Q3] =  -w3/2*sampleTime;
//	A_k.pData[X_Q0*A_COLS + X_W1] =  0;
//	A_k.pData[X_Q0*A_COLS + X_W2] =  0;
//	A_k.pData[X_Q0*A_COLS + X_W3] =  0;
//	A_k.pData[X_Q0*A_COLS + X_B]  =  0;
	
	//Row 2 d 1/2*(Omega1*q0 + Omega3*q2 - Omega2*q3)
	A_k.pData[X_Q1*A_COLS + X_Q0] =  w1/2*sampleTime;
//	A_k.pData[X_Q1*A_COLS + X_Q1] =  1;
	A_k.pData[X_Q1*A_COLS + X_Q2] =  w3/2*sampleTime;
	A_k.pData[X_Q1*A_COLS + X_Q3] = -w2/2*sampleTime;
//	A_k.pData[X_Q1*A_COLS + X_W1] =  0;
//	A_k.pData[X_Q1*A_COLS + X_W2] =  0;
//	A_k.pData[X_Q1*A_COLS + X_W3] =  0;
//	A_k.pData[X_Q1*A_COLS + X_B]  =  0;
	
	//Row 3 d 1/2*(Omega2*q0 - Omega3*q1 + Omega1*q3)
	A_k.pData[X_Q2*A_COLS + X_Q0] =  w2/2*sampleTime;
	A_k.pData[X_Q2*A_COLS + X_Q1] = -w3/2*sampleTime;
//	A_k.pData[X_Q2*A_COLS + X_Q2] =  1;
	A_k.pData[X_Q2*A_COLS + X_Q3] =  w1/2*sampleTime;
//	A_k.pData[X_Q2*A_COLS + X_W1] =  0;
//	A_k.pData[X_Q2*A_COLS + X_W2] =  0;
//	A_k.pData[X_Q2*A_COLS + X_W3] =  0;
//	A_k.pData[X_Q2*A_COLS + X_B]  =  0;
	
	//Row 4 d 1/2*(Omega3*q0 + Omega2*q1 - Omega1*q2)
	A_k.pData[X_Q3*A_COLS + X_Q0] =  w3/2*sampleTime;
	A_k.pData[X_Q3*A_COLS + X_Q1] =  w2/2*sampleTime;
	A_k.pData[X_Q3*A_COLS + X_Q2] = -w1/2*sampleTime;
//	A_k.pData[X_Q3*A_COLS + X_Q3] =  1;
//	A_k.pData[X_Q3*A_COLS + X_W1] =  0;
//	A_k.pData[X_Q3*A_COLS + X_W2] =  0;
//	A_k.pData[X_Q3*A_COLS + X_W3] =  0;
//	A_k.pData[X_Q3*A_COLS + X_B]  =  0;
	
	//Row 5 d ((I2-I3)/I1)*Omega1*Omega2
//	A_k.pData[X_W1*A_COLS + X_Q0] =  0;
//	A_k.pData[X_W1*A_COLS + X_Q1] =  0;
//	A_k.pData[X_W1*A_COLS + X_Q2] =  0;
//	A_k.pData[X_W1*A_COLS + X_Q3] =  0;
//	A_k.pData[X_W1*A_COLS + X_W1] =  1;
	A_k.pData[X_W1*A_COLS + X_W2] =  w3*(I2-I3)/I1*sampleTime;
	A_k.pData[X_W1*A_COLS + X_W3] =  w2*(I2-I3)/I1*sampleTime;
//	A_k.pData[X_W1*A_COLS + X_B]  =  0;
	
	//Row 6 d ((I3-I1)/I2)*Omega1*Omega3
//	A_k.pData[X_W2*A_COLS + X_Q0] =  0;
//	A_k.pData[X_W2*A_COLS + X_Q1] =  0;
//	A_k.pData[X_W2*A_COLS + X_Q2] =  0;
//	A_k.pData[X_W2*A_COLS + X_Q3] =  0;
	A_k.pData[X_W2*A_COLS + X_W1] =  w3*(I3-I1)/I2*sampleTime;
//	A_k.pData[X_W2*A_COLS + X_W2] =  1;
	A_k.pData[X_W2*A_COLS + X_W3] =  w1*(I3-I1)/I2*sampleTime;
//	A_k.pData[X_W2*A_COLS + X_B]  =  0;
	
	//Row 7 d ((I1-I2)/I3)*Omega1*Omega2
//	A_k.pData[X_W3*A_COLS + X_Q0] =  0;
//	A_k.pData[X_W3*A_COLS + X_Q1] =  0;
//	A_k.pData[X_W3*A_COLS + X_Q2] =  0;
//	A_k.pData[X_W3*A_COLS + X_Q3] =  0;
	A_k.pData[X_W3*A_COLS + X_W1] =  w2*(I1-I2)/I3*sampleTime;
	A_k.pData[X_W3*A_COLS + X_W2] =  w1*(I1-I2)/I3*sampleTime;
//	A_k.pData[X_W3*A_COLS + X_W3] =  1;
//	A_k.pData[X_W3*A_COLS + X_B]  =  0;
	
	//Row 7 0
//	A_k.pData[X_B*A_COLS + X_Q0] =  0;
//	A_k.pData[X_B*A_COLS + X_Q1] =  0;
//	A_k.pData[X_B*A_COLS + X_Q2] =  0;
//	A_k.pData[X_B*A_COLS + X_Q3] =  0;
//	A_k.pData[X_B*A_COLS + X_W1] =  0;
//	A_k.pData[X_B*A_COLS + X_W2] =  0;
//	A_k.pData[X_B*A_COLS + X_W3] =  0;
//	A_k.pData[X_B*A_COLS + X_B]  =  1;

//Compute transpose
array_transpose(&A_k, &A_k_t);	
}

/*******************************************************************/

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
									data_t T3)
{
	X_ap_k.pData[X_Q0] =  q0 + sampleTime*(-(w1 * q1) - (w2 * q2) - (w3 * q3))/2;
	X_ap_k.pData[X_Q1] =  q1 + sampleTime*((w1 * q0) + (w3 * q2) - (w2 * q3))/2;
	X_ap_k.pData[X_Q2] =  q2 + sampleTime*((w2 * q0) - (w3 * q1) + (w1 * q3))/2;
	X_ap_k.pData[X_Q3] =  q3 + sampleTime*((w3 * q0) + (w2 * q1) - (w1 * q2))/2;
	X_ap_k.pData[X_W1] =  w1 + (sampleTime/I1) * ((I2 - I3)*w2*w3 + T1);
	X_ap_k.pData[X_W2] =  w2 + (sampleTime/I2) * ((I3 - I1)*w1*w3 + T2);
	X_ap_k.pData[X_W3] =  w3 + (sampleTime/I3) * ((I1 - I2)*w1*w2 + T3);
	X_ap_k.pData[X_B] = B; 
}

/*******************************************************************/

void kalmanInv	 (data_t q0,
									data_t q1,
									data_t q2,
									data_t q3,
									data_t w1,
									data_t w2,
									data_t w3, 
									data_t B)
{
	Y_k.pData[Z_B1] = Z_k.pData[Z_B1] - B*(q1*q3-q0*q2);
	Y_k.pData[Z_B2] = Z_k.pData[Z_B2] - B*(q2*q3+q0*q1);
	Y_k.pData[Z_B3] = Z_k.pData[Z_B3] - (q2*q2 + q3*q3)/2;
	Y_k.pData[Z_W1] = Z_k.pData[Z_W1] - w1;
	Y_k.pData[Z_W2] = Z_k.pData[Z_W2] - w2;
	Y_k.pData[Z_W3] = Z_k.pData[Z_W3] - w3;
}

//****************************************************************** 

void kalmanInit(void)
{
 //Matrix initialization
	array_init(&X_k, X_ROWS, X_COLS, X_k_data);							//State vector
	array_init(&X_t1_k, X_T1_ROWS, X_T1_COLS, X_t1_k_data);	//State temp 1 vector
	array_init(&X_t2_k, X_T2_ROWS, X_T2_COLS, X_t2_k_data);	//State temp 2 vector
	array_init(&X_ap_k, X_AP_ROWS, X_AP_COLS, X_ap_k_data);	//State vector apriori estimate
	array_init(&Z_k, Z_ROWS, Z_COLS, Z_k_data);							//Measurement vector
	array_init(&Y_k, Y_ROWS, Y_COLS, Y_k_data);							//Measurement innovation
	array_init(&Qph1_k, QPH1_ROWS, QPH1_COLS, Qph1_k_data);	//System noise covariance matrix of phase 1
	array_init(&Rph1_k, RPH1_ROWS, RPH1_COLS, Rph1_k_data);	//Measurement noise covariance matrix of phase 1
	array_init(&Qph2_k, QPH2_ROWS, QPH2_COLS, Qph2_k_data);	//System noise covariance matrix of phase 2
	array_init(&Rph2_k, RPH2_ROWS, RPH2_COLS, Rph2_k_data); //Measurement noise covariance matrix of phase 2
	array_init(&A_k, A_ROWS, A_COLS, A_k_data);							//Global jacobian of f
	array_init(&A_k_t, A_COLS, A_ROWS, A_k_t_data);					//Global jacobian of h transpositnion
	array_init(&H_k, H_ROWS, H_COLS, H_k_data);							//Global jacobian of h
	array_init(&H_k_t, H_COLS, H_ROWS, H_k_t_data);					//Global jacobian of h transpositnion
	array_init(&P_k, P_ROWS, P_COLS, P_k_data);							//Global state pretiction covariance matrix
	array_init(&Kg_k, KG_ROWS, KG_COLS, Kg_k_data);					//Global kalman gain
	array_init(&S_k, S_ROWS, S_COLS, S_k_data);							//Scaling matrix
	array_init(&Temp_1, TEMP_1_ROWS, TEMP_1_COLS, Temp_1_data);		  //Temp matrix 1
	array_init(&Temp_2, TEMP_2_ROWS, TEMP_2_COLS, Temp_2_data);		  //Temp matrix 2
	array_init(&Temp_3, TEMP_3_ROWS, TEMP_3_COLS, Temp_3_data);		  //Temp matrix 3
	array_init(&Temp_4, TEMP_4_ROWS, TEMP_4_COLS, Temp_4_data);		  //Temp matrix 4 [NOI*DIM]
	array_init(&Temp_5, TEMP_5_ROWS, TEMP_5_COLS, Temp_5_data);		  //Temp matrix 5 [NOI*NOI]
	array_init(&Temp_6, TEMP_6_ROWS, TEMP_6_COLS, Temp_6_data);		  //Temp matrix 6 [NOI*NOI]
	array_init(&Temp_7, TEMP_7_ROWS, TEMP_7_COLS, Temp_7_data);		  //Temp matrix 7 [NOI*NOI]
	
	//Copy initial vaules to transposition matrices
	array_transpose(&H_k, &H_k_t);													//Global jacobian of h transpositnion
	array_transpose(&A_k, &A_k_t);													//Global jacobian of f transpositnion
	array_transpose(&W_k, &W_k_t);													//Wk transpositnion
	array_transpose(&V_k, &V_k_t);													//Vk transpositnion
	
	//State Vector initialization 
	kalmanInitializeStateVec();

}


//******************************************************************

void kalmanGetSensors(void)
{
	Z_k.pData[Z_B1] = 0;   //B1
	Z_k.pData[Z_B2] = 0;   //B2
	Z_k.pData[Z_B3] = 0;   //B3
	Z_k.pData[Z_W1] = 0;   //w1
	Z_k.pData[Z_W2] = 0;   //w2
	Z_k.pData[Z_W3] = 0;   //w3
}



//*********************************************************************

void kalmanInitializeStateVec(void)
{
	
	//Request data from sensors 
	kalmanGetSensors();
	
	//X_B - is used in further calculations
	//q0k = sqrt((B3 +Bk)/(2*Bk));
	sqrtf(Z_k.pData[Z_B1]*Z_k.pData[Z_B1] +
			 Z_k.pData[Z_B2]*Z_k.pData[Z_B2] +
			 Z_k.pData[Z_B3]*Z_k.pData[Z_B3],
			&(X_k.pData[X_B]));
	
	//Q0
	//Q0k = sqrt((B3 +Bk)/(2*Bk))
	sqrtf(Z_k.pData[Z_B3] + X_k.pData[X_B],
			&(X_k.pData[X_Q0]));
	X_k.pData[X_Q0] /= 2*X_k.pData[X_B];
	
	//Q1
	//q1K =B2 * sqrt(1/(2*Bk * (Bk + B3)))
	sqrtf(1/(2*(X_k.pData[X_B])*(X_k.pData[X_B] + Z_k.pData[Z_B3])),
			&(X_k.pData[X_Q1]));
	X_k.pData[X_Q1] *= Z_k.pData[Z_B2];
	
	//Q2
	//q2k = -B1 * sqrt(1/(2*Bk * (Bk + B3)))
	sqrtf(1/(2*X_k.pData[X_B]*(X_k.pData[X_B] + Z_k.pData[Z_B3])),
			&(X_k.pData[X_Q2]));
	X_k.pData[X_Q2] *= -Z_k.pData[Z_B1];
	
	//Q3
	//q3k = 0
	X_k.pData[X_Q3] = 0; 
	
	//W1
	//w1k = W1 
	X_k.pData[X_W1] = Z_k.pData[Z_W1];
	
	//W2
	//w2k = W2
	X_k.pData[X_W2] = Z_k.pData[Z_W2];
	
	//W3
	//w3k = W3
	X_k.pData[X_W3] = Z_k.pData[Z_W3];
}

//*********************************************************************
void kalmanStep(void)
{
	//Get state vector apriori estimate ---------------------------------
	kalmanF(X_k.pData[X_Q0],
					X_k.pData[X_Q1],
					X_k.pData[X_Q2],
					X_k.pData[X_Q3],
					X_k.pData[X_W1],
					X_k.pData[X_W2],
					X_k.pData[X_W3],
					X_k.pData[X_B],
					0,
					0,
					0);
	
	//Get F jacobian ----------------------------------------------------
	jacobianAph1(X_k.pData[X_Q0],
							 X_k.pData[X_Q1],
							 X_k.pData[X_Q2],
							 X_k.pData[X_Q3],
							 X_k.pData[X_W1],
							 X_k.pData[X_W2],
							 X_k.pData[X_W3],
							 X_k.pData[X_B]);
							 
	//Get P estimate apriori --------------------------------------------
	//Pkm = Ak * Pk * (Ak') + Wk * Qph1 * Wk'
	array_mult(&A_k, 		&P_k, 		&Temp_1); //  Temp_1 = (Ak * Pk)    	 [n*n]*[n*n] = [n*n] : Temp1
	array_mult(&Temp_1, &A_k_t,	  &Temp_2); //  Temp_2 = (Temp_1*Ak') 	 [n*n]*[n*n] = [n*n] : Temp2 - first part of the sum
	array_mult(&W_k, 		&Qph1_k,	&Temp_1); //	Temp_1 = (Wk * Qph1) 		 [n*n]*[n*n] = [n*n] : Temp1
	array_mult(&Temp_1, &Qph1_k,	&Temp_3); //	Temp_3 = (Temp_1*Qph1')  [n*n]*[n*n] = [n*n] : Temp3 - second part of the sum
	array_add(&Temp_1, &Temp_2,   &P_ap_k); //  P_ap_k = Temp_1 + Temp_3 [n*n] + [n*n]			 : P_ap_k
	
	
	//Get H jacobian ----------------------------------------------------
	jacobianHph1(X_k.pData[X_Q0],
							 X_k.pData[X_Q1],
							 X_k.pData[X_Q2],
							 X_k.pData[X_Q3],
							 X_k.pData[X_W1],
							 X_k.pData[X_W2],
							 X_k.pData[X_W3],
							 X_k.pData[X_B]);
	


	//Computing Kalman gain----------------------------------------------
	//(Hk * Pkm * (Hk') + Vk * Rph1 * (Vk')) = K1
	array_mult(&H_k, 		&P_ap_k, 		&Temp_4);  //  Temp_4 = (Hk * P_ap_k)[m*n]*[n*n] = [m*n] : Temp4 
	array_mult(&Temp_4, &H_k_t,			&Temp_5);  //  Temp_5 = (Temp_4* H_k)[m*n]*[n*m] = [m*m] : Temp5 -  first part of the sum
	array_mult(&V_k, 		&Rph1_k, 		&Temp_6);  //  Temp_6 = (Vk * Rph1_k)[m*n]*[n*n] = [m*m] : Temp6
	array_mult(&Temp_6, &V_k_t, 		&Temp_7);  //  Temp_7 = (Temp_5 * Vk')[m*n]*[n*n] = [m*m]: Temp7 -  second part of the sum
	 array_add(&Temp_5, &Temp_7, 		&Temp_6);  //  Temp_6 = Temp_5 + Temp_7 [m*m]+[m*m] = [m*m]
	
	// K2 = S*K1*S 
	array_mult(&S_k, 		&Temp_6, 		&Temp_5);  //  Temp_5 = S*Temp_6 [m*m]*[m*m] = [m*m]
	array_mult(&Temp_5, &S_k, 			&Temp_6);  //	 Temp_6 = Temp_5*S [m*m]*[m*m] = [m*m]
	
	// K3 = S*inv(K2)*S 
	array_inv(&Temp_6, &Temp_5);							 //	 Temp_5 = Temp_6^-1
	array_mult(&S_k, 		&Temp_5, 		&Temp_6);  //  Temp_6 = S*Temp_5 [m*m]*[m*m] = [m*m]
	array_mult(&Temp_6, &S_k, 			&Temp_5);  //	 Temp_5 = Temp_6*S [m*m]*[m*m] = [m*m]
	
	// Kg_k = Pkm * (Hk') * K3
	array_mult(&P_ap_k, &H_k_t, 		&Temp_6);  //  Temp_6 = P_ap_k*H_k'  [n*n]*[n*m] = [n*m]
	array_mult(&Temp_6, &Temp_5, 		&Kg_k);    //	 Kg_k = Temp_6*Temp_5[n*m]*[m*m] = [n*m]
	
	
	
	//Update state vector---------------------------------------------
	//xk = xkm + Kk * (zk - h(q0k, q1k, q2k, q3k, w1k, w2k, w3k, Bk));
	//                 ______________________________________________ - innovation - Y_k vecotr
	kalmanGetSensors();				 //Update Z_k
	kalmanInv(X_k.pData[X_Q0], //Update Y_k
						X_k.pData[X_Q1],
						X_k.pData[X_Q2],
						X_k.pData[X_Q3],
						X_k.pData[X_W1],
						X_k.pData[X_W2],
						X_k.pData[X_W3],
						X_k.pData[X_B]);
						
	array_mult(&Kg_k,  &Y_k, 		&X_t1_k);      //	 X_t1_k = Kg*Y_k[n*m]*[m*1] = [n*1]
	X_t2_k.pData[X_T2_Q0] = X_k.pData[X_Q0];   //  Create copy of state vector in temp 2
	X_t2_k.pData[X_T2_Q1] = X_k.pData[X_Q1];
	X_t2_k.pData[X_T2_Q2] = X_k.pData[X_Q2];
	X_t2_k.pData[X_T2_Q3] = X_k.pData[X_Q3];
	X_t2_k.pData[X_T2_W1] = X_k.pData[X_W1];
	X_t2_k.pData[X_T2_W2] = X_k.pData[X_W2];
	X_t2_k.pData[X_T2_W3] = X_k.pData[X_W3];
	X_t2_k.pData[X_T2_B] = X_k.pData[X_B];
	
	array_add(&X_t1_k,  &X_t2_k, 		&X_k); 		// X_k = X_t1_k + X_t2_k [n*1] + [n*1] = [n*1]
	
	

	//Update covariance matrix-------------------------------------- 
	//Pk = Pkm- Kk*Hk* Pkm;
	array_mult(&Kg_k, 		&H_k, 		&Temp_1);  //  Temp_1 = (Kg_k*H_k)			[n*m]*[m*n] = [n*n]
	array_mult(&Temp_1, 	&P_ap_k, 	&Temp_2);  //  Temp_2 = (Temp_1*P_ap_k) [n*n]*[n*n] = [n*n]
	array_sub( &P_ap_k, 	&Temp_2,	&P_k);		 //	 Pk = P_ap_k - Temp_2 		[n*m]-[n*n] = [n*n]
}
