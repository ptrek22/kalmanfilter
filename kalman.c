#define ARM_MATH_CM0PLUS
#include "kalman.h"

void jacobianHph1 (array_t* dest,
									data_t q0,
									data_t q1,
									data_t q2,
									data_t q3,
									data_t w1,
									data_t w2,
									data_t w3, 
									data_t B)

{	//Row 1 d B(q1q3-q0q2)
	dest->pData[0*_DIM_ + 0] = -2*q2*B;
	dest->pData[0*_DIM_ + 1] =  2*q3*B;
	dest->pData[0*_DIM_ + 2] = -2*q0*B;
	dest->pData[0*_DIM_ + 3] =  2*q1*B;
//	dest->pData[0*_DIM_ + 4] =  0;
//	dest->pData[0*_DIM_ + 5] =  0;
//	dest->pData[0*_DIM_ + 6] =  0;
	dest->pData[0*_DIM_ + 7] =  2*(q1*q3-q0*q2);
	
	//Row 2 d B(q2q3+q0q1)
	dest->pData[1*_DIM_ + 0] =  2*q1*B;
	dest->pData[1*_DIM_ + 1] =  2*q0*B;
	dest->pData[1*_DIM_ + 2] =  2*q3*B;
	dest->pData[1*_DIM_ + 3] =  2*q2*B;
//	dest->pData[1*_DIM_ + 4] =  0;
//	dest->pData[1*_DIM_ + 5] =  0;
//	dest->pData[1*_DIM_ + 6] =  0;
	dest->pData[1*_DIM_ + 7] =  2*(q2*q3+q0*q1);
	
	//Row 3 d B(q0^2 - q1^2 + q2^2 + q3^2)/2
	dest->pData[2*_DIM_ + 0] =  2*q0*B;
	dest->pData[2*_DIM_ + 1] = -2*q1*B;
	dest->pData[2*_DIM_ + 2] = -2*q2*B;
	dest->pData[2*_DIM_ + 3] =  2*q3*B;
//  dest->pData[2*_DIM_ + 4] =  0;
//  dest->pData[2*_DIM_ + 5] =  0;
//  dest->pData[2*_DIM_ + 6] =  0;
	dest->pData[2*_DIM_ + 7] =  (q0*q0 - q1*q1 - q2*q2 + q3*q3);
	
	/*Constant Part*/
	
	//Row 4 d Omega1
//	dest->pData[3*_DIM_ + 0] =  0;
//	dest->pData[3*_DIM_ + 1] =  0;
//	dest->pData[3*_DIM_ + 2] =  0;
//	dest->pData[3*_DIM_ + 3] =  0;
//	dest->pData[3*_DIM_ + 4] =  0;
//	dest->pData[3*_DIM_ + 5] =  1;
//	dest->pData[3*_DIM_ + 6] =  0;
//	dest->pData[3*_DIM_ + 7] =  0;
	
	//Row 5 d Omega2
//	dest->pData[4*_DIM_ + 0] =  0;
//	dest->pData[4*_DIM_ + 1] =  0;
//	dest->pData[4*_DIM_ + 2] =  0;
//	dest->pData[4*_DIM_ + 3] =  0;
//	dest->pData[4*_DIM_ + 4] =  0;
//	dest->pData[4*_DIM_ + 5] =  0;
//	dest->pData[4*_DIM_ + 6] =  1;
//	dest->pData[4*_DIM_ + 7] =  0;
	
	//Row 6 d Omega3 
//	dest->pData[5*_DIM_ + 0] =  0;
//	dest->pData[5*_DIM_ + 1] =  0;
//	dest->pData[5*_DIM_ + 2] =  0;
//	dest->pData[5*_DIM_ + 3] =  0;
//	dest->pData[5*_DIM_ + 4] =  0;
//	dest->pData[5*_DIM_ + 5] =  0;
//	dest->pData[5*_DIM_ + 6] =  0;
//	dest->pData[5*_DIM_ + 7] =  1;
}

/*****************************************************/

void jacobianAph1 (array_t* dest,
									data_t q0,
									data_t q1,
									data_t q2,
									data_t q3,
									data_t w1,
									data_t w2,
									data_t w3, 
									data_t B)

{	
	//Row 1 d  1/2*(-Omega1*q1 - Omega2*q2 - Omega3*q3)
//	dest->pData[0*_DIM_ + 0] = 	1;
	dest->pData[0*_DIM_ + 1] =  -w1/2*sampleTime;
	dest->pData[0*_DIM_ + 2] =  -w2/2*sampleTime;
	dest->pData[0*_DIM_ + 3] =  -w3/2*sampleTime;
//	dest->pData[0*_DIM_ + 4] =  0;
//	dest->pData[0*_DIM_ + 5] =  0;
//	dest->pData[0*_DIM_ + 6] =  0;
//	dest->pData[0*_DIM_ + 7] =  0;
	
	//Row 2 d 1/2*(Omega1*q0 + Omega3*q2 - Omega2*q3)
	dest->pData[1*_DIM_ + 0] =  w1/2*sampleTime;
//	dest->pData[1*_DIM_ + 1] =  1;
	dest->pData[1*_DIM_ + 2] =  w3/2*sampleTime;
	dest->pData[1*_DIM_ + 3] = -w2/2*sampleTime;
//	dest->pData[1*_DIM_ + 4] =  0;
//	dest->pData[1*_DIM_ + 5] =  0;
//	dest->pData[1*_DIM_ + 6] =  0;
//	dest->pData[1*_DIM_ + 7] =  0;
	
	//Row 3 d 1/2*(Omega2*q0 - Omega3*q1 + Omega1*q3)
	dest->pData[2*_DIM_ + 0] =  w2/2*sampleTime;
	dest->pData[2*_DIM_ + 1] = -w3/2*sampleTime;
//	dest->pData[2*_DIM_ + 2] =  1;
	dest->pData[2*_DIM_ + 3] =  w1/2*sampleTime;
//	dest->pData[2*_DIM_ + 4] =  0;
//	dest->pData[2*_DIM_ + 5] =  0;
//	dest->pData[2*_DIM_ + 6] =  0;
//	dest->pData[2*_DIM_ + 7] =  0;
	
	//Row 4 d 1/2*(Omega3*q0 + Omega2*q1 - Omega1*q2)
	dest->pData[3*_DIM_ + 0] =  w3/2*sampleTime;
	dest->pData[3*_DIM_ + 1] =  w2/2*sampleTime;
	dest->pData[3*_DIM_ + 2] = -w1/2*sampleTime;
//	dest->pData[3*_DIM_ + 3] =  1;
//	dest->pData[3*_DIM_ + 4] =  0;
//	dest->pData[3*_DIM_ + 5] =  0;
//	dest->pData[3*_DIM_ + 6] =  0;
//	dest->pData[3*_DIM_ + 7] =  0;
	
	//Row 5 d ((I2-I3)/I1)*Omega1*Omega2
//	dest->pData[4*_DIM_ + 0] =  0;
//	dest->pData[4*_DIM_ + 1] =  0;
//	dest->pData[4*_DIM_ + 2] =  0;
//	dest->pData[4*_DIM_ + 3] =  0;
//	dest->pData[4*_DIM_ + 4] =  1;
	dest->pData[4*_DIM_ + 5] =  w3*(I2-I3)/I1*sampleTime;
	dest->pData[4*_DIM_ + 6] =  w2*(I2-I3)/I1*sampleTime;
//	dest->pData[4*_DIM_ + 7] =  0;
	
	//Row 6 d ((I3-I1)/I2)*Omega1*Omega3
//	dest->pData[5*_DIM_ + 0] =  0;
//	dest->pData[5*_DIM_ + 1] =  0;
//	dest->pData[5*_DIM_ + 2] =  0;
//	dest->pData[5*_DIM_ + 3] =  0;
	dest->pData[5*_DIM_ + 4] =  w3*(I3-I1)/I2*sampleTime;
//	dest->pData[5*_DIM_ + 5] =  1;
	dest->pData[5*_DIM_ + 6] =  w1*(I3-I1)/I2*sampleTime;
//	dest->pData[5*_DIM_ + 7] =  0;
	
	//Row 7 d ((I1-I2)/I3)*Omega1*Omega2
//	dest->pData[6*_DIM_ + 0] =  0;
//	dest->pData[6*_DIM_ + 1] =  0;
//	dest->pData[6*_DIM_ + 2] =  0;
//	dest->pData[6*_DIM_ + 3] =  0;
	dest->pData[6*_DIM_ + 4] =  w2*(I1-I2)/I3*sampleTime;
	dest->pData[6*_DIM_ + 5] =  w1*(I1-I2)/I3*sampleTime;
//	dest->pData[6*_DIM_ + 6] =  1;
//	dest->pData[6*_DIM_ + 7] =  0;
	
	//Row 7 0
//	dest->pData[7*_DIM_ + 0] =  0;
//	dest->pData[7*_DIM_ + 1] =  0;
//	dest->pData[7*_DIM_ + 2] =  0;
//	dest->pData[7*_DIM_ + 3] =  0;
//	dest->pData[7*_DIM_ + 4] =  0;
//	dest->pData[7*_DIM_ + 5] =  0;
//	dest->pData[7*_DIM_ + 6] =  0;
//	dest->pData[7*_DIM_ + 7] =  1;
}

/*******************************************************************/

void kalmanF			  (array_t* dest,
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
									data_t T3)
{
	dest->pData[0] =  q0 + sampleTime*(-(w1 * q1) - (w2 * q2) - (w3 * q3))/2;
	dest->pData[1] =  q1 + sampleTime*((w1 * q0) + (w3 * q2) - (w2 * q3))/2;
	dest->pData[2] =  q2 + sampleTime*((w2 * q0) - (w3 * q1) + (w1 * q3))/2;
	dest->pData[3] =  q3 + sampleTime*((w3 * q0) + (w2 * q1) - (w1 * q2))/2;
	dest->pData[4] =  w1 + (sampleTime/I1) * ((I2 - I3)*w2*w3 + T1);
	dest->pData[5] =  w2 + (sampleTime/I2) * ((I3 - I1)*w1*w3 + T2);
	dest->pData[6] =  w3 + (sampleTime/I3) * ((I1 - I2)*w1*w2 + T3);
	dest->pData[7] = B; 
}

/*******************************************************************/

void kalmanH 			(array_t* dest,
									data_t q0,
									data_t q1,
									data_t q2,
									data_t q3,
									data_t w1,
									data_t w2,
									data_t w3, 
									data_t B)
{
	dest->pData[0] = B*(q1*q3-q0*q2);
	dest->pData[1] = B*(q2*q3+q0*q1);
	dest->pData[2] = B*(q0*q0 - q1*q1 - q2*q2 + q3*q3)/2;
	dest->pData[4] = w1;
	dest->pData[5] = w2;
	dest->pData[6] = w3;
}

//****************************************************************** 

void kalmanInit()
{
 //Matrix initialization
	array_init(&X_k, X_ROWS, X_COLS, X_k_data);							//State vector
	array_init(&Z_k, Z_ROWS, Z_COLS, Z_k_data);							//Measurement vector
	array_init(&Qph1_k, QPH1_ROWS, QPH1_COLS, Qph1_k_data);	//System noise covariance matrix of phase 1
	array_init(&Rph1_k, RPH1_ROWS, RPH1_COLS, Rph1_k_data);	//Measurement noise covariance matrix of phase 1
	array_init(&Qph2_k, QPH2_ROWS, QPH2_COLS, Qph2_k_data);	//System noise covariance matrix of phase 2
	array_init(&Rph2_k, RPH2_ROWS, RPH2_COLS, Rph2_k_data); //Measurement noise covariance matrix of phase 2
	array_init(&A_k, A_ROWS, A_COLS, A_k_data);							//Global jacobian of f
	array_init(&H_k, H_ROWS, H_COLS, H_k_data);							//Global jacobian of h 
	array_init(&P_k, P_ROWS, P_COLS, P_k_data);							//Global state pretiction covariance matrix
}


//******************************************************************

void getSensors()
{
	Z_k.pData[0] = 0;   //B1
	Z_k.pData[1] = 0;   //B2
	Z_k.pData[2] = 0;   //B3
	Z_k.pData[3] = 0;   //w1
	Z_k.pData[4] = 0;   //w2
	Z_k.pData[5] = 0;   //w3
}


