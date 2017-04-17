#define ARM_MATH_CM0PLUS
#include "kalman.h"

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

void kalmanInit()
{
 //Matrix initialization
	array_init(&X_k, X_ROWS, X_COLS, X_k_data);							//State vector
	array_init(&X_ap_k, X_AP_ROWS, X_AP_COLS, X_ap_k_data);				//State vector apriori estimate
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
	array_init(&Temp_1, TEMP_1_ROWS, TEMP_1_COLS, Temp_1_data);		  //Temp matrix 1
	array_init(&Temp_2, TEMP_2_ROWS, TEMP_2_COLS, Temp_2_data);		  //Temp matrix 2
	array_init(&Temp_3, TEMP_3_ROWS, TEMP_3_COLS, Temp_3_data);		  //Temp matrix 3
	array_init(&Temp_4, TEMP_4_ROWS, TEMP_4_COLS, Temp_4_data);		  //Temp matrix 4 [NOI*NOI]
	array_init(&Temp_5, TEMP_5_ROWS, TEMP_5_COLS, Temp_5_data);		  //Temp matrix 5 [NOI*NOI]
	array_init(&Temp_6, TEMP_6_ROWS, TEMP_6_COLS, Temp_6_data);		  //Temp matrix 4 [NOI*NOI]
	
	//Copy initial vaules to transposition matrices
	array_transpose(&H_k, &H_k_t);													//Global jacobian of h transpositnion
	array_transpose(&A_k, &A_k_t);													//Global jacobian of f transpositnion
	array_transpose(&W_k, &W_k_t);													//Global jacobian of f transpositnion
	
	//State Vector initialization 
	kalmanInitializeStateVec();

}


//******************************************************************

void kalmanGetOtputs()
{
	Z_k.pData[Z_B1] = 0;   //B1
	Z_k.pData[Z_B2] = 0;   //B2
	Z_k.pData[Z_B3] = 0;   //B3
	Z_k.pData[Z_W1] = 0;   //w1
	Z_k.pData[Z_W2] = 0;   //w2
	Z_k.pData[Z_W3] = 0;   //w3
}

//******************************************************************

void updateStateVec(
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
	X_k.pData[X_Q0] = 0;   //Q0
	X_k.pData[X_Q1] = 0;   //Q1
	X_k.pData[X_Q2] = 0;   //Q2
	X_k.pData[X_Q3] = 0;   //Q2
	X_k.pData[X_W1] = 0;   //w1
	X_k.pData[X_W2] = 0;   //w2
	X_k.pData[X_W3] = 0;   //w3
	X_k.pData[X_B] 	= 0;
}


//*********************************************************************

void kalmanInitializeStateVec()
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
void kalmanStep()
{
	//Get state vector apriori estimate
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
	
	//Get F jacobian
	jacobianAph1(X_k.pData[X_Q0],
							 X_k.pData[X_Q1],
							 X_k.pData[X_Q2],
							 X_k.pData[X_Q3],
							 X_k.pData[X_W1],
							 X_k.pData[X_W2],
							 X_k.pData[X_W3],
							 X_k.pData[X_B]);
							 
	//Get P estimate apriori	
	//Pkm = Ak * Pk * (Ak') + Wk * Qph1 * Wk'
	array_mult(&A_k, 		&P_k, 		&Temp_1); //  Temp_1 = (Ak * Pk)    	 [n*n]*[n*n] = [n*n]
	array_mult(&Temp_1, &A_k_t,	  &Temp_2); //  Temp_2 = (Temp_1*Ak') 	 [n*n]*[n*n] = [n*n] - first part of the sum
	array_mult(&W_k, 		&Qph1_k,	&Temp_1); //	Temp_1 = (Wk * Qph1) 		 [n*n]*[n*n] = [n*n]
	array_mult(&Temp_1, &Qph1_k,	&Temp_3); //	Temp_3 = (Temp_1*Qph1')  [n*n]*[n*n] = [n*n] - second part of the sum
	array_add(&Temp_1, &Temp_2,   &P_ap_k); //  P_ap_k = Temp_1 + Temp_3 [n*n] + [n*n]
	
	
	//Get H jacobian
	jacobianHph1(X_k.pData[X_Q0],
							 X_k.pData[X_Q1],
							 X_k.pData[X_Q2],
							 X_k.pData[X_Q3],
							 X_k.pData[X_W1],
							 X_k.pData[X_W2],
							 X_k.pData[X_W3],
							 X_k.pData[X_B]);
							 
	//Computing Kalman gain
	
}
