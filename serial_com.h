#ifndef _SERIAL_COM_H_
#define _SERIAL_COM_H_

#include "uart.h"
#define ARM_MATH_CM0PLUS
#include "arm_math.h"     

typedef float32_t data_t;
typedef arm_matrix_instance_f32  array_t ;
#define _NOFBITS_ (sizeof(data_t))
	
union float2char
{
	data_t f;
	char b[_NOFBITS_];
};

extern union float2char f2c;



void sendInParts(data_t toConvert);
data_t reshapeParts(char* toConvert);
void sendVec(data_t* p, int size);
void recVec(data_t* p, int size);
void transmitMatrix(array_t* p);

#endif