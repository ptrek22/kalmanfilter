#include "uart.h"
#include "kalman.h"

#define _NOFBITS_ (sizeof(data_t))
	
union float2char
{
	data_t f;
	char b[_NOFBITS_];
};



void sendInParts(data_t toConvert);
data_t reshapeParts(char* toConvert);
void sendVec(data_t* p, int size);
void recVec(data_t* p, int size);
void transmitMatrix(array_t* p);
