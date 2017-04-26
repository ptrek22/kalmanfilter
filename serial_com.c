#include "serial_com.h"

union float2char f2c;

void sendInParts(data_t toConvert)
{
	f2c.f = toConvert;
	for (int i = 0; i < _NOFBITS_;  i++)
		uart_putchar(UART1, (char)(f2c.b[i]));
}

data_t reshapeParts(char* toConvert)
{
	for (int i = 0; i < _NOFBITS_;  i++)
		f2c.b[i] = toConvert[i];
	return f2c.f;
}

void sendVec(data_t* p, int size)
{
	for (int i = 0; i< size; i++)
	{
		sendInParts(p[i]);
	}
}


void recVec(data_t* p, int size)
{
	char buf[_NOFBITS_];
	for (int i = 0; i< size; i++)
	{
		for(int j = 0; j < _NOFBITS_; j++)
			buf[j] = uart_getchar(UART1);
		p[i] = reshapeParts(buf);
	}
}

void transmitMatrix(array_t* p)
{
	uart_putchar(UART1, (char)(p->numRows));
	uart_putchar(UART1, (char)(p->numCols));
	sendVec(p->pData,(p->numRows)*(p->numCols));
	
}
