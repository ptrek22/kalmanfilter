#include "system_MKL25Z4.h"             // Keil::Device:Startup
#include "MKL25Z4.h"                    // Device header
//#include "kalman.h"
#include "serial_com.h"
#include "uart.h"
#include "stdio.h"

void uart_gpio_init(void)
{
	SIM_SCGC5 |= SIM_SCGC5_PORTC_MASK;
	PORTC_PCR3 |= (3<<8); //Alt 3
	PORTC_PCR4 |= (3<<8); //Alt 3
	uart_init(UART1_BASE_PTR, DEFAULT_SYSTEM_CLOCK/2, 9600);
}
	

int main()
{	
	uart_gpio_init();
	kalmanInit();

	recVec(Z_k.pData,_NOI_);
	recVec(X_k.pData,_DIM_);
	
	while(1)
	{
		kalmanStep();
		transmitMatrix(&X_k);
	}
//	transmitMatrix(&A_k);
//	transmitMatrix(&P_k);
//	transmitMatrix(&W_k);
//	transmitMatrix(&Qph1_k);
//	transmitMatrix(&P_ap_k);
//	transmitMatrix(&Temp_8);
//	transmitMatrix(&Temp_5);
//	transmitMatrix(&Kg_k);
	

}
