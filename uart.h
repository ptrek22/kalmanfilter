#include "MKL25Z4.h"                    // Device header
#include <stdio.h>

//----------------------------------
//UART PINS IN FRDM KL25Z

//UART 0
//RX - PTA14, PTB16, PTD6 (ALT 4)
//TX - PTA15, PTB17, PTD7 (ALT 4)


//UART 1
//RX - PTA18, PTC3, PTE1 (ALT 3)
//TX - PTA19, PTC4, PTE0 (ALT 3)

//UART 2
//RX - PTE22, (ALT4) PTD2, PTD4 (ALT 3)
//TX - PTE23, (ALT4) PTD3, PTD5 (ALT 3)




void uart_init(UART_MemMapPtr uartPtr, int sysclk, int baud);

int uart_putchar (UART_MemMapPtr channel, char ch);

char uart_getchar (UART_MemMapPtr channel);

void uart_enable_recieve_interrupt(UART_MemMapPtr uartPtr);

void uart_disable_recieve_interrupt(UART_MemMapPtr uartPtr);

void uart_enable_transfer_interrupt(UART_MemMapPtr uartPtr);

void uart_disable_transfer_interrupt(UART_MemMapPtr uartPtr);


