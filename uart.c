#include "uart.h"
#ifndef _UART_H_
#define _UART_H_

/*
Zródlem zegara dla UART1 i UART2 jest BUSCLOCK, jest on generowany przez podzial  czestotliwosci SYSTEM_CLOCK przez wartosc zawarta w rejestrze SIM_CLKDIV1
*/

//------------------------------------------------------------------------------------------------------------------------------

void uart_init(UART_MemMapPtr uartPtr, int sysclk, int baud)
{
	register uint32_t ubd;
	uint8_t temp;
	
	/* Doprowadzenie zegara do wybranego UART */
	//if(uartPtr == UART0_BASE_PTR)					//Tu wystepuje niezgodnosc typów, gdyz UART0 ma osobna strukture realizujaca dostep do rejestrów
	//	SIM_SCGC4 |= SIM_SCGC4_UART1_MASK;	//nie mniej powinno to dzialac, bo kolejnosc uzywanych tutaj rejestrów jest zgodna (UART0 ma rejestry dodatkowe)
	//a jednak nie...
  if(uartPtr == UART1_BASE_PTR)		//SIM - System Integration Module
		SIM_SCGC4 |= SIM_SCGC4_UART1_MASK;  //System Clock Gating Control Register, odpoiwadaja za doprowadzanie zegara
	else if(uartPtr == UART2_BASE_PTR)
		SIM_SCGC4 |= SIM_SCGC4_UART2_MASK;
	else return;


	/* Nie nalezy zmieniac ustawien gdy UART pracuje, nalezy go wylaczyc */
	UART_C2_REG(uartPtr) &= ~(UART_C2_TE_MASK | UART_C2_RE_MASK); //uartPtr -> C2

	/* Przywrócenie ustawien domyslnych 
		8 bit, brak kontroli parzystosci */
	UART_C1_REG(uartPtr) = 0;
		
	/* predkosc transferu */
	ubd = (uint16_t)((sysclk)/(baud*16));

	/* 4 najmlodsze bity rejestru BDH i rejestr BDL (w sumie 13 bitów - SBR) ustalaja podzial czestotliwosci zegarowania  UART (UART0 moze posiadac inne
		zródla zegara (rej SOPT2)) -  zegar wewnetrzny jest dodatkowo dzelony przez 16 */
		
	/*przechowanie starej wartosci BDH z wyjatkiem 4 najmlodszych bitów które tworza 4 najstarsze bity SBR */

	temp = UART_BDH_REG(uartPtr) & ~(UART_BDH_SBR_MASK);

	/*ustalenie predkosci transmisji*/

	UART_BDL_REG(uartPtr) = (uint8_t) UART_BDL_SBR_MASK & ubd;
	UART_BDH_REG(uartPtr) = (uint8_t) ((ubd >> 8) & UART_BDH_SBR_MASK) | temp ;

	/* Uruchomienie transmisji i odbioru */
	UART_C2_REG(uartPtr) |= (UART_C2_TE_MASK | UART_C2_RE_MASK); 

}

//------------------------------------------------------------------------------------------------------------------------------

int uart_putchar (UART_MemMapPtr channel, char ch)
{
	/*Transmit Data Register Empty Flag - ustawiana gdy bufor danych do wyslania UART_D_REG jest wolny */
	/*Czekaj az bedzie miejsce*/ 
	while(!(UART_S1_REG(channel) & UART_S1_TDRE_MASK)) continue;
	/* Wysylanie znaku */
	return UART_D_REG(channel) = (uint8_t)ch;
}

//------------------------------------------------------------------------------------------------------------------------------

char uart_getchar (UART_MemMapPtr channel)
{
	/* Receive Data Register Full Flag - ustawiana po odebraniu danych, zerowana po odczycie z  UART_D_REG*/
	while (!(UART_S1_REG(channel) & UART_S1_RDRF_MASK)) continue;
	/* Return the 8-bit data from the receiver */
	return UART_D_REG(channel);
}

//------------------------------------------------------------------------------------------------------------------------------

void uart_enable_recieve_interrupt(UART_MemMapPtr uartPtr)
{
	UART_C2_REG(uartPtr) |=  UART_C2_RIE_MASK; //Receiver Interrupt Enable
}

//------------------------------------------------------------------------------------------------------------------------------

void uart_disable_recieve_interrupt(UART_MemMapPtr uartPtr)
{
	UART_C2_REG(uartPtr) &=  ~UART_C2_RIE_MASK; //Receiver Interrupt Enable
}

//------------------------------------------------------------------------------------------------------------------------------

void uart_enable_transmit_interrupt(UART_MemMapPtr uartPtr)
{
	UART_C2_REG(uartPtr) |=  UART_C2_TCIE_MASK; //Transmission complete - alternatywnie TIE - zglasza przerwanie gdy bufor jest pusty 
}

//------------------------------------------------------------------------------------------------------------------------------

void uart_enable_disable_interrupt(UART_MemMapPtr uartPtr)
{
	UART_C2_REG(uartPtr) &=  ~UART_C2_TCIE_MASK; //Transmission complete - alternatywnie TIE - zglasza przerwanie gdy bufor jest pusty 
}

//------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------Copied from retarget.c of some example code -------------------------------

struct __FILE 
{
	int handle;
	/* Add whatever you need here */ 
};
FILE __stdout;
FILE __stdin;

//------------------------------------------------------------------------------------------------------------------------------

int fputc(int c, FILE *f)
{
	return uart_putchar(UART1_BASE_PTR, c);
}
//------------------------------------------------------------------------------------------------------------------------------
int getkey (void)
{
  return uart_getchar(UART1_BASE_PTR);
}
//------------------------------------------------------------------------------------------------------------------------------
int fgetc (FILE *f)
{
  return (getkey ());
}




#endif
