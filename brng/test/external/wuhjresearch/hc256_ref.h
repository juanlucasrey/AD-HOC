/* This program gives the reference implementation of HC-256 for 32-bit platform

   The docuement of HC-256 is available at
   1) Hongjun Wu. ``A New Stream Cipher HC-256.'' Fast Software Encryption -- FSE 2004, LNCS 3017, pp. 226-244, Springer-Verlag 2004.
   2) eSTREAM website:  http://www.ecrypt.eu.org/stream/hcp3.html

   -----------------------------------------
   Performance of this non-optimized implementation:

   Microprocessor: Intel CORE 2 processor (Core 2 Duo Mobile P9400 2.53GHz) 
   Operating System: 32-bit Debian 5.0 (Linux kernel 2.6.26-2-686)
   Speed of encrypting long message:
   1) 7.0 cycle/byte   compiler: Intel C++ compiler 11.1   compilation option: icc -O2 
   2) 6.5 cycles/byte  compiler: gcc 4.3.2                 compilation option: gcc -O3  

   Microprocessor: Intel CORE 2 processor (Core 2 Quad Q6600 2.4GHz) 
   Operating System: 32-bit Windows Vista Business
   Speed of encrypting long message:
   1) 7.3 cycles/byte  compiler: Intel C++ compiler 11.1    compilation option: icl /O2
   2) 8.9 cycles/byte  compiler: Microsoft Visual C++ 2008  compilation option: release

   ------------------------------------------
   Written by: Hongjun Wu
   Last Modified: December 15, 2009  
*/

#include <string.h>

typedef unsigned char uint8;
typedef unsigned long long uint64;

/*for LP64, "int" is 32-bit integer, while "long" is 64-bit integer*/ 
#if defined(_LP64) 
    typedef unsigned int uint32; 
#else 
    typedef unsigned long uint32; 
#endif 

typedef struct { 
      uint32 P[1024];
      uint32 Q[1024];
      uint32 counter2048;     /*counter2048 = i mod 2048 */
      uint32 keystreamword;   /*a 32-bit keystream word*/
} HC256_State; 

#define ROTR32(x,n)   ( ((x) >> (n))  | ((x) << (32 - (n))) )

#define f1(x)    (ROTR32((x),7) ^ ROTR32((x),18) ^ ((x) >> 3))
#define f2(x)    (ROTR32((x),17) ^ ROTR32((x),19) ^ ((x) >> 10))

#define g1(x,y)  ((ROTR32((x),10)^ROTR32((y),23))+state->Q[((x)^(y))&0x3ff])
#define g2(x,y)  ((ROTR32((x),10)^ROTR32((y),23))+state->P[((x)^(y))&0x3ff])

/*function h1*/ 
uint32 h1(HC256_State *state, uint32 u) {
      uint32 tem; 			
      uint8  a,b,c,d;			
      a = (uint8) ((u));		
      b = (uint8) ((u) >> 8);	
      c = (uint8) ((u) >> 16);	
      d = (uint8) ((u) >> 24);	
      tem = state->Q[a]+state->Q[256+b]+state->Q[512+c]+state->Q[768+d];
      return (tem);
}

/*function h2*/
uint32 h2(HC256_State *state, uint32 u) {
      uint32 tem; 			
      uint8  a,b,c,d;			
      a = (uint8) ((u));		
      b = (uint8) ((u) >> 8);	
      c = (uint8) ((u) >> 16);	
      d = (uint8) ((u) >> 24);	
      tem = state->P[a]+state->P[256+b]+state->P[512+c]+state->P[768+d];
      return (tem);
}

/* one step of HC-256: 
   state is updated;
   a 32-bit keystream word is generated and stored in "state->keystreamword";
*/
void OneStep(HC256_State *state) 
{
      uint32 i,i3, i10, i12, i1023;

      i   = state->counter2048 & 0x3ff;
      i3  = (i - 3) & 0x3ff;
      i10 = (i - 10) & 0x3ff;
      i12 = (i - 12) & 0x3ff;
      i1023 = (i - 1023) & 0x3ff;

      if (state->counter2048 < 1024) {
            state->P[i] = state->P[i] + state->P[i10] + g1(state->P[i3],state->P[i1023]);
            state->keystreamword = h1(state,state->P[i12]) ^ state->P[i];
      }
      else {
            state->Q[i] = state->Q[i] + state->Q[i10] + g2(state->Q[i3],state->Q[i1023]);
            state->keystreamword = h2(state,state->Q[i12]) ^ state->Q[i];
      }
      state->counter2048 = (state->counter2048+1) & 0x7ff;
}

/*this function initialize the state using 256-bit key and 256-bit IV*/
void Initialization(HC256_State *state, uint8 *key, uint8 *iv) 
{
      uint32 W[2560],i;

      /*expand the key and iv into the state*/
      for (i = 0; i < 8;  i++) W[i] = ((uint32*)key)[i];
      for (i = 8; i < 16; i++) W[i] = ((uint32*)iv)[i-8]; 

      for (i = 16; i < 2560; i++) W[i] = f2(W[i-2]) + W[i-7] + f1(W[i-15]) + W[i-16]+i; 

      for (i = 0; i < 1024; i++)  state->P[i] = W[i+512];
      for (i = 0; i < 1024; i++)  state->Q[i] = W[i+1536];

      state->counter2048 = 0;

      /*run the cipher for 4096 steps without generating output*/
      for (i = 0; i < 4096; i++)  OneStep(state);
}	


/* this function encrypts a message*/
void EncryptMessage(HC256_State *state, uint8 *message, uint8 *ciphertext, uint64 msglength)
{
      uint64 i;
      uint32 j;
      
      /*encrypt a message, each time 4 bytes are encrypted*/
      for (i = 0; (i+4) <= msglength; i += 4, message += 4, ciphertext += 4) {
            /*generate 32-bit keystream and store it in state.keystreamword*/
            OneStep(state);   
            /*encrypt 32 bits of the message*/
            ((uint32*)ciphertext)[0] = ((uint32*)message)[0] ^ state->keystreamword;   
      }
      /*encrypt the last message block if the message length is not multiple of 4 bytes*/
      if ((msglength & 3) != 0) {
            OneStep(state);
            for (j = 0; j < (msglength & 3); j++) {
                  *(ciphertext+j) = *(message+j) ^ ((uint8*)&state->keystreamword)[j];
            }
      }    
}


/* this function encrypts a message,
   there are four inputs to this function: a 256-bit key, a 256-bit iv, a message, the message length in bytes
   one output from this function: ciphertext
*/
void HC256(uint8 *key, uint8 *iv, uint8 *message, uint8 *ciphertext, uint64 msglength)
{
      HC256_State state;

      /*initializing the state*/
      Initialization(&state,key,iv);

      /*encrypt a message*/
      EncryptMessage(&state,message,ciphertext,msglength);
}
