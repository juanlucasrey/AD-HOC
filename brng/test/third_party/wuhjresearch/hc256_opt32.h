/* This program gives the optimized implementation of HC-256 for 32-bit platform

   The docuement of HC-256 is available at:
   1) Hongjun Wu. ``A New Stream Cipher HC-256.'' Fast Software Encryption -- FSE 2004, LNCS 3017, pp. 226-244, Springer-Verlag 2004.
   2) eSTREAM website:  http://www.ecrypt.eu.org/stream/hcp3.html

   -----------------------------------------
   Performance:

   Microprocessor: Intel CORE 2 processor (Core 2 Duo Mobile P9400 2.53GHz) 
   Operating System: 32-bit Debian 5.0 (Linux kernel 2.6.26-2-686)
   Speed of encrypting long message:
   1) 3.3 cycle/byte   compiler: Intel C++ compiler 11.1   compilation option: icc -O2 
   2) 4.9 cycles/byte  compiler: gcc 4.3.2                 compilation option: gcc -O3  

   Microprocessor: Intel CORE 2 processor (Core 2 Quad Q6600 2.4GHz) 
   Operating System: 32-bit Windows Vista Business
   Speed of encrypting long message:
   3) 3.5 cycles/byte  compiler: Intel C++ compiler 11.1    compilation option: icl /O2
   4) 4.4 cycles/byte  compiler: Microsoft Visual C++ 2008  compilation option: release

   ------------------------------------------
   In this simplified optimization program, loop unrolling is applied to the description of HC-256 directly.
   16 steps are executed in each loop.

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

/*define data alignment for different C compilers*/
#if defined(__GNUC__)
      #define DATA_ALIGN16(x) x __attribute__ ((aligned(16)))
#else
      #define DATA_ALIGN16(x) __declspec(align(16)) x
#endif

typedef struct { 
      DATA_ALIGN16(uint32 P[1024]);
      DATA_ALIGN16(uint32 Q[1024]);
      DATA_ALIGN16(uint32 keystream[16]);   /*the output of the stream cipher: 16 32-bit keystream word*/
      uint32 counter2048;                   /*counter2048 = i mod 2048 */
} HC256_State; 

/*this function right rotates a 32-bit word x by n positions*/ 
#define ROTR32(x,n)  (((x)>>(n))|((x)<<(32-(n))))

/*the h1 function in HC-256*/
#define h1(x,y) {              \
      a = (uint8)(x);          \
      b = (uint8)((x) >> 8);   \
      c = (uint8)((x) >> 16);  \
      d = (uint8)((x) >> 24);  \
      (y) = state->Q[a]+state->Q[256+b]+state->Q[512+c]+state->Q[768+d]; \
}

/*the h2 function in HC-256*/
#define h2(x,y) {              \
      a = (uint8)(x);          \
      b = (uint8)((x) >> 8);   \
      c = (uint8)((x) >> 16);  \
      d = (uint8)((x) >> 24);  \
      (y) = state->P[a]+state->P[256+b]+state->P[512+c]+state->P[768+d]; \
}

/* one step in HC-256 when table P is updated:
   an element of table P is updated;
   a 32-bit keystream word is generated
*/ 
#define step_P(m0,m1023,m3,m10,m12,s) {         \
      tem0 = ROTR32((m1023),23);                \
      tem1 = ROTR32((m3),10);                   \
      tem2 = ((m3) ^ (m1023)) & 0x3ff;          \
      (m0) += (m10)+(tem0^tem1)+state->Q[tem2]; \
      h1((m12),tem3);                           \
      (s) = tem3 ^ (m0) ;                       \
}       

/* one step in HC-256 when table Q is updated:
   an element of table Q is updated;
   a 32-bit keystream word is generated
*/ 
#define step_Q(m0,m1023,m3,m10,m12,s) {         \
      tem0 = ROTR32((m1023),23);                \
      tem1 = ROTR32((m3),10);                   \
      tem2 = ((m3) ^ (m1023)) & 0x3ff;          \
      (m0) += (m10)+(tem0^tem1)+state->P[tem2]; \
      h2((m12),tem3);                           \
      (s) = tem3 ^ (m0) ;                       \
}       

/* this function computes 16 steps of HC-256
   the state is updated, 
   16 32-bit keystream words are generated and stored in the array state->keystream
*/
void SixteenSteps(HC256_State *state)  
{
      uint32 tem0,tem1,tem2,tem3;
      uint8 a,b,c,d;
      uint32 cc,dd,ee;

      cc = state->counter2048 & 0x3ff;
      dd = (cc+16)&0x3ff;
      ee = (cc-16)&0x3ff;

      if (state->counter2048 < 1024)	
      {   		
            step_P(state->P[cc+0], state->P[cc+1], state->P[ee+13],state->P[ee+6], state->P[ee+4], state->keystream[0]);
            step_P(state->P[cc+1], state->P[cc+2], state->P[ee+14],state->P[ee+7], state->P[ee+5], state->keystream[1]);
            step_P(state->P[cc+2], state->P[cc+3], state->P[ee+15],state->P[ee+8], state->P[ee+6], state->keystream[2]);
            step_P(state->P[cc+3], state->P[cc+4], state->P[cc+0], state->P[ee+9], state->P[ee+7], state->keystream[3]);
            step_P(state->P[cc+4], state->P[cc+5], state->P[cc+1], state->P[ee+10],state->P[ee+8], state->keystream[4]);
            step_P(state->P[cc+5], state->P[cc+6], state->P[cc+2], state->P[ee+11],state->P[ee+9], state->keystream[5]);
            step_P(state->P[cc+6], state->P[cc+7], state->P[cc+3], state->P[ee+12],state->P[ee+10],state->keystream[6]);
            step_P(state->P[cc+7], state->P[cc+8], state->P[cc+4], state->P[ee+13],state->P[ee+11],state->keystream[7]);
            step_P(state->P[cc+8], state->P[cc+9], state->P[cc+5], state->P[ee+14],state->P[ee+12],state->keystream[8]);
            step_P(state->P[cc+9], state->P[cc+10],state->P[cc+6], state->P[ee+15],state->P[ee+13],state->keystream[9]);
            step_P(state->P[cc+10],state->P[cc+11],state->P[cc+7], state->P[cc+0], state->P[ee+14],state->keystream[10]);
            step_P(state->P[cc+11],state->P[cc+12],state->P[cc+8], state->P[cc+1], state->P[ee+15],state->keystream[11]);
            step_P(state->P[cc+12],state->P[cc+13],state->P[cc+9], state->P[cc+2], state->P[cc+0], state->keystream[12]);
            step_P(state->P[cc+13],state->P[cc+14],state->P[cc+10],state->P[cc+3], state->P[cc+1], state->keystream[13]);
            step_P(state->P[cc+14],state->P[cc+15],state->P[cc+11],state->P[cc+4], state->P[cc+2], state->keystream[14]);
            step_P(state->P[cc+15],state->P[dd+0], state->P[cc+12],state->P[cc+5], state->P[cc+3], state->keystream[15]);
      }
      else				    
      {
            step_Q(state->Q[cc+0], state->Q[cc+1], state->Q[ee+13],state->Q[ee+6], state->Q[ee+4], state->keystream[0]);
            step_Q(state->Q[cc+1], state->Q[cc+2], state->Q[ee+14],state->Q[ee+7], state->Q[ee+5], state->keystream[1]);
            step_Q(state->Q[cc+2], state->Q[cc+3], state->Q[ee+15],state->Q[ee+8], state->Q[ee+6], state->keystream[2]);
            step_Q(state->Q[cc+3], state->Q[cc+4], state->Q[cc+0], state->Q[ee+9], state->Q[ee+7], state->keystream[3]);
            step_Q(state->Q[cc+4], state->Q[cc+5], state->Q[cc+1], state->Q[ee+10],state->Q[ee+8], state->keystream[4]);
            step_Q(state->Q[cc+5], state->Q[cc+6], state->Q[cc+2], state->Q[ee+11],state->Q[ee+9], state->keystream[5]);
            step_Q(state->Q[cc+6], state->Q[cc+7], state->Q[cc+3], state->Q[ee+12],state->Q[ee+10],state->keystream[6]);
            step_Q(state->Q[cc+7], state->Q[cc+8], state->Q[cc+4], state->Q[ee+13],state->Q[ee+11],state->keystream[7]);
            step_Q(state->Q[cc+8], state->Q[cc+9], state->Q[cc+5], state->Q[ee+14],state->Q[ee+12],state->keystream[8]);
            step_Q(state->Q[cc+9], state->Q[cc+10],state->Q[cc+6], state->Q[ee+15],state->Q[ee+13],state->keystream[9]);
            step_Q(state->Q[cc+10],state->Q[cc+11],state->Q[cc+7], state->Q[cc+0], state->Q[ee+14],state->keystream[10]);
            step_Q(state->Q[cc+11],state->Q[cc+12],state->Q[cc+8], state->Q[cc+1], state->Q[ee+15],state->keystream[11]);
            step_Q(state->Q[cc+12],state->Q[cc+13],state->Q[cc+9], state->Q[cc+2], state->Q[cc+0], state->keystream[12]);
            step_Q(state->Q[cc+13],state->Q[cc+14],state->Q[cc+10],state->Q[cc+3], state->Q[cc+1], state->keystream[13]);
            step_Q(state->Q[cc+14],state->Q[cc+15],state->Q[cc+11],state->Q[cc+4], state->Q[cc+2], state->keystream[14]);
            step_Q(state->Q[cc+15],state->Q[dd+0], state->Q[cc+12],state->Q[cc+5], state->Q[cc+3], state->keystream[15]);
      }
      state->counter2048 = (state->counter2048 + 16) & 0x7ff;
}

/*The following defines the initialization functions*/

/*the functions used for expanding the key and iv*/
#define f1(x)  (ROTR32((x),7) ^ ROTR32((x),18) ^ ((x) >> 3))
#define f2(x)  (ROTR32((x),17) ^ ROTR32((x),19) ^ ((x) >> 10))
#define f(a,b,c,d) (f2((a)) + b + f1((c)) + d) 

/*update one element in table P*/
#define update_P(m0,m1023,m3,m10) {              \
      tem0 = ROTR32((m1023),23);                 \
      tem1 = ROTR32((m3),10);                    \
      tem2 = ((m3) ^ (m1023)) & 0x3ff;           \
      (m0) += (m10)+(tem0^tem1)+state->Q[tem2];  \
}       

/*update one element in table P*/
#define update_Q(m0,m1023,m3,m10) {              \
      tem0 = ROTR32((m1023),23);                 \
      tem1 = ROTR32((m3),10);                    \
      tem2 = ((m3) ^ (m1023)) & 0x3ff;           \
      (m0) += (m10)+(tem0^tem1)+state->P[tem2];  \
}       

/*update the state for 16 steps, without generating keystream*/
void UpdateSixteenSteps(HC256_State *state)  
{
      uint32 tem0,tem1,tem2;
      uint32 cc,dd,ee;

      cc = state->counter2048 & 0x3ff;
      dd = (cc+16)&0x3ff;
      ee = (cc-16)&0x3ff;

      if (state->counter2048 < 1024)	
      {   		
            update_P(state->P[cc+0], state->P[cc+1], state->P[ee+13],state->P[ee+6]);
            update_P(state->P[cc+1], state->P[cc+2], state->P[ee+14],state->P[ee+7]);
            update_P(state->P[cc+2], state->P[cc+3], state->P[ee+15],state->P[ee+8]);
            update_P(state->P[cc+3], state->P[cc+4], state->P[cc+0], state->P[ee+9]);
            update_P(state->P[cc+4], state->P[cc+5], state->P[cc+1], state->P[ee+10]);
            update_P(state->P[cc+5], state->P[cc+6], state->P[cc+2], state->P[ee+11]);
            update_P(state->P[cc+6], state->P[cc+7], state->P[cc+3], state->P[ee+12]);
            update_P(state->P[cc+7], state->P[cc+8], state->P[cc+4], state->P[ee+13]);
            update_P(state->P[cc+8], state->P[cc+9], state->P[cc+5], state->P[ee+14]);
            update_P(state->P[cc+9], state->P[cc+10],state->P[cc+6], state->P[ee+15]);
            update_P(state->P[cc+10],state->P[cc+11],state->P[cc+7], state->P[cc+0]);
            update_P(state->P[cc+11],state->P[cc+12],state->P[cc+8], state->P[cc+1]);
            update_P(state->P[cc+12],state->P[cc+13],state->P[cc+9], state->P[cc+2]);
            update_P(state->P[cc+13],state->P[cc+14],state->P[cc+10],state->P[cc+3]);
            update_P(state->P[cc+14],state->P[cc+15],state->P[cc+11],state->P[cc+4]);
            update_P(state->P[cc+15],state->P[dd+0], state->P[cc+12],state->P[cc+5]);
      }
      else				    
      {
            update_Q(state->Q[cc+0], state->Q[cc+1], state->Q[ee+13],state->Q[ee+6]);
            update_Q(state->Q[cc+1], state->Q[cc+2], state->Q[ee+14],state->Q[ee+7]);
            update_Q(state->Q[cc+2], state->Q[cc+3], state->Q[ee+15],state->Q[ee+8]);
            update_Q(state->Q[cc+3], state->Q[cc+4], state->Q[cc+0], state->Q[ee+9]);
            update_Q(state->Q[cc+4], state->Q[cc+5], state->Q[cc+1], state->Q[ee+10]);
            update_Q(state->Q[cc+5], state->Q[cc+6], state->Q[cc+2], state->Q[ee+11]);
            update_Q(state->Q[cc+6], state->Q[cc+7], state->Q[cc+3], state->Q[ee+12]);
            update_Q(state->Q[cc+7], state->Q[cc+8], state->Q[cc+4], state->Q[ee+13]);
            update_Q(state->Q[cc+8], state->Q[cc+9], state->Q[cc+5], state->Q[ee+14]);
            update_Q(state->Q[cc+9], state->Q[cc+10],state->Q[cc+6], state->Q[ee+15]);
            update_Q(state->Q[cc+10],state->Q[cc+11],state->Q[cc+7], state->Q[cc+0]);
            update_Q(state->Q[cc+11],state->Q[cc+12],state->Q[cc+8], state->Q[cc+1]);
            update_Q(state->Q[cc+12],state->Q[cc+13],state->Q[cc+9], state->Q[cc+2]);
            update_Q(state->Q[cc+13],state->Q[cc+14],state->Q[cc+10],state->Q[cc+3]);
            update_Q(state->Q[cc+14],state->Q[cc+15],state->Q[cc+11],state->Q[cc+4]);
            update_Q(state->Q[cc+15],state->Q[dd+0], state->Q[cc+12],state->Q[cc+5]);
      }
      state->counter2048 = (state->counter2048 + 16) & 0x7ff;
}

/*initialization of the cipher, the key and iv are used to update the state*/
void Initialization(HC256_State *state, uint8 *key, uint8 *iv) 
{   
      uint32 i;

      /*expand the key and iv into P and Q*/
      for (i = 0; i < 8; i++)   state->P[i] = ((uint32*)key)[i];
      for (i = 8; i < 16; i++)  state->P[i] = ((uint32*)iv)[i-8]; 

      for (i = 16; i < 528; i++)  
            state->P[i] = f(state->P[i-2],state->P[i-7],state->P[i-15],state->P[i-16])+i;
      for (i = 0; i < 16; i++)    
            state->P[i] = state->P[i+512];
      for (i = 16; i < 1024; i++) 
            state->P[i] = f(state->P[i-2],state->P[i-7],state->P[i-15],state->P[i-16])+512+i;

      for (i = 0;  i < 16;  i++)  
            state->Q[i] = state->P[1024-16+i];
      for (i = 16; i < 32;  i++)  
            state->Q[i] = f(state->Q[i-2],state->Q[i-7],state->Q[i-15],state->Q[i-16])+1520+i;
      for (i = 0;  i < 16;  i++)  
            state->Q[i] = state->Q[i+16];
      for (i = 16; i < 1024;i++)  
            state->Q[i] = f(state->Q[i-2],state->Q[i-7],state->Q[i-15],state->Q[i-16])+1536+i;

      /*initialize counter2048*/ 
      state->counter2048 = 0;

      /*run the cipher 4096 steps without generating keystream*/	
      for (i = 0; i < 256; i++) UpdateSixteenSteps(state);
}	

/* encrypt a message
   three inputs to this function: cipher state, message, the message length in bytes
   one output:  ciphertext
*/
void EncryptMessage(HC256_State *state, uint8 *message, uint8 *ciphertext, uint64 msglength)
{
      uint64 i;
      uint32 j;

      /*encrypt a message, each time 64 bytes are encrypted*/
      for (i = 0; (i+64) <= msglength; i += 64, message += 64, ciphertext += 64) {
            /*generate 16 32-bit keystream and store it in state.keystream*/
            SixteenSteps(state);   
            /*encrypt 64 bytes of the message*/
            for (j = 0; j < 16; j++) ((uint32*)ciphertext)[j] = ((uint32*)message)[j] ^ state->keystream[j];   
      }

      /*encrypt the last message block if the message length is not multiple of 64 bytes*/
      if ((msglength & 0x3f) != 0) {
            SixteenSteps(state);
            for (j = 0; j < (msglength & 0x3f); j++) {
                  *(ciphertext+j) = *(message+j) ^ *( ((uint8*)state->keystream) +j);
            }
      }    
}

/* this function encrypts a message,
   there are four inputs to this function: a 256-bit key, a 256-bit iv, a message, the message length in bytes
   one output from this function: ciphertext
*/
void HC256(uint8 *key, uint8 *iv, uint8 *message,uint8 *ciphertext, uint64 msglength)
{
      HC256_State state;

      /*initializing the state*/
      Initialization(&state,key,iv);

      /*encrypt a message*/
      EncryptMessage(&state,message,ciphertext,msglength);
}
