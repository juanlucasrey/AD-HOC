/*
 * & is bitwise AND, ^ is bitwise XOR, a<<b shifts a by b
 * ind(mm,x) is bits 2..9 of x, or (floor(x/4) mod 256)*4
 * in rngstep barrel(a) was replaced with a^(a<<13) or such
 */
typedef  unsigned int  u4;      /* unsigned four bytes, 32 bits */
typedef  unsigned char u1;      /* unsigned one  byte,  8  bits */
#define ind(mm,x)  (*(u4 *)((u1 *)(mm) + ((x) & (255<<2))))
#define rngstep(mix,a,b,mm,m,m2,r,x) \
{ \
  x = *m;  \
  a = (a^(mix)) + *(m2++); \
  *(m++) = y = ind(mm,x) + a + b; \
  *(r++) = b = ind(mm,y>>8) + x; \
}

/* static */ void isaac(u4 *mm,u4 *rr,u4 *aa,u4 *bb,u4 *cc)
// u4 *mm;      /* Memory: array of SIZE ALPHA-bit terms */
// u4 *rr;      /* Results: the sequence, same size as m */
// u4 *aa;      /* Accumulator: a single value */
// u4 *bb;      /* the previous result */
// u4 *cc;      /* Counter: one ALPHA-bit value */
{
  register u4 a,b,x,y,*m,*m2,*r,*mend;
  m=mm; r=rr;
  a = *aa; b = *bb + (++*cc);
  for (m = mm, mend = m2 = m+128; m<mend; )
  {
    rngstep( a<<13, a, b, mm, m, m2, r, x);
    rngstep( a>>6 , a, b, mm, m, m2, r, x);
    rngstep( a<<2 , a, b, mm, m, m2, r, x);
    rngstep( a>>16, a, b, mm, m, m2, r, x);
  }
  for (m2 = mm; m2<mend; )
  {
    rngstep( a<<13, a, b, mm, m, m2, r, x);
    rngstep( a>>6 , a, b, mm, m, m2, r, x);
    rngstep( a<<2 , a, b, mm, m, m2, r, x);
    rngstep( a>>16, a, b, mm, m, m2, r, x);
  }
  *bb = b; *aa = a;
}
