#define m1 2147483647
#define m2 2147462579
#define norm 4.656612873077393e-10
#define mask12 511      // 111111111 (9)
#define mask13 16777215 // 111111111111111111111111 (24)
#define mask21 65535    // 1111111111111111 (16)
#define SEED 12345

unsigned long x10 = SEED, x11 = SEED, x12 = SEED, x20 = SEED, x21 = SEED,
              x22 = SEED;

double MRG31k3p() {
    register unsigned long y1, y2; /* For intermediate results */
    /* First component */
    y1 = (((x11 & mask12) << 22) + (x11 >> 9)) +
         (((x12 & mask13) << 7) + (x12 >> 24));
    if (y1 > m1) {
        y1 -= m1;
    }

    y1 += x12;

    if (y1 > m1) {
        y1 -= m1;
    }
    x12 = x11;
    x11 = x10;
    x10 = y1;

    /* Second component */
    y1 = ((x20 & mask21) << 15) + 21069 * (x20 >> 16);
    if (y1 > m2) {
        y1 -= m2;
    }
    y2 = ((x22 & mask21) << 15) + 21069 * (x22 >> 16);
    if (y2 > m2) {
        y2 -= m2;
    }
    y2 += x22;
    if (y2 > m2) {
        y2 -= m2;
    }
    y2 += y1;
    if (y2 > m2) {
        y2 -= m2;
    }
    x22 = x21;
    x21 = x20;
    x20 = y2;

    /* Combinaison */
    if (x10 <= x20) {
        return ((x10 - x20 + m1) * norm);
    } else {
        return ((x10 - x20) * norm);
    }
}

// Juan Lucas Rey: add initialiser
void InitMRG31k3p (unsigned long *init)
{
    x12 = init[0];
    x11 = init[1];
    x10 = init[2];
    x22 = init[3];
    x21 = init[4];
    x20 = init[5];
}
