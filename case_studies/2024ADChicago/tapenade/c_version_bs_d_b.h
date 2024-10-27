#ifndef C_VERSION_BS_D_B_INCLUDED
#define C_VERSION_BS_D_B_INCLUDED

void call_price_d_b(double S, double *Sb, double Sd, double *Sdb, double K,
                    double *Kb, double Kd, double *Kdb, double v, double *vb,
                    double vd, double *vdb, double T, double *Tb, double Td,
                    double *Tdb, double *call_price, double *call_priceb,
                    double call_price_db);

#endif
