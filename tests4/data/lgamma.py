# from mpmath import *
# mp.dps = 40
import mpmath
mpmath.mp.dps = 40


f = open('lgamma.csv', 'w+')
f.truncate()
start = -2.5
end = -2.4
increment = 0.01

val = start
while val < end:
    if val != 1.:
        print(val, end = ",", file=f)
        print(mpmath.loggamma(val), end = ",", file=f)
        print(mpmath.psi(0, val), end = ",", file=f)
        print(mpmath.psi(1, val), end = ",", file=f)
        print(mpmath.psi(2, val), end = ",", file=f)
        print(mpmath.psi(3, val), end = ",", file=f)
        print(mpmath.psi(4, val), end = ",", file=f)
        print(mpmath.psi(5, val), file=f)
    val = round(val + increment, 4)
