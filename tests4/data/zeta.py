import mpmath
mpmath.mp.dps = 40

f = open('zeta.csv', 'w+')
f.truncate()
start = -110
end = -100
increment = 0.01

val = start
while val < end:
    if val != 1.:
        print(val, end = ",", file=f)
        print(mpmath.zeta(val), end = ",", file=f)
        print(mpmath.zeta(val, derivative=1), end = ",", file=f)
        print(mpmath.zeta(val, derivative=2), end = ",", file=f)
        print(mpmath.zeta(val, derivative=3), end = ",", file=f)
        print(mpmath.zeta(val, derivative=4), end = ",", file=f)
        print(mpmath.zeta(val, derivative=5), end = ",", file=f)
        print(mpmath.zeta(val, derivative=6), file=f)
    val = round(val + increment, 4)
