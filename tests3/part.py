from sympy import diff, log, sin, cos 
from sympy.abc import x,y
expr=sin(log(y) * cos(x))

derx3 = diff(expr,x,3)
print(derx3)

derx2y1 = diff(expr,x,2, y, 1)
print(derx2y1)

derx1y2 = diff(expr,x,1, y, 2)
print(derx1y2)

dery3 = diff(expr,y,3)
print(dery3)

valx3 = derx3.subs(x, 0.3).subs(y, 0.5)
print(valx3)

valx2y1 = derx2y1.subs(x, 0.3).subs(y, 0.5)
print(valx2y1)

valx1y2 = derx1y2.subs(x, 0.3).subs(y, 0.5)
print(valx1y2)

valy3 = dery3.subs(x, 0.3).subs(y, 0.5)
print(valy3)


print("x1y2")
f1 = diff(expr,x,1)
print(f1)

f2 = diff(f1,y,1)
print(f2)

f3 = diff(f2,y,1)
print(f3)