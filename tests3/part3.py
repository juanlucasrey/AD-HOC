from sympy import diff, log, sin, cos 
from sympy.abc import x,y
expr=sin(y * x)

# derx4 = diff(expr,x,4)
# print(derx4)

# derx3y1 = diff(expr,x,3, y, 1)
# print(derx3y1)

# derx2y2 = diff(expr,x,2, y, 2)
# print(derx2y2)

# derx1y3 = diff(expr,x,1, y, 3)
# print(derx1y3)

# dery4 = diff(expr,y,4)
# print(dery4)

dery4 = diff(expr,y,4)
print(dery4)

dery3x1 = diff(expr, y, 3, x, 1)
print(dery3x1)

dery2x2 = diff(expr, y, 2, x, 2)
print(dery2x2)

dery1x3 = diff(expr, y, 1, x, 3)
print(dery1x3)

derx4 = diff(expr, x, 4)
print(derx4)

valy4 = dery4.subs(x, 0.3).subs(y, 0.5)
print(valy4)

valy3x1 = dery3x1.subs(x, 0.3).subs(y, 0.5)
print(valy3x1)

valy2x2 = dery2x2.subs(x, 0.3).subs(y, 0.5)
print(valy2x2)

valy1x3 = dery1x3.subs(x, 0.3).subs(y, 0.5)
print(valy1x3)

valx4 = derx4.subs(x, 0.3).subs(y, 0.5)
print(valx4)


sec1 = -y*y*3*cos(x*y)
print(sec1.subs(x, 0.3).subs(y, 0.5))
print(4*sec1.subs(x, 0.3).subs(y, 0.5))

sec2 = y*y*x*y*sin(x*y)
print(sec2.subs(x, 0.3).subs(y, 0.5))
print(4*sec2.subs(x, 0.3).subs(y, 0.5))

# sec3 = -6*cos(x*y)
# print(sec3.subs(x, 0.3).subs(y, 0.5))
# print(4*sec3.subs(x, 0.3).subs(y, 0.5))
# valx2y1 = derx2y1.subs(x, 0.3).subs(y, 0.5)
# print(valx2y1)

sec4 = y*y*x*y
print(sec4.subs(x, 0.3).subs(y, 0.5))
print(4*sec4.subs(x, 0.3).subs(y, 0.5))

# valx1y2 = derx1y2.subs(x, 0.3).subs(y, 0.5)
# print(valx1y2)

# valy3 = dery3.subs(x, 0.3).subs(y, 0.5)
# print(valy3)


# print("x1y2")
# f1 = diff(expr,x,1)
# print(f1)

# f2 = diff(f1,y,1)
# print(f2)

# f3 = diff(f2,y,1)
# print(f3)