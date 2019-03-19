import matplotlib.pyplot as plt
import numpy as np
from math import pi,sin,cos,exp

def f(x,y):
    return 32*sin(2*x)*cos(2*x)+8*y*sin(2*x)

def w(x):
    return exp(-4*cos(2*x))-4*cos(2*x)+1

def norm(org,y,n):
	s=0
	for i in range(n):
		b=abs(y[i]-org[i])
		if b>s: s=b
	return s

def rozw(x,y,h,b):
    wyn=[y]
    arg=[x]
    while(x<b):
        y=y+h*f(x,y)
        x=x+h
        wyn.append(y)
        arg.append(x)
    return (arg,wyn)

(a,b)=(-pi/4,3*pi/2)
print("N: (h=(a-b)/(N-1))")
n=int(input())
h=(b-a)/(n-1)
print(h)
(x,y)=rozw(a,w(a),h,b)

xe=np.linspace(a,b,num=n,endpoint=True)
ye=[]
for i in range(n): ye.append(w(xe[i]))
err=norm(ye,y,n)
print(err)

org=[]
ile=300
xx=np.linspace(a,b,num=ile,endpoint=True)
for i in range(ile): org.append(w(xx[i]))
plt.scatter(xx,org,s=2)
plt.scatter(x,y,s=3,c='r')
plt.show()
