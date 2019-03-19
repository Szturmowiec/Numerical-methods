import matplotlib.pyplot as plt
import numpy as np
from math import pi,sin,cos

def thomas(a,b,c,d,n):
	beta=[(0)]*n
	gamma=[(0)]*n
	beta[0]=((-1)*c[0])/b[0]
	gamma[0]=d[0]/b[0]

	for i in range(1,n):
		beta[i]=((-1)*c[i])/(a[i]*beta[i-1]+b[i])
		gamma[i]=(d[i]-a[i]*gamma[i-1])/(a[i]*beta[i-1]+b[i])

	xt=[(0)]*n
	xt[n-1]=gamma[n-1]
	for i in range(n-2,-1,-1):
		xt[i]=beta[i]*xt[i+1]+gamma[i]
	return xt

def f(x):
    return 8*x

def w(x):
    return -2*sin(2*x)+2*x

def norm(org,y,n):
	s=0
	for i in range(n):
		b=abs(y[i]-org[i])
		if b>s: s=b
	return s

def rozw(aa,bb,n,h):
    a=[-1]*n
    b=[]
    c=[-1]*n
    d=[]
    x=aa
    arg=[x]
    q=-4

    a[0]=0
    c[0]=0
    d.append(w(aa))
    b.append(1)

    for i in range(1,n-1):
        x=x+h
        arg.append(x)
        b.append(2+h*h*q)
        d.append(-h*h*f(x))

    arg.append(bb)
    b.append(1)
    d.append(w(bb))
    c[n-1]=0
    a[n-1]=0
    wyn=thomas(a,b,c,d,n)
    return (arg,wyn)

(a,b)=(0,(2*pi+2)/2)
print("N: (h=(a-b)/(N-1))")
n=int(input())
h=(b-a)/(n-1)
(x,y)=rozw(a,b,n,h)

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
