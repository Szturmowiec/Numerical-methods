import numpy as np
import matplotlib.pyplot as plt
from math import pi,sqrt,cos,sin

def norm(f,y,n):
	s=0
	for i in range(n):
		s+=(f[i]-y[i])**2
	s=sqrt(s/n)
	return s

def ff(x):
	return x*sin((0.8*pi)/x)

def apro(x,y,m,n):
	v=x.copy()
	for i in range(n+1): v[i]=(v[i]-0.1)*2*pi/0.7

	a=[]
	b=[]
	s=0
	for i in range(n): s+=y[i]
	a.append(s*2/(n+1))

	for k in range(1,m+1):
		s=0
		s2=0
		for i in range(n+1):
			s+=y[i]*cos(k*v[i])
			s2+=y[i]*sin(k*v[i])
		a.append((2/(n+1))*s)
		b.append((2/(n+1))*s2)
	return (a,b)

def f(a,b,x):
	w=0
	for k in range(len(b)):
		w+=a[k+1]*cos((k+1)*((x-0.1)*2*pi/0.7))+b[k]*sin((k+1)*((x-0.1)*2*pi/0.7))
	return w+a[0]/2

print("Number of nodes:")
pkt=int(input())
(a,b)=(0.1,0.8)

x=np.linspace(a,b,num=pkt,endpoint=True)
y=x*np.sin((0.8*pi)/x)

print("Polynomial degree:")
m=int(input())

(c,d)=apro(x,y,m,pkt-1)

ile=1000
xx=np.linspace(a,b,num=ile,endpoint=True)
yy=[]
org=[]

for i in range(ile):
	yy.append(f(c,d,xx[i]))
	org.append(ff(xx[i]))

s=norm(org,yy,ile)
print("Norm:")
print(s)

plt.scatter(xx,org,s=2)
plt.scatter(xx,yy,s=2,c='r')
plt.scatter(x,y,s=10,c='black')
plt.ylim((-0.6,0.4))
plt.show()