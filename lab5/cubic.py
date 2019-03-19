import numpy as np
import matplotlib.pyplot as plt
from math import pi,sqrt

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

def spline(x,y,n,bc):
	h=[]
	b=[]
	for i in range(n):
		h.append(x[i+1]-x[i])
		b.append((y[i+1]-y[i])/h[i])

	v=[]
	u=[]
	for i in range(1,n):
		v.append(2*(h[i-1]+h[i]))
		u.append(6*(b[i]-b[i-1]))

	a=h.copy()
	c=h.copy()

	if bc=="natural":
		a[0]=0
		c[len(c)-1]=0
		del c[0]
		del a[len(a)-1]

		z=thomas(a,v,c,u,n-1)
		z.append(0)
		z.insert(0,0)
	else:
		a.insert(0,0)
		c.append(0)
		a[len(a)-1]=1
		v.append(2)
		v.insert(0,2)
		c[0]=1

		last=((-6)/(x[len(x)-1]-x[len(x)-2])*((y[len(y)-1]-y[len(y)-2])/(x[len(x)-1]-x[len(x)-2])))
		first=(6/(x[1]-x[0])*((y[1]-y[0])/(x[1]-x[0])))
		u.append(last)
		u.insert(0,first)
		z=thomas(a,v,c,u,n+1)

	return(x,y,z,h)

def s(x,y,z,h,arg,i):
	w=(z[i+1]/(6*h[i]))*((arg-x[i])**3)+(z[i]/(6*h[i]))*((x[i+1]-arg)**3)
	w+=(arg-x[i])*((y[i+1]/h[i])-((z[i+1]/6)*h[i]))
	w+=(x[i+1]-arg)*((y[i]/h[i])-((h[i]/6)*z[i]))
	return w

def norm(f,y,n):
	s=0
	for i in range(n):
		s+=(f[i]-y[i])**2
	s=sqrt(s)
	return s

print("Choose number of nodes:")
pkt=int(input())
(a,b)=(0.1,0.8)

x=np.linspace(a,b,num=pkt,endpoint=True)
y=x*np.sin((0.8*pi)/x)

print("Choose boundary condition:")
print("1) Natural (second derivatives at the edges of the interval are equal to zero")
print("2) Clamped (first derivatives at the edges of the interval are equal to zero)")
task=int(input())

if task==1: (x,y,z,h)=spline(x,y,pkt-1,"natural")
else: (x,y,z,h)=spline(x,y,pkt-1,"clamped")

ile=1000/(pkt-1)
ile=int(round(ile))
xx=[]
yy=[]
org=[]
for i in range(pkt-1):
	for j in range(ile):
		arg=x[i]+j*(x[i+1]-x[i])/ile
		xx.append(arg)
		yy.append(s(x,y,z,h,arg,i))
		org.append(arg*np.sin((0.8*pi)/arg))
	arg=x[i+1]
	xx.append(arg)
	yy.append(s(x,y,z,h,arg,i))
	org.append(arg*np.sin((0.8*pi)/arg))

s=norm(org,yy,ile*(pkt-1))
print("Norm:")
print(s)

plt.scatter(xx,org,s=2)
plt.scatter(xx,yy,s=2,c='r')
plt.scatter(x,y,s=10,c='black')
plt.show()