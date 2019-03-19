import numpy as np
import matplotlib.pyplot as plt
from math import pi,sqrt

def spline(x,y,n,bc):
	z=[0]*n
	if bc!="natural":
		z[0]=(y[1]-y[0])/(x[1]-x[0])
	for i in range(0,n-1):
		z[i+1]=-z[i]+2*((y[i+1]-y[i])/(x[i+1]-x[i]))
	return(x,y,z)

def s(x,y,z,arg,i):
	return y[i]+z[i]*(arg-x[i])+((z[i+1]-z[i])/(2*(x[i+1]-x[i])))*((arg-x[i])**2)

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
print("1) First derivative in the boundary point is equal to zero")
print("2) First spline is linear")
task=int(input())
if task==1: bc="natural"
else: bc=""
(x,y,z)=spline(x,y,pkt,bc)

ile=1000/(pkt-1)
ile=int(round(ile))
xx=[]
yy=[]
org=[]
for i in range(pkt-1):
	for j in range(ile):
		arg=x[i]+j*(x[i+1]-x[i])/ile
		xx.append(arg)
		yy.append(s(x,y,z,arg,i))
		org.append(arg*np.sin((0.8*pi)/arg))
	arg=x[i+1]
	xx.append(arg)
	yy.append(s(x,y,z,arg,i))
	org.append(arg*np.sin((0.8*pi)/arg))

s=norm(org,yy,ile*(pkt-1))
print("Norm:")
print(s)

plt.scatter(xx,org,s=2)
plt.scatter(xx,yy,s=2,c='r')
plt.scatter(x,y,s=10,c='black')
plt.show()