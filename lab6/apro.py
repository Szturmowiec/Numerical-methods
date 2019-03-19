import numpy as np
import matplotlib.pyplot as plt
from math import pi,sqrt,cos

def solv(A,b,n):
	for i in range(n):
		A[i].append(b[i])

	for i in range(n-1):
		maks=abs(A[i][i])
		maxid=i
		for j in range(i+1,n):
			if abs(A[j][i])>maks:
				maks=abs(A[j][i])
				maxid=j

		for k in range(i,n+1):
			tmp=A[maxid][k]
			A[maxid][k]=A[i][k]
			A[i][k]=tmp

		d=A[i][i]
		for k in range(i,n+1):
			A[i][k]=A[i][k]/d

		for z in range(i+1,n):
			g=A[z][i]
			for k in range(i,n+1):
				A[z][k]-=g*A[i][k]

	x=[0]*n
	for i in range(n-1,-1,-1):
		for j in range(n-1,i,-1):
			A[i][n]-=x[j]*A[i][j]
		x[i]=A[i][n]/A[i][i]
	return x

def czeb(a,b,n):
	x=[]
	for k in range(n):
		root=cos(pi*((2*(k+1)-1)/(2*n)))
		w=((b-a)*root+a+b)/2
		x.append(w)
	return x

def norm(f,y,n):
	s=0
	for i in range(n):
		s+=(f[i]-y[i])**2
	s=sqrt(s/n)
	return s

def apro(x,y,n,m):
	A=[[0]*(m+1) for i in range(m+1)]
	for k in range(m+1):
		for j in range(m+1):
			s=0
			for i in range(n+1):
				s+=x[i]**(j+k)
			A[k][j]=s

	b=[0]*(m+1)
	for k in range(m+1):
		s=0
		for i in range(n+1):
			s+=y[i]*x[i]**k
		b[k]=s

	return solv(A,b,m+1)

def f(c,x):
	w=0
	for i in range(len(c)):
		w+=c[i]*x**(i)
	return w

print("Number of nodes:")
pkt=int(input())
(a,b)=(0.1,0.8)

x=np.linspace(a,b,num=pkt,endpoint=True)
x2=czeb(a,b,pkt)
x2.sort(key=float)
y=x*np.sin((0.8*pi)/x)
y2=[]
for i in range(pkt):
	y2.append(x2[i]*np.sin((0.8*pi)/x2[i]))

print("Polynomial degree:")
m=int(input())

print("Choose nodes:")
print("1) Equal")
print("2) Chebyshev")
task=int(input())

if task==1: c=apro(x,y,pkt-1,m)
else: c=apro(x2,y2,pkt-1,m)

ile=1000
xx=np.linspace(a,b,num=ile,endpoint=True)
yy=[]
org=[]

for i in range(ile):
	yy.append(f(c,xx[i]))
	org.append(xx[i]*np.sin((0.8*pi)/xx[i]))

s=norm(org,yy,ile)
print("Norm:")
print(s)

plt.scatter(xx,org,s=2)
plt.scatter(xx,yy,s=2,c='r')
if task==1: plt.scatter(x,y,s=10,c='black')
else: plt.scatter(x2,y2,s=10,c='black')
plt.ylim((-0.6,0.4))
plt.show()
