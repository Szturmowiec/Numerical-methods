from numpy import *
from math import sqrt

def matrixx(n):
	A=[[float(0)]*n for i in range(n)]
	for i in range(1,n+1):
		for j in range(1,n+1):
			if i==j: A[i-1][j-1]=float(10)
			else: A[i-1][j-1]=float(1)/((abs(float(i)-float(j)))+float(5))
	return A

def euklides2(x,n):
	s=float(0)
	for i in range(0,n): s=s+x[i]*x[i]
	s=sqrt(s)
	return s

def calc_b(A,x,n):
	b=[0]*n
	for i in range(n):
		for k in range(n):
			b[i]+=A[i][k]*x[k]
	return b

def sor(A,b,w,n,p,task):
	D=[[float(0)]*n for i in range(n)]
	Lw=[[float(0)]*n for i in range(n)]
	Uw=[[float(0)]*n for i in range(n)]
	for i in range(0,n):
		for j in range(0,n):
			if i==j: D[i][j]=A[i][j]
			elif i>j: Lw[i][j]=A[i][j]*w
			else: Uw[i][j]=A[i][j]*w
	x=[float(0)]*n

	DL=[[float(0)]*n for i in range(n)]
	DU=[[float(0)]*n for i in range(n)]
	DLw=[[0]*n for i in range(n)]
	for i in range(0,n):
		for j in range(0,n):
			DL[i][j]=D[i][j]+Lw[i][j]
			DU[i][j]=(float(1)-w)*D[i][j]-Uw[i][j]

	DL=asarray(DL)
	DU=asarray(DU)
	DLw=asarray(DLw)
	b=asarray(b)
	DL=linalg.inv(DL)
	DLw=DL.dot(w)

	S1=DLw.dot(b)
	S2=DL.dot(DU)
	it=0
	stop=float(0)

	while(stop>=p or it==0):
		x2=S2.dot(x)+S1
		if task==1:
			dif=x2-x
			stop=euklides2(dif,n)
		else:
			bi=calc_b(A,x,n)-b
			stop=euklides2(bi,n)

		x=x2
		it+=1
	return (x,it)

def euklides(xd,xw,n):
	s=float(0)
	for i in range(0,n): s=s+(xd[i]-xw[i])*(xd[i]-xw[i])
	s=sqrt(s)
	return s

print("Matrix size:")
n=int(input())
print("Iterations limiter:")
p=float(input())

A=matrixx(n)
x=[float(-2)]*n
#for i in range(0,n):
#	if i%2==0: x[i]=float(1)
#	else: x[i]=float(-1)

b=calc_b(A,x,n)

print("Choose iterations limiter method:")
print("1) ||x_(i+1)-x_(i)||<p")
print("2) ||Ax_(i)-b||<p")
task=int(input())
print("Choose w parameter:")
w=float(input())
(xw,it)=sor(A,b,w,n,p,task)
e=euklides(x,xw,n)

print("Euklide's norm of difference between x given and calculated:")
print(e)
print("Number of iterations:")
print(it)
