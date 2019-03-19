from decimal import *
from numpy import *
from scipy import linalg
from sympy import *

def matrixx(n):
	A=[[Decimal(0)]*n for i in range(n)]
	for i in range(1,n+1):
		for j in range(1,n+1):
			if i==j: A[i-1][j-1]=Decimal(10)
			else: A[i-1][j-1]=Decimal(1)/((abs(Decimal(i)-Decimal(j)))+Decimal(5))
	return A

def euklides2(x,n):
	s=Decimal(0)
	for i in range(0,n): s=s+x[i]*x[i]
	s=s.sqrt()
	return s

def calc_b(A,x,n):
	b=[Decimal(0)]*n
	for i in range(n):
		for k in range(n):
			b[i]+=A[i][k]*x[k]
	return b

def spekt(M,n):
	#l=symbols('l')
	#I=eye(n)
	#A=Matrix([[0]*n for i in range(n)])
	#k=0
	#for i in range(0,n):
	#	for j in range(0,n):
	#		A[k]=M[i][j]
	#		k+=1

	#e=Eq((A-l*I).det(),0)
	#print(e)
	#S=solve(e)
	#s=abs(S[0])
	for i in range(0,n):
		for j in range(0,n):
			M[i][j]=float(M[i][j])
	S=linalg.eigvals(M)

	s=abs(S[0])
	for i in range(1,len(S)):
		if abs(S[i])>s: s=abs(S[i])
	return s

def jacob(A,b,n,p,task):
	N=[Decimal(0)]*n
	LU=[[Decimal(0)]*n for i in range(n)]
	M=[[Decimal(0)]*n for i in range(n)]
	for i in range(0,n):
		for j in range(0,n):
			if i==j:
				N[i]=Decimal(1)/A[i][j]
			else: LU[i][j]=A[i][j]

	for j in range(0,n):
		for i in range(0,n):
			M[j][i]-=LU[j][i]*N[j]

	Nb=[Decimal(0)]*n
	for i in range(0,n):
		Nb[i]=N[i]*b[i]

	stop=Decimal(0)
	x=[Decimal(0)]*n
	dif=[Decimal(0)]*n
	x2=[Decimal(0)]*n
	it=0

	while(stop>=p or it==0):
		print(stop)
		for i in range(0,n):
			x2[i]=Nb[i]
			for j in range(0,n):
				x2[i]+=M[i][j]*x[j]

		if task==1:
			for i in range(0,n):
				dif[i]=x2[i]-x[i]
			stop=euklides2(dif,n)
		else:
			bi=calc_b(A,x,n)
			for i in range(0,n):
				bi[i]-=b[i]
			stop=euklides2(bi,n)

		for i in range(0,n): x[i]=x2[i]
		it+=1
	return (x,it,M)

def euklides(xd,xw,n):
	s=Decimal(0)
	for i in range(0,n): s=s+(xd[i]-xw[i])*(xd[i]-xw[i])
	s=s.sqrt()
	return s

print("Matrix size:")
n=int(input())
print("Iterations limiter:")
pp=Decimal(input())
print("Choose precision:")
p=int(input())
getcontext().prec=p

A=matrixx(n)
x=[Decimal(2)]*n
#for i in range(0,n):
#	if i%2==0: x[i]=Decimal(1)
#	else: x[i]=Decimal(-1)

b=calc_b(A,x,n)

print("Choose iterations limiter method:")
print("1) ||x_(i+1)-x_(i)||<p")
print("2) ||Ax_(i)-b||<p")
task=int(input())

xw=[Decimal(0)]*n
(xw,it,M)=jacob(A,b,n,pp,task)
e=euklides(x,xw,n)
s=spekt(M,n)

print("Euklide's norm of difference between x given and calculated:")
print(e)
print("Number of iterations:")
print(it)
print("Spectral radius:")
print(s)
