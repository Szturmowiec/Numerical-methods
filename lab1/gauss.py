from decimal import *
import numpy as np

def matrix1(n):
	A=[[Decimal(0)]*n for i in range(n)]
	for i in range(1,n+1):
		for j in range(1,n+1):
			if i==1: A[i-1][j-1]=Decimal(1)
			else: A[i-1][j-1]=Decimal(1)/(Decimal(i)+Decimal(j)-Decimal(1))
	return A

def matrix2(n):
	A=[[Decimal(0)]*n for i in range(n)]
	for i in range(1,n+1):
		for j in range(1,n+1):
			if j>=i: A[i-1][j-1]=Decimal(2)*Decimal(i)/Decimal(j)
			else: A[i-1][j-1]=Decimal(A[j-1][i-1])
	return A

def euklides(xd,xw,n):
	s=Decimal(0)
	for i in range(0,n): s=s+(xd[i]-xw[i])*(xd[i]-xw[i])
	s=s.sqrt()
	return s

def war(A,n):
	B=[[0]*n for i in range(n)]
	for i in range(0,n):
		for j in range(0,n):
			B[i][j]=np.float64(A[i][j])

	I=np.linalg.inv(B)
	nA=0
	nI=0
	nAtmp=nA
	nItmp=nI
	for i in range(0,n):
		for j in range(0,n):
			nAtmp+=abs(B[i][j])
			nItmp+=abs(I[i][j])
		if nAtmp>nA: nA=nAtmp
		if nItmp>nI: nI=nItmp
		nAtmp=0
		nItmp=0
	return nA*nI

def gauss(A,b,n):
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

	x=[Decimal(0)]*n
	for i in range(n-1,-1,-1):
		for j in range(n-1,i,-1):
			A[i][n]-=x[j]*A[i][j]
		x[i]=A[i][n]/A[i][i]
	return x

def calc_b(A,x,n):
	b=[Decimal(0)]*n
	for i in range(n):
		for k in range(n):
			b[i]+=A[i][k]*x[k]
	return b

print("matrix size:")
n=int(input())
print("Choose precision:")
p=int(input())
getcontext().prec=p

print("Choose task:")
print("1)\n2)")
t=int(input())
if (t==1): A=matrix1(n)
else: A=matrix2(n)

x=[Decimal(0)]*n
for i in range(0,n):
	if i%2==0: x[i]=Decimal(1)
	else: x[i]=Decimal(-1)

b=calc_b(A,x,n)
w=war(A,n)

x2=[Decimal(0)]*n
x2=gauss(A,b,n)
e=euklides(x,x2,n)

print("Euklide's norm of difference between x given and calculated:")
print(e)
print("Condition number:")
print(w)