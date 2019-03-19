from decimal import *
from time import *

def matrix(n):
	a=[Decimal(0)]*n
	b=[Decimal(0)]*n
	c=[Decimal(0)]*n
	for i in range(1,n+1):
		for j in range(1,n+1):
			if i==j: b[i-1]=Decimal(-4)*Decimal(i)-Decimal(5)
			if j==i+1: c[i-1]=Decimal(i)
			if j==i-1 and i>1: a[i-1]=Decimal(4)/Decimal(i)
	return (a,b,c)

def euklides(xd,xw,n):
	s=Decimal(0)
	for i in range(0,n): s=s+(xd[i]-xw[i])*(xd[i]-xw[i])
	s=s.sqrt()
	return s

def thomas(a,b,c,d,n):
	beta=[Decimal(0)]*n
	gamma=[Decimal(0)]*n
	beta[0]=(Decimal(-1)*c[0])/b[0]
	gamma[0]=d[0]/b[0]

	for i in range(1,n):
		beta[i]=(Decimal(-1)*c[i])/(a[i]*beta[i-1]+b[i])
		gamma[i]=(d[i]-a[i]*gamma[i-1])/(a[i]*beta[i-1]+b[i])

	xt=[Decimal(0)]*n
	xt[n-1]=gamma[n-1]
	for i in range(n-2,-1,-1):
		xt[i]=beta[i]*xt[i+1]+gamma[i]
	return xt

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

	xg=[0]*n
	for i in range(n-1,-1,-1):
		for j in range(n-1,i,-1):
			A[i][n]-=xg[j]*A[i][j]
		xg[i]=A[i][n]/A[i][i]
	return xg

def calc_d(a,b,c,x,n):
	d=[Decimal(0)]*n
	d[0]=b[0]*x[0]+c[0]*x[1]
	for i in range(1,n-1):
		d[i]=a[i]*x[i-1]+b[i]*x[i]+c[i]*x[i+1]
	d[n-1]=a[n-1]*x[n-2]+b[n-1]*x[n-1]
	return d

def fullmatrix(a,b,c,n):
	A=[[Decimal(0)]*n for i in range(n)]
	for i in range(1,n+1):
		for j in range(1,n+1):
			if i==j: A[i-1][j-1]=b[i-1]
			if j==i+1: A[i-1][j-1]=c[i-1]
			if j==i-1 and i>1: A[i-1][j-1]=a[i-1]
	return A

print("matrix size:")
n=int(input())
print("Choose precision:")
p=int(input())
getcontext().prec=p

(a,b,c)=matrix(n)

x=[Decimal(0)]*n
for i in range(0,n):
	if i%2==0: x[i]=Decimal(1)
	else: x[i]=Decimal(-1)

d=calc_d(a,b,c,x,n)

start=Decimal(perf_counter())
xt=thomas(a,b,c,d,n)
end=Decimal(perf_counter())

tt=end-start
e2=euklides(x,xt,n)

A=fullmatrix(a,b,c,n)

start=Decimal(perf_counter())
xg=gauss(A,d,n)
end=Decimal(perf_counter())

tg=end-start
e1=euklides(x,xg,n)

print("Euklide's norm of difference between x given and calculated using gaussian elimination:")
print(e1)
print("Euklide's norm of difference between x given and calculated using Thomas method:")
print(e2)

print("Calculation time when using gaussian elimination:")
print(tg)
print("Calculation time when using Thomas method:")
print(tt)