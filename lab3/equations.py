from sympy import *
from math import sqrt

def ff(s):
	f=[]
	f.append(s[0]**2-4*s[1]**2+s[2]**3-1)
	f.append(2*s[0]**2+4*s[1]**2-3*s[2])
	f.append(s[0]**2-2*s[1]+s[2]**2-1)
	return f

def euklides(x,n):
	s=0
	for i in range(0,n): s=s+x[i]*x[i]
	s=sqrt(s)
	return s

def solv(f,s,n,start,p,task):
	J=[[0]*n for i in range(n)]
	for i in range(0,n):
		for j in range(0,n):
			J[i][j]=f[i].diff(s[j])

	x=start
	stop=0
	it=0
	while stop>=p or it==0:
		w=[[1]*n for i in range(n)]
		jj=[[1]*n for i in range(n)]

		for i in range(0,n):
			w[i]=f[i].evalf(subs={s[0]: x[0],s[1]: x[1],s[2]: x[2]})
		for i in range(0,n):
			for j in range(0,n):
				jj[i][j]=J[i][j].evalf(subs={s[0]: x[0],s[1]: x[1],s[2]: x[2]})

		jj=Matrix(jj)
		w=Matrix(w)
		x=Matrix(x)
		jj=jj.inv()
		x1=x-jj*w

		if task==1: stop=euklides(x1-x,n)
		else:
			k=[[1]*n for i in range(n)]
			for i in range(0,n):
				k[i]=f[i].evalf(subs={s[0]: x1[0],s[1]: x1[1],s[2]: x1[2]})
			stop=euklides(k,n)

		x=x1
		it+=1
	return (x1,it)

s=[]
s.append(symbols('x1'))
s.append(symbols('x2'))
s.append(symbols('x3'))
f=ff(s)
print("System of equations:")
print("x1^2-4x2^2+x3^3=1")
print("2x1^2+4x2^2-3x3=0")
print("x1^2-2x2+x3^2=1")

print("Starting vector:")
start=list(map(float,input().split(' ')))
print("p parameter:")
p=0.0000001
print(p)
print("")

print("Iterations limiting criteria: ||x^(i+1)-x^(i)||<p")
x,it=solv(f,s,3,start,p,1)
print("Solution:")
print(list(x))
print("Number of iterations:")
print(it)
print("")

print("Iterations limiting criteria: ||F(x^(i))||<p")
x,it=solv(f,s,3,start,p,2)
print("Solution:")
print(list(x))
print("Number of iterations:")
print(it)
