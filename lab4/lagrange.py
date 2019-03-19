from sympy import *
import math as m
import matplotlib.pyplot as plt

def czeb(n,a,b):
	x=[]
	for k in range(n):
		root=m.cos(pi*((2*(k+1)-1)/(2*n)))
		w=((b-a)*root+a+b)/2
		x.append(w)
	return x

def ff(x):
	if x==symbols('x'): return x*sin((0.8*m.pi)/x)
	return x*m.sin((0.8*m.pi)/x)

def inter(x,sym):
	w=0
	for i in range(len(x)):
		l=1
		for j in range(len(x)):
			if i!=j: l*=(sym-x[j])/(x[i]-x[j])
		w+=ff(x[i])*l
	return w

sym=symbols('x')
(a,b)=(0.1,0.8)
print("Number of nodes:")
n=int(input())
x1=[]
for i in range(n): x1.append(a+((b-a)/(n-1)*i))
x2=czeb(n,a,b)

w=inter(x1,sym)
w2=inter(x2,sym)

ile=200

y1=[]
xx1=[]
for i in range(ile):
	xx1.append(a+i*((b-a)/ile))
	y1.append(w.evalf(subs={sym: xx1[i]}))
xx1.append(b)
y1.append(w.evalf(subs={sym: xx1[ile]}))



y2=[]
xx2=[]
for i in range(ile):
	xx2.append(a+i*((b-a)/ile))
	y2.append(w2.evalf(subs={sym: xx2[i]}))
xx2.append(b)
y2.append(w2.evalf(subs={sym: xx2[ile]}))

yorg=[]
xorg=[]
for i in range(ile):
	xorg.append(a+i*((b-a)/ile))
	yorg.append(ff(xorg[i]))
xorg.append(b)
yorg.append(ff(xorg[ile]))

yr=[]
yc=[]
for i in range(n):
	yr.append(ff(x1[i]))
	yc.append(ff(x2[i]))

print("Choose nodes:")
print("1) Equal")
print("2) Chebyshev")
task=int(input())

maks=0
for i in range(ile):
	if task==1:
		if abs(y1[i]-yorg[i])>maks: maks=abs(y1[i]-yorg[i])
	else:
		if abs(y2[i]-yorg[i])>maks: maks=abs(y2[i]-yorg[i])
print("Difference between interpolated function and interpolating polynomial using maximum norm:")
print(maks)

plt.scatter(xorg,yorg,s=5,color='b')
if task==1:
	plt.scatter(xx1,y1,s=5,color='r')
	plt.scatter(x1,yr,s=10,color="black")
else:
	plt.scatter(xx2,y2,s=5,color='r')
	plt.scatter(x2,yc,s=10,color="black")

axes=plt.gca()
axes.set_xlim([0.1,0.8])
axes.set_ylim([-1,1])
plt.show()