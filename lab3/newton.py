from sympy import *

def ff(x):
	f=50*x*E**(-50)-50*E**(-10*x)+1/500
	f1=diff(f,x)
	f2=diff(f1,x)
	return f,f1,f2

def newton(f,f1,f2,a,b,p,task,nr):
	if (f2.evalf(subs={x: a})>0 and f.evalf(subs={x: a})>0) or (f2.evalf(subs={x: a})<0 and f.evalf(subs={x: a})<0): xx=a
	else: xx=b

	it=0
	stop=0
	while it==0 or stop>=p:
		xx2=xx-f.evalf(subs={x: xx})/f1.evalf(subs={x: xx})
		it+=1
		if task==1: stop=abs(xx2-xx)
		else: stop=abs(f.evalf(subs={x: xx2}))
		xx=xx2
	return (xx,it)

x=symbols('x')
(f,f1,f2)=ff(x)
print("Function:")
print(f)
print("p parameter:")
p=float(input())
print("")
for i in range(0,1):
	#print("Number of experiment:")
	#print(i+1)
	#print("")
	print("Iterations limiting criteria: |x^(i+1)-x^(i)|<p")
	x1,it=newton(f,f1,f2,-0.5,1.5,p,1,i)
	print("Root in an interval <-0.5,1.5>:")
	print(x1)
	print("Number of iterations:")
	print(it)
	print("")

	print("Iterations limiting criteria: |f(x^(i)|<p")
	x2,it=newton(f,f1,f2,-0.5,1.5,p,2,i)
	print("Root in an interval <-0.5,1.5>:")
	print(x2)
	print("Number of iterations:")
	print(it)
	if i!=5: print("")
