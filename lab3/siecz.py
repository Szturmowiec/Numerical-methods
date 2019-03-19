from math import exp

def f(x):
	wyn=50*x*exp(-50)-50*exp(-10*x)+1/500
	return wyn

def siecz(a,b,p,task):
	x0=a+0.2
	x1=b-0.2

	it=0
	stop=0
	while it==0 or stop>=p:
		if f(x1)==f(x0):
			print("Cannot continue calculations, because f(x^(i+1))=f(x^(i))")
			break
		x2=x1-(f(x1)*(x1-x0))/(f(x1)-f(x0))

		it+=1
		if task==1: stop=abs(x2-x1)
		else: stop=abs(f(x2))
		x0=x1
		x1=x2
	return (x1,it)

print("Function:")
print("50*x*exp(-50) + 0.002 - 50*exp(-10*x)")
print("p parameter:")
p=float(input())
print("")

print("Iterations limiting criteria: |x^(i+1)-x^(i)|<p")
x1,it=siecz(-0.5,1.5,p,1)
print("Root in an interval <-0.5,1.5>:")
print(x1)
print("Number of iterations:")
print(it)
print("")

print("Iterations limiting criteria: |f(x^(i)|<p")
x2,it=siecz(-0.5,1.5,p,2)
print("Root in an interval <-0.5,1.5>:")
print(x2)
print("Number of iterations:")
print(it)