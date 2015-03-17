from numpy import *
from pylab import *


I=complex(0.,1.)

L=1. #length
T=1. #time
x0 = L/4.
sig0=L/20.
lamb=1. #lambda
J=500  # number spatial steps
eps =L/float(J) #spatial step
delt=2.*eps**2/lamb 
N=int(T/delt) #time iterations
k0 = pi/7./eps


#potential
V=zeros(J+1,dtype=complex) 
a1=int((J+1)*(.5-.032))
a2=int((J+1)*(.5+.032))
V[a1:a2]=(1*(50*pi)**2)*ones(a2-a1,dtype=complex)


#initial condition
x=linspace(0,1,J+1)
psi = exp(I*k0*x)*exp(-(x-x0)**2/(2.*sig0**2))


#Build e vector
E=zeros(J+1)*complex(0.,0.)
E[1]=2.+eps**2*V[1]-I*lamb
for i in range(2,J+1):
	E[i]=2.+eps**2*V[i]-I*lamb-(1/E[i-1])


#initializing vectors
Omega=zeros(J+1)*complex(0.,0.)
Psi = zeros(J+1)*complex(0.,0.)
f = zeros(J+1)*complex(0.,0.)
Psi[1:-1]=psi[1:-1]



ion()

p,=plot(x,abs(Psi))
r,=plot(x,real(Psi))
q,=plot(x,V)
#show()
#the iteration
for n in range(0,N):
	Omega[1:-1]=-Psi[2:]-Psi[0:-2]+(I*lamb+2.)*Psi[1:-1]+eps**2*(Psi[1:-1]*V[1:-1])
	f[1]=Omega[1]
	for i in range(2,J+1):
		f[i]=Omega[i]+f[i-1]/E[i-1]
	for i in range(2,J+1):
		Psi[-i]=(Psi[-i+1]-f[-i])/E[-i]
	if i%100==0:
		p.set_ydata(abs(Psi))
		r.set_ydata(real(Psi))
		ylim(-1,1)	
		draw()

ioff()
show()































