#!/usr/bin/env python
'''
Created Jan 24, 2019 by Casey Berger
Last edited June 4, 2020 by Casey Berger

Notes:
computes the analytical solutions available to us for the noninteracting, nonrotating case
'''

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import cmath
from pylab import meshgrid

def lambda_k(dtau, dim, m, mu, k):
	if dim == 1:
		k0 = k[0]
		kx = k[1] 
		l_Im = np.exp(dtau*mu)*np.sin(k0)
		#l_Re = 1. + dtau*float(dim)/float(m) - np.exp(dtau*mu)*np.cos(k0) 
		l_Re = 1. - np.exp(dtau*mu)*np.cos(k0)
		l_Re = l_Re + (dtau/float(m))*(1 - np.cos(kx))
		return complex(l_Re,l_Im)
	elif dim == 2:
		k0 = k[0]
		kx = k[1] 
		ky = k[2]
		l_Im = np.exp(dtau*mu)*np.sin(k0)
		#l_Re = 1. + dtau*float(dim)/float(m) - np.exp(dtau*mu)*np.cos(k0) 
		l_Re = 1. - np.exp(dtau*mu)*np.cos(k0) 
		l_Re = l_Re + (dtau/float(m))*(1 - np.cos(kx))
		l_Re = l_Re + (dtau/float(m))*(1 - np.cos(ky))
		return complex(l_Re,l_Im)
	elif dim == 3:
		k0 = k[0]
		kx = k[1] 
		ky = k[2]
		kz = k[3]
		l_Im = np.exp(dtau*mu)*np.sin(k0)
		#l_Re = 1. + dtau*float(dim)/float(m) - np.exp(dtau*mu)*np.cos(k0) 
		l_Re = 1. - np.exp(dtau*mu)*np.cos(k0) 
		l_Re = l_Re + (dtau/float(m))*(1 - np.cos(kx))
		l_Re = l_Re + (dtau/float(m))*(1 - np.cos(ky))
		l_Re = l_Re + (dtau/float(m))*(1 - np.cos(kz))
		return complex(l_Re,l_Im)

def dD_dmu(k0,mu,dtau):
	dDdmu = complex((-np.exp(mu*dtau)*np.cos(k0)),(np.exp(mu*dtau)*np.sin(k0)))
	return dDdmu

def dD_dtau(dim,k0,k,mu,dtau,m):
	if dim == 1:
		k0 = k[0]
		kx = k[1] 
		dD_Im = mu*np.exp(mu*dtau)*np.sin(k0)
		dD_Re = -1.*mu*np.exp(mu*dtau)*np.cos(k0) + (1. - np.cos(kx))/m
		return complex(dD_Re,dD_Im)
	elif dim == 2:
		k0 = k[0]
		kx = k[1] 
		ky = k[2]
		dD_Im = mu*np.exp(mu*dtau)*np.sin(k0)
		dD_Re = -1.*mu*np.exp(mu*dtau)*np.cos(k0) + (2. - np.cos(kx) - np.cos(ky))/m
		return complex(dD_Re,dD_Im)
	elif dim == 3:
		k0 = k[0]
		kx = k[1] 
		ky = k[2]
		kz = k[3]
		dD_Im = mu*np.exp(mu*dtau)*np.sin(k0)
		dD_Re = -1.*mu*np.exp(mu*dtau)*np.cos(k0) + (3. - np.cos(kx) - np.cos(ky) - np.cos(kz))/m
		return complex(dD_Re,dD_Im)
'''
def U_k(Nx, Nt, dim, space, t, k):
	volume = float(Nt)*float(np.power(Nx+1,dim))
	if dim == 1:
		k0 = k[0]
		kx = k[1] 
		x = space[0]
		Ux = np.sqrt(np.power(2.,dim))/np.sqrt(float(volume))*np.sin(kx*x)
		U = [Ux*np.cos(k0*t),Ux*np.sin(k0*t)]
		return U
	if dim == 2:
		k0 = k[0]
		kx = k[1]
		ky = k[2]
		x = space[0]
		y = space[1]
		Ux = np.sqrt(np.power(2.,dim))/np.sqrt(float(volume))*np.sin(kx*x)*np.sin(ky*y)
		U = [Ux*np.cos(k0*t),Ux*np.sin(k0*t)]
		return U
	if dim == 3:
		k0 = k[0]
		kx = k[1]
		ky = k[2]
		kz = k[3]
		x = space[0]
		y = space[1]
		z = space[2]
		Ux = np.sqrt(np.power(2.,dim))/np.sqrt(float(volume))*np.sin(kx*x)*np.sin(ky*y)*np.sin(kz*z)
		U = [Ux*np.cos(k0*t),Ux*np.sin(k0*t)]
		return U

def Udagger_k(Nx, Nt, dim, space, t, k):
	volume = float(Nt)*float(np.power(Nx+1,dim))
	if dim == 1:
		k0 = k[0]
		kx = k[1] 
		x = space[0]
		Ux = np.sqrt(np.power(2.,dim))/np.sqrt(float(volume))*np.sin(kx*x)
		U = [Ux*np.cos(k0*t),-1.*Ux*np.sin(k0*t)]
		return U
	if dim == 2:
		k0 = k[0]
		kx = k[1]
		ky = k[2]
		x = space[0]
		y = space[1]
		Ux = np.sqrt(np.power(2.,dim))/np.sqrt(float(volume))*np.sin(kx*x)*np.sin(ky*y)
		U = [Ux*np.cos(k0*t),-1.*Ux*np.sin(k0*t)]
		return U
	if dim == 3:
		k0 = k[0]
		kx = k[1]
		ky = k[2]
		kz = k[3]
		x = space[0]
		y = space[1]
		z = space[2]
		Ux = np.sqrt(np.power(2.,dim))/np.sqrt(float(volume))*np.sin(kx*x)*np.sin(ky*y)*np.sin(kz*z)
		U = [Ux*np.cos(k0*t),-1.*Ux*np.sin(k0*t)]
		return U
'''
def noninteracting_phisq(dim,m,mu,Nx,Nt,dtau):
	NxMin = 1 #previously was 1
	NxMax = Nx+2 #previously was Nx+2
	volume = float(Nt)*float(np.power(Nx+1,dim))
	if dim==1:
		phisq = complex(0.,0.)
		count = 0
		for n0 in range(0,Nt):
			k0 = 2.*np.pi*float(n0)/float(Nt)
			for nx in range(NxMin,NxMax):
				kx = np.pi*float(nx)/float(Nx+1)
				k = [k0,kx]
				D = lambda_k(dtau, dim, m, mu, k)
				phisq += D.conjugate()/(abs(D)*abs(D))
				count += 1
		#print("%s diagonals in Dkk"%count)
		return phisq/volume
	elif dim==2:
		phisq = complex(0.,0.)
		count = 0
		for n0 in range(0,Nt):
			k0 = 2.*np.pi*float(n0)/float(Nt)
			for nx in range(NxMin,NxMax):
				kx = np.pi*float(nx)/float(Nx+1)
				for ny in range(NxMin,NxMax):
					ky = np.pi*float(ny)/float(Nx+1)
					k = [k0,kx,ky]
					D = lambda_k(dtau, dim, m, mu, k)
					phisq += D.conjugate()/(abs(D)*abs(D))
					count += 1
		#print("%s diagonals in Dkk"%count)
		return phisq/volume
	elif dim==3:
		phisq = complex(0.,0.)
		count = 0
		for n0 in range(0,Nt):
			k0 = 2.*np.pi*float(n0)/float(Nt)
			for nx in range(NxMin,NxMax):
				kx = np.pi*float(nx)/float(Nx+1)
				for ny in range(NxMin,NxMax):
					ky = np.pi*float(ny)/float(Nx+1)
					for nz in range(NxMin,NxMax):
						kz = np.pi*float(nz)/float(Nx+1)
						k = [k0,kx,ky,kz]
						D = lambda_k(dtau, dim, m, mu, k)
						phisq += D.conjugate()/(abs(D)*abs(D))
						count += 1
		#print("%s diagonals in Dkk"%count)
		#print("Volume = %s"%(Nt*Nx*Nx*Nx))
		return phisq/volume
	else:
		print("invalid dimension, D = %s"%dim)
		sys.exit()

def noninteracting_density(dim,m,mu,Nx,Nt,dtau):
	NxMin = 1 #previously was 1
	NxMax = Nx+2 #previously was Nx+2
	volume = float(Nt)*float(np.power(Nx+1,dim))
	if dim==1:
		density = complex(0.,0.)
		count = 0
		for n0 in range(0,Nt):
			k0 = 2.*np.pi*float(n0)/float(Nt)
			for nx in range(NxMin,NxMax):
				kx = np.pi*float(nx)/float(Nx+1)
				k = [k0,kx]
				D = lambda_k(dtau, dim, m, mu, k)
				dDdmu = dD_dmu(k0,mu,dtau)
				density += (D.conjugate()*dDdmu)/(abs(D)*abs(D))
				count += 1
				#print(str(D)+"\t"+str(dDdmu)+"\t"+str(abs(D))+"\t"+str(density/(Nx*Nt))+"\n")
				#print(str(D)+"\t"+str(dDdmu)+"\t"+str(density)+"\n")
		#print("%s diagonals in Dkk"%count)
		#return density/(float(Nt)*float(np.power(Nx,d)))
		return -1.*density/volume
	elif dim==2:
		density = complex(0.,0.)
		count = 0
		for n0 in range(0,Nt):
			k0 = 2.*np.pi*float(n0)/float(Nt)
			for nx in range(NxMin,NxMax):
				kx = np.pi*float(nx)/float(Nx+1)
				for ny in range(NxMin,NxMax):
					ky = np.pi*float(ny)/float(Nx+1)
					k = [k0,kx,ky]
					D = lambda_k(dtau, dim, m, mu, k)
					dDdmu = dD_dmu(k0,mu,dtau)
					density += (D.conjugate()*dDdmu)/(abs(D)*abs(D))
					count += 1
		#print("%s diagonals in Dkk"%count)
		return -1.*density/volume
	elif dim==3:
		density = complex(0.,0.)
		count = 0
		for n0 in range(0,Nt):
			k0 = 2.*np.pi*float(n0)/float(Nt)
			for nx in range(NxMin,NxMax):
				kx = np.pi*float(nx)/float(Nx+1)
				for ny in range(NxMin,NxMax):
					ky = np.pi*float(ny)/float(Nx+1)
					for nz in range(NxMin,NxMax):
						kz = np.pi*float(nz)/float(Nx+1)
						k = [k0,kx,ky,kz]
						D = lambda_k(dtau, dim, m, mu, k)
						dDdmu = dD_dmu(k0,mu,dtau)
						density += (D.conjugate()*dDdmu)/(abs(D)*abs(D))
						#print(str(D)+"\t"+str(dDdmu)+"\t"+str(abs(D))+"\t"+str(density/(Nx*Nx*Nx*Nt))+"\n")
						count += 1
		#print("%s diagonals in Dkk"%count)
		#print("Volume = %s"%(Nt*Nx*Nx*Nx))
		return -1.*density/volume
	else:
		print("invalid dimension, D = %s"%dim)
		sys.exit()

def noninteracting_energy(dim,m,mu,Nx,Nt,dtau):
	NxMin = 1 #previously was 1
	NxMax = Nx+2 #previously was Nx+2
	volume = float(Nt)*float(np.power(Nx+1,dim))
	if dim==1:
		energy = complex(0.,0.)
		count = 0
		for n0 in range(0,Nt):
			k0 = 2.*np.pi*float(n0)/float(Nt)
			for nx in range(NxMin,NxMax):
				kx = np.pi*float(nx)/float(Nx+1)
				k = [k0,kx]
				D = lambda_k(dtau, dim, m, mu, k)
				dDdtau = dD_dtau(dim,k0,k,mu,dtau,m)
				energy += (D.conjugate()*dDdtau)/(abs(D)*abs(D))
				count += 1
		return energy/volume
	elif dim==2:
		energy = complex(0.,0.)
		count = 0
		for n0 in range(0,Nt):
			k0 = 2.*np.pi*float(n0)/float(Nt)
			for nx in range(NxMin,NxMax):
				kx = np.pi*float(nx)/float(Nx+1)
				for ny in range(NxMin,NxMax):
					ky = np.pi*float(ny)/float(Nx+1)
					k = [k0,kx,ky]
					D = lambda_k(dtau, dim, m, mu, k)
					dDdtau = dD_dtau(dim,k0,k,mu,dtau,m)
					energy += (D.conjugate()*dDdtau)/(abs(D)*abs(D))
					count += 1
		return energy/volume
	elif dim==3:
		energy = complex(0.,0.)
		count = 0
		for n0 in range(0,Nt):
			k0 = 2.*np.pi*float(n0)/float(Nt)
			for nx in range(NxMin,NxMax):
				kx = np.pi*float(nx)/float(Nx+1)
				for ny in range(NxMin,NxMax):
					ky = np.pi*float(ny)/float(Nx+1)
					for nz in range(NxMin,NxMax):
						kz = np.pi*float(nz)/float(Nx+1)
						k = [k0,kx,ky,kz]
						D = lambda_k(dtau, dim, m, mu, k)
						dDdtau = dD_dtau(dim,k0,k,mu,dtau,m)
						energy += (D.conjugate()*dDdtau)/(abs(D)*abs(D))
						count += 1
		return energy/volume
	else:
		print("invalid dimension, D = %s"%dim)
		sys.exit()