from sage.all import *
from Tema1 import *

import random
import math

def descomposicionMR(n):
	b = 0;
	t = n;
	while(t % 2 == 0):
		b = b+1
		t = t/2
	return([b,t])

# p es el posible primo
# e es el numero de elementos de {2, ..., n-1} con 
# los que probaremos el algoritmo
def MillerRabin(p,e):
	if p == 2:
		return [True]
	i = 1
	d = descomposicionMR(p-1)
	while(i<=e):
		i=i+1
		a = random.randint(2, p-1)
		y = MillerRabinR(p,a,d)
		if y == False:
			print "Falla el test con el numero", a
			return ([False, d])
		print "Pasa el test con el numero", a
	return([True,d])

def MillerRabinR(p,a,d):
	y = []
	y.append(ExpRapida(a, d[1], p))
	if y[0] == 1:
		return True
	if y[0] == p-1:
		return True
	for i in range(d[0]):
		y.append(ExpRapida(y[i], 2, p))
		if y[i+1] == p-1:
			return True
	return False
	
# primo_grande(500000, 5000000)
# Encuentra todos los posibles primos fuertes entre 
#   2*min + 1 y 2*max + 1
def primo_grande(min, max, e):
	a = min
	while(a <= max):
		p = MillerRabin(2*a + 1, e)
		# El metodo de MillerRabin ha asegurado que
		#  a es probablemente primo
		#  ademas, a-1 = 2*t
		if (p[0] and p[1][0] == 1):
			# Comprobamos que t es primo con una fiabilidad de 
			#  2^(-17) \approx 10^(-6)
			p = MillerRabin(p[1][1], e)
			if p[0]:
				print ("Primo fuerte:", a, 2*a+1)
		a = a+1

# Intenta calcular un divisor de n
def Pollard(n):
	i = 0
	x = [random.randint(2,n-1)]
	y = [x[0]]
	print(i, x[i], y[i])
	while True:
		i = i+1
		x.append(Modulo(ExpRapida(x[i-1],2,n) + 1, n))
		y.append(Modulo(ExpRapida(ExpRapida(y[i-1],2,n)+1,2,n)+1,n))
		a = gcd(abs(x[i]-y[i]), n)
		print(i, x[i], y[i], a)
		if 1 < a < n:
			return a
		elif a == n:
			return False

# Devuelve un listado con los posibles primos entre
# 1 y (|_ sqrt(n) _|+1)
# Para ver si es primo ejecutamos el test de MillerRabin
def ListaPrimos(n):
	tope = int(math.sqrt(n))+1
	i = 2;
	primos = [];
	while(i<= tope):
		if(MillerRabin(i, 5)[0] == True):
			primos.append(i)
		i=i+1
	return primos

# Calculo de la m
def mPollard(n):
	m = 1
	tope = int(math.sqrt(n))+1
	primos = ListaPrimos(n)
	for p in primos:
		i = 1
		while(p**i <= tope):
			i = i+1
		m = m * p**(i-1)
	return m

# Algoritmo p-1 de Pollard
# e es el numero de condiciones iniciales que usaremos
def Pollard1(n, e):
	i = 1
	m = mPollard(n)
	#print m
	while(i<=e):
		a = random.randint(2, n-1)
		print "Condicion inicial", a
		d = gcd(a,n)
		if d > 1:
			return d
		x = ExpRapida(a, m, n)
		d = gcd(x-1, n)
		#print (a,d)
		if d > 1 and d != n:
			return d
		i=i+1
	return False
