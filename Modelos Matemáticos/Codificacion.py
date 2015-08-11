from sage.all import *

from Laboratorio.Tema1 import *
from Laboratorio.Tema2 import *

#
# C = Codi(p, m)
# C.gf(polinomio(C.x))
#

class Codi:
	def __init__(self, p, m):
		self.p = p
		self.m = m
		self.n = p**m-1
		R = PolynomialRing(GF(p, 'x'), 'x')
		self.x = R.gen()

	def gf(self, modulo):
		self.GF = GF(self.p**self.m, 'alpha', modulus=modulo)
		self.alpha = self.GF.gen()
		self.Polynomial = PolynomialRing(self.GF, 'y')
		self.y = self.Polynomial.gen()

        def LogaritmeZech(self):
                print("n \t Ln(n)")
                for k, a in enumerate(self.GF):
                    for j, b in enumerate(self.GF):
                        if a == b+1:
                            print(str(j) + "\t" + str(k))
                            
        def ordre(self, beta):
            	t = 1
           	while beta**t != 1:
                	t += 1
            	return t
            
        def polMinim(self, beta):
            	r = self.ordre(beta)
            	m = 1
           	while (self.p**m-1) % r != 0:
                	m +=1
            	prod = 1
            	for i in range(0, m):
                	prod = prod*(self.y - beta**(self.p**i))
            	return prod
            
        def codificar(self, polinomio):
            	Cociente = self.Polynomial.quotient(self.polinomi, 'y')
            	polinomio = (self.y**self.polinomi.degree()) * polinomio
            	polinomio = polinomio - Cociente(polinomio).lift()
            	return polinomio
            
	def polSindrome(self, polinomio):
		self.sindrome = 0
		for i in range(1,self.d):
			self.sindrome += polinomio(self.alpha**i)*self.y**(i-1)

	def polError(self):
		A = EuclidesEstes(self.y**(self.d-1), self.sindrome)
		i = 0
		while A[1][i].degree() >= int((self.d-1)/2):
			i+=1
		matriz = matrix([[1,0],[0,1]])
		for k in range(i):
			matriz = matrix([[A[4][k], 1],[1,0]]) * matriz
		self.localitzador = matriz[0,1]
		self.difflocalitzador = diff(self.localitzador)
		self.avaluador = (-1)**2*A[1][i]

	def exponente(self,beta):
		for i in range(self.n):
			if self.alpha**i == beta**(-1):
				return i

	def corregirerrores(self, polinomio):
		self.polSindrome(polinomio)
		if self.sindrome == 0:
			return polinomio
		self.polError()
		error = []
		error2 = []
		for i in range(self.n):
			if self.localitzador(self.alpha**i) == 0:
				error.append(self.exponente(self.alpha**i))
				error2.append(i)
		corregir = []
		print error2
		for i in range(len(error)):
			corregir.append(self.avaluador(self.alpha**error2[i])/self.difflocalitzador(self.alpha**error2[i]))
		for i in range(len(corregir)):
			polinomio -= corregir[i]*self.y**error[i]
		print error
		return polinomio

	def decodificar(self, polinomio):
		pol = self.corregirerrores(polinomio)
		valor = []
		coeficientes = pol.coeffs()
		max = len(coeficientes)
		for i in range(self.n-self.k, len(coeficientes)):
			valor.append(coeficientes[i])
		for i in range(len(valor), self.k):
			valor.append(0)
		return valor

class CodiBCH(Codi):
        def __init__(self, p, m, t):
            	Codi.__init__(self,p,m)
            	self.t = t
		self.d = 2*t+1

        def polSindrome(self, polinomio):
            	self.sindrome = polinomio.quo_rem(self.polinomi)[1]
            	return self.sindrome
        
        def polGenerador(self):
            	self.polinomi = 1
            	for i in range(1, 2*self.t, 2):
                	self.polinomi *= self.polMinim(self.alpha**i)
            	self.k = self.n - self.polinomi.degree()
            	return self.polinomi
	
class CodiRS(Codi):
    	def __init__(self, p, m, d):
        	Codi.__init__(self,p,m)
        	self.d = d
        	self.k = self.n-(d-1)
        
   	def polGenerador(self):
        	self.polinomi = 1
        	for i in range(self.d-1):
            		self.polinomi *= (self.y-self.alpha**(1+i))
        	return self.polinomi

	def binary_to_message(self, vector):
		partes = len(vector)/self.m
		info = []
		for i in range(partes):
			# trabajando con la subsucesion
			# a_{partes*i} ... a_{partes*i + m-1}
			elemento = 0
			for k in range(self.m*i, self.m*(i+1)):
				elemento += vector[k]*(self.alpha**(k- self.m*i))
			info.append(elemento)
		return info

	# a_0 + a_1x + a_2x^2 + ... a_nx^n --> [a_0, a_1, ..., a_n]
	def message_to_binary(self,polinomio):			
		vector=[]
		for i in range(len(polinomio)):
			vector2 = []
			elemento = self.GF(polinomio[i]).polynomial()
			for j in range(self.m-1, -1, -1):
				elemento2 = (self.alpha**j).polynomial()
				q = elemento.quo_rem(elemento2)
				if (q[0].coeffs() == []):
					vector2[0:0] = [0]
				else:
					vector2[0:0] = [q[0].coeffs()[0]]
			for k in range(len(vector2)):
				vector.append(vector2[k])
		return vector
