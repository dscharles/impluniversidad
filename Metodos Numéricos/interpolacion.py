# coding: utf-8

from sage.all import *

class Polinomio:
    def __init__(self, puntos):
        self.puntos = puntos
        
    def __call__(self, x):
        i = 1
        x = n(x)
        if x <= self.puntos[0][0]:
            return self.puntos[0][1]
        while i < len(self.puntos):
            if self.puntos[i-1] <= x <= self.puntos[i][0]:
                return n(self.puntos[i-1][1]+(x-self.puntos[i-1][0])/(self.puntos[i][0]-self.puntos[i-1][0])*(self.puntos[i][1]-self.puntos[i-1][1]))
            i = i+1
        return self.puntos[-1][1]
        
    def __float__(self, x):
        return self(x)

class PolinomioLagrange:
	"""
	Genera el unico polinomio de grado n tal que
		p(x_i) = f(x_i)   i = 0, ..., n

	CONSTRUCTOR:
	* PolinomioLagrange()
	Genera el objeto vacio. Despues se puede anyadir puntos mediante el metodo add_point.

	* PolinomioLagrange(valores)
	valores -- Una lista con los puntos por los que el polinomio debe pasar. 
	Por ejemplo: [ [0, 0], [1,1] ] genera la recta que pasa por los puntos (0,0) y (1,1)

	EJEMPLO:
	>>> val = [ [i, i] for i in range(2) ]; val
	[ [0,0], [1,1] ]
	>>> p = PolinomioLagrange(val)
	>>> p(0.5)
	0.5

	EJEMPLO:
	>>> p = PolinomioLagrange()
	>>> p(0)
	None
	>>> p.add_point(0,0)
	>>> p.add_point(1,1)
	>>> p.add_point(2,8)
	>>> p(10)
	280
	>>> p.add_point(3,27)
	>>> p(10)
	1000
	>>> p.poli
	x^3
	>>> p.diff()
	3*x^2
	"""
	def __init__(self, valores = []):
		self.v = valores
		self.dict = dict()
		self.poli = None
		R = PolynomialRing(RR, 'x')
		self.x = R.gen()

	def _tostr_(self, v):
		a = ''
		for p in v:
			a = a + '-' + str(p)
		return a

	def diferencias_divididas(self, val):
		"""
		Calcula las diferencias divididas de [ [x_1, f(x_1)], ..., [x_k, f(x_k)] ]
		"""
		try:
			return self.dict[self._tostr_(val)]
		except:
			if len(val) == 1:
				self.dict[self._tostr_(val)] = val[0][1]
			else:
				self.dict[self._tostr_(val)] = (self.diferencias_divididas(val[1:len(val)]) - self.diferencias_divididas(val[0:len(val)-1]))/(val[-1][0] - val[0][0])
			return self.dict[self._tostr_(val)]

	def generar_polinomio(self):
		"""
		Genera el polinomio de grado n que pasa por los puntos (x_i, f(x_i))  i = 0, ..., n
		Si despues se quiere acceder al polinomio, se guarda en la variable poli
		"""
		self.poli = 0
		for i in range(len(self.v)):
			poli2 = n(self.diferencias_divididas(self.v[0:i+1]))
			for j in range(i):
				poli2 *= self.x-self.v[j][0]
			self.poli = self.poli + poli2

	def __call__(self, x):
		if len(self.v) < 2:
			return None
		if self.poli == None:
			self.generar_polinomio()
		return self.poli(x)

	def __float__(self,x):
		return self(x)

	def diff(self, x = None):
		"""
		Devuelve la derivada del polinomio. Si se especifica el punto x, se devuelve la funcion evaluada en ese punto.
		Devuelve None si no hay los suficientes puntos.
		"""
		if len(self.v) < 2:
			return None
		if self.poli == None:
			self.generar_polinomio()
		if x != None:
			return diff(self.poli)(x)
		return diff(self.poli)

	def add_point(self, x, fx):
		"""
		Anyade el punto (x, fx) a la lista de puntos por el que tiene que pasar el polinomio. Se tiene que volver a generar el polinomio ya que por defecto no es generado.
		"""
		self.v.append([x,fx])
		self.poli = None

class PolinomioHermite:
	"""
	Genera el unico polinomio tal que
		p(x_i) = f_i     i = 0, ..., n
		p'(x_i) = d_i    i = 0, ..., n

	CONSTRUCTORES:
	* PolinomioHermite()
	Genera el objeto vacio. Despues se puede anyadir los puntos mediante el metodo add_point()

	* PolinomioHermite(valores)
	valores -- Una lista con los puntos por los que debe pasar el polinomio, y su derivada.
	Por ejemplo: [ [0,0,0], [1,1,2] ] genera la funcion x^2

	EJEMPLO:
	>>> val = [ [i, i^2, 2*i] for i in range(3) ]; val
	[ [0,0,0], [1,1,2], [2,4,4] ]
	>>> p = PolinomioHermite(val)
	>>> p(10)
	100

	EJEMPLO:
	>>> p = PolinomioHermite()
	>>> p.add_point(0,0,0)
	>>> p.add_point(1,1,2)
	>>> p.generar_polinomio()
	>>> p.poli
	x^2
	>>> p(10)
	100
	>>> p.diff()
	2*x
	"""
	def __init__(self, valores = []):
		self.v = []
		for p in valores:
			self.v.append(p)
			self.v.append(p)
		self.dict = dict()
		self.poli = None
		R = PolynomialRing(RR, 'x')
		self.x = R.gen()
	
	def _tostr_(self, v):
		a = ''
		for p in v:
			a = a + '-' + str(p)
		return a

	def add_point(self, x, fx, dx):
		"""
		Anyade el punto (x, fx) por el que tiene que pasar el polinomio, ademas obliga a que la derivada en x sea dx
		"""
		self.v.append([x, fx, dx])
		self.v.append([x, fx, dx])
		self.poli = None

	def diferencias_divididas(self, val):
		try:
			return self.dict[self._tostr_(val)]
		except:
			if len(val) == 1:
				self.dict[self._tostr_(val)] = val[0][1]
			elif len(val) == 2 and val[0] == val[1]:
				self.dict[self._tostr_(val)] = val[0][2]
			else:
				self.dict[self._tostr_(val)] = (self.diferencias_divididas(val[1:len(val)]) - self.diferencias_divididas(val[0:len(val)-1]))/(val[-1][0] - val[0][0])
			return self.dict[self._tostr_(val)]

	def generar_polinomio(self):
		"""
		Genera el polinomio interpolador. Si después se quiere hacer uso del polinomio se guarda en la variable poli
		"""
		self.poli = 0
		for i in range(len(self.v)):
			poli2 = self.diferencias_divididas(self.v[0:i+1])
			for j in range(i):
				poli2 = poli2*(self.x-self.v[j][0])
			self.poli = self.poli + poli2
	
	def __call__(self, x):
		if len(self.v) < 4:
			return None
		if self.poli == None:
			self.generar_polinomio()
		return self.poli(x)

	def __float__(self, x):
		return self(x)

	def diff(self):
		"""
		Devuelve la derivada del polinomio. Si se especifica el punto x, se devuelve la funcion evaluada en ese punto.
		Devuelve None si no hay suficientes puntos.
		"""
		if len(self.v) < 4:
			return None
		if self.poli == None:
			self.generar_polinomio()
		if x != None:
			return diff(self.poli)(x)
		return diff(self.poli)

class Spline:
	"""
	Genera un Spline a partir de una lista de puntos por los que tiene que pasar.

	El constructor acepta una lista de puntos de R^2 por los que tiene que pasar la funcion, 
	el formato es el siguiente
		[ [x_0, f(x_0)], [x_2, f(x_2)], ..., [x_n, f(x_n)] ]

	Despues de haber creado el Spline se puede anyadir puntos mediante el metodo add_point

	EJEMPLO:
	>>> val = [ [-1, 1], [0,0], [1,1] ]  # f(x) = x^2
	>>> S = Spline(val)
	>>> S(-1) == S(1) == 1
	True
	>>> S(0)
	0
	>>> S.diff(0)
	0
	>>> S(0.5)
	0.3125
	>>> S.add_point(0.5, 1.5)
	>>> S(0.5)
	1.5
	"""
	def __init__(self, valores = []):
		self.N = len(valores) -1
		self.Sp = matrix(RR,self.N, 6)
		self.puntos = [ p[0] for p in valores ]
		self.imagenes = [ p[1] for p in valores ]
		self.rellenado = False

	def add_point(self, x, fx):
		"""
		Anyade el punto (x, fx) a la lista de puntos por los que tiene que pasar el Spline.
		"""
		i = 0
		while i < self.N:
			if self.puntos[i] < x < self.puntos[i+1]:
				self.puntos[i+1:i+1] = [x]
				self.imagenes[i+1:i+1] = [fx]
				self.N = self.N + 1
				self.rellenado = False
				i = self.N
			i = i+1

	def _fill_(self):
		if self.N < 2:
			return None
		puntos = self.puntos
		imagenes = self.imagenes
		self.rellenado = True
		self.Sp = matrix(RR, self.N, 6)
		for i in range(self.N):
			self.Sp[i,0] = puntos[i]			# Inicio intervalo
			self.Sp[i,1] = imagenes[i]			# Imagen del extremo inicial
			self.Sp[i,5] = puntos[i+1] - puntos[i]		# Longitud del intervalo
		A = matrix(RR, self.N+1)
		A[0,0] = 1
		A[self.N, self.N] = 1
		for i in range(1, self.N):
			A[i,i] = 2*(self.Sp[i-1,5] + self.Sp[i,5])
			A[i,i-1] = self.Sp[i-1,5]
			A[i,i+1] = self.Sp[i,5]
		b = []
		b.append(0)
		for i in range(1, self.N):
			b.append((3/self.Sp[i,5])*(imagenes[i+1] - imagenes[i]) - (3/self.Sp[i-1,5])*(imagenes[i] - imagenes[i-1]))
		b.append(0)
		c = A.solve_right(vector(b))
		for i in range(self.N):
			self.Sp[i,2] = (imagenes[i+1]-imagenes[i])/(self.Sp[i, 5])-self.Sp[i,5]*(c[i+1]+2*c[i])/(3)
			self.Sp[i,3] = c[i]
			self.Sp[i,4] = (c[i+1]-c[i])/(3*self.Sp[i,5])

	def __call__(self, punto):
		if self.N < 2:
			return None
		if not self.rellenado:
			self._fill_()
		if(punto < self.Sp[0,0]):
			return self.Sp[0,1]
		if(punto > self.Sp[self.N-1, 0] + self.Sp[self.N-1,5]):
			return (self.Sp[self.N-1,1] + 
					self.Sp[self.N-1,2]*self.Sp[self.N-1,5] + 
					self.Sp[self.N-1,3]*self.Sp[self.N-1,5]**2 + 
					self.Sp[self.N-1,4]*self.Sp[self.N-1,5]**3 )
		j = 0
		while True:
			if(punto <= self.Sp[j, 0] + self.Sp[j, 5]):
				return (self.Sp[j,1] + self.Sp[j,2]*(punto - self.Sp[j,0]) + self.Sp[j,3]*((punto - self.Sp[j,0])**2) + self.Sp[j,4]*((punto - self.Sp[j,0]))**3 )
			j = j+1

	def __float__(self, punto):
		return self(punto)

	def diff(self, punto):
		if self.N < 2:
			return None
		if not self.rellenado:
			self._fill_()
		if(punto < self.Sp[0,0]):
			return selfw.Sp[0,2]
		if(punto > self.Sp[self.N-1, 0] + self.Sp[self.N-1,5]):
			return (self.Sp[self.N-1,2] + 
					2*self.Sp[self.N-1,3]*self.Sp[self.N-1,5] + 
					3*self.Sp[self.N-1,4]*self.Sp[self.N-1,5]**2 )
		j = 0
		while True:
			if(punto <= self.Sp[j, 0] + self.Sp[j, 5]):
				return (self.Sp[j,2] + 2*self.Sp[j,3]*(punto - self.Sp[j,0]) + 3*self.Sp[j,4]*(punto - self.Sp[j,0])**2 )
			j = j+1
