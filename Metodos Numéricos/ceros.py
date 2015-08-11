from sage.all import *

def biseccio(f, a, b, e=1e-6, NMAX = infinity):
	"""Devuelve el cero de f, si existe, que hay en el intervalo [a,b] con un error mas pequenyo que e.
	El algoritmo hace, como maximo, NMAX iteraciones.

	EJEMPLO:
	>>> f(x) = x
	>>> biseccio(f, -2, 5)
	2.98023223876953e-7
	>>> biseccio(f, -2, 5, 0.1)
	-0.00390625"""
	i = 0
	while i < NMAX:
		if(f(a) == 0):
			return a
		elif(f(b) == 0):
			return b
		elif(b-a < e):
			return (a+b)/2
		else:
			c = (a+b)/2
			if(f(a)*f(c) < 0):
				b = c
			else:
				a = c
		i = i+1
	return (a+b)/2

class PuntsFixes:
	"""
	Esta clase intentara buscar un punto fijo de una funcion.

	CONSTRUCTOR:
	* PuntsFixes(f)
	Queremos buscar un punto fijo de f
	
	EJEMPLO:
	>>> f(x) = x^2   # Punto fijo en x = 0
	>>> A = PuntsFixes(f)
	>>> A.cercar_zero(0.9)
	3.73391848741029e-24
	>>> A.cercar_zero_llistat(0.9)
	[0.900000000000000, 0.810000000000000, 0.656100000000000,
	0.430467210000000, 0.185302018885184, 0.0343368382029252,
	0.00117901845777386, 1.39008452377146e-6, 1.93233498322892e-12,
	3.73391848741029e-24]
	"""
	def __init__(self, f = None):
		self.iteracio = f

	def cercar_puntfix(self, x, e=1e-6, NMAX = infinity):
		"""
		Busca un punto fijo de la funcion de iteracion.
		
		Para ello, a partir del x_0 que recibe como parametro construye la sucesion {x_i} con la siguiente regla
			x_k = f(x_{k-1})
		Hasta que llega un momento en que la distancia entre un iterado y el anterior es mas pequenya que epsilon.
		"""
		i = 0
		xk = x
		while i < NMAX:
			xk1 = self.iteracio(xk)
			if(abs(xk1-xk) < e):
				return xk1
			else:
				xk = xk1
			i = i+1
		return xk

	def cercar_puntfix_llistat(self, x, e=1e-6, NMAX = infinity):
		"""
		Ver cercar_zero(). 
		Este metodo devuelve todas las iteraciones por las que ha pasado.
		"""
		i = 0
		xk = [x]
		while i < NMAX:
			xk.append(self.iteracio(xk[i]))
			if(abs(xk[i+1]-xk[i]) < e):
				return xk
			else:
				i = i+1
		return xk

class Newton(PuntsFixes):
	"""
	Esta clase intentara buscar un cero de la funcion f.
	Para ello generamos la funcion de iteracion
		g(x) = x - m*f(x)/f'(x)
	Y hacemos uso de la clase PuntsFixes para buscar un punto fijo de g.

	CONSTRUCTORES:
	* Newton(f)
	Buscaremos un cero de la funcion f.
	* Newton(f, Df = g)
	Buscaremos un cero de la funcion f, pero se le puede especificar que use como derivada la funcion g.
	* Newton(f, m = k)
	Buscaremos un cero de la funcion f con la constante multiplicativa m como k.

	EJEMPLO:
	>>> f(x) = (x-1)^2
	>>> A = Newton(f)
	>>> A.cercar_zero(0)
	0.999999046325684
	>>> len(A.cercar_zero_llistat(0))
	21
	>>> A = Newton(f, m = 2)
	>>> A.cercar_zero(0)
	1
	>>> len(A.cercar_zero_llistat(0))
	3
	"""
	def __init__(self, f, Df = None, m = 1):
		self.funcion = f
		self.m = m
		if Df != None:
			self.derivada = Df
		else:
			self.derivada = f.diff()

	def iteracio(self, x):
		if (self.funcion(x) == 0):
			return x
		return n(x - self.m*self.funcion(x)/self.derivada(x))

	def cercar_zero(self, x):
		"""
		Busca un cero de la funcion f.

		Vease PuntsFixes.cercar_puntfix
		"""
		return self.cercar_puntfix(x)

	def cercar_zero_llistat(self, x):
		"""
		Busca un cero de la funcion f.
		Devuelve la lista de iterados.

		Vease PuntsFixes.cercar_puntfix_llistat
		"""
		return self.cercar_puntfix_llistat(x)

class Secant:
	"""
	Esta clase intentara buscar un cero de la funcion f mediante el metodo de la secante.

	CONSTRUCTORES:
	* Secant(f)
	Recibimos como parametro la funcion de la que queremos encontrar el 0

	EJEMPLO:
	>>> f(x) = (x-1)^2
	>>> A = Secant(f)
	>>> A.cercar_zero(f, 2,3)
	1.00000107498100
	"""
	def __init__(self, f):
		self.funcion = f

	def iteracio(self, xn, x1):
		return n(x1 - self.funcion(x1)*(xn-x1)/(self.funcion(xn)-self.funcion(x1)))
	
	def cercar_zero(self, x, y, e=1e-6, NMAX = infinity):
		"""
		Busca un cero de la funcion f.
		Recibe como parametro dos puntos iniciales, el error esperado y el maximo de iteraciones permitidos.
		"""
		i = 0
		x1 = x
		xn = y
		xk = y
		while i < NMAX:
			if(self.funcion(xk) == 0):
				return xk
			xk1 = n(self.iteracio(xn, x1))
			if(abs(xk1-xk) < e):
				return xk1
			else:
				x1 = xn
				xn = xk1
				xk = xk1
			i = i+1
		return xk

	def cercar_zero_llistat(self, x, y, e=1e-6, NMAX = infinity):
		"""
		Busca un cero de la funcion f.
		Devuelve una lista con todos los iterados.
		"""
		i = 1
		xk = [x, y]
		while i < NMAX+1:
			if(self.funcion(xk[i]) == 0):
				return xk
			xk.append(n(self.iteracio(xk[i], xk[i-1])))
			if(abs(xk[i+1]-xk[i]) < e):
				return xk
			else:
				i = i+1
		return xk

class NNewton:
    def __init__(self, f, variables):
        self.funciones = f
        self.variables = variables
        
    def derivada(self, punto):
        A = matrix(RR, len(self.funciones))
        for i in range(len(self.funciones)):
            valor = jacobian(self.funciones[i], self.variables)
            for j in range(len(self.variables)):
                valor = valor.subs(self.variables[j] == punto[j])
            A[i] = valor
        return A
        
    def funcion(self, punto):
        resultado = []
        for i in range(len(self.funciones)):
            valor = self.funciones[i]
            for j in range(len(self.variables)):
                valor = valor.subs(self.variables[j] == punto[j])
            resultado.append(valor())
        return vector(resultado)
        
    def iteracio(self, punto):
        DH = self.derivada(punto)
        H = -self.funcion(punto)
        return DH.solve_left(H)
        
    def cercar_zero(self, punto, delta = 1e-6, NMAX = Infinity):
	try:
            punto = vector(punto)
	except:
	    punto = vector([punto])
        d = self.iteracio(punto)
        i = 0
        while i < NMAX and d.norm() > delta:
            punto = punto + d
            d = self.iteracio(punto)
            i = i+1
        return punto
