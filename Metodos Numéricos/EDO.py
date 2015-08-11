from sage.all import *
from Metodos.ceros import NNewton

def flecha(punto, pendiente, max, min, longitud):
    angulo = arctan(pendiente)
    ratio = abs(pendiente)/max
    if ratio < 0.1:
        ratio = 0.1
    ve = ratio*vector([cos(angulo), sin(angulo)])
    final = vector(punto) + longitud*ve
    return line([punto, final])

def retrato_fase(f, inicio = [-1, -1], final = [1,1], div = 10):
    p = point( [(inicio[0] + final[0])/2, (inicio[1] + final[1])/2] )
    derivadas = []
    minim = Infinity
    maxim = 0
    longx = (final[0]-inicio[0])/div
    longy = (final[1]-inicio[1])/div
    ancho = min(longx, longy)
    for i in range(div):
        for j in range(div):
            x = inicio[0] + (i+0.5)*longx
            y = inicio[1] + (j+0.5)*longy
            try:
                der = f(x,y)
                derivadas.append([x,y,der])
                minim = min(abs(der), minim)
                maxim = max(abs(der), maxim)
            except:
                pass
    for pareja in derivadas:
        p = p + flecha([pareja[0],pareja[1]], pareja[2], maxim, minim, ancho)
    return p

class RungeKutta:
    def __init__(self, A, b, c):
        self.A = A
        self.b = b
        self.c = c
        for i in range(len(b)):
            for j in range(i,len(b)):
                if A[i][j] != 0:
                    print "El metodo no funcionara, la matriz tiene que ser estrictamente diagonal inferior"
                    
    def orden(self):
        suma = 0
        N = len(self.b)
        for i in range(N):
            suma = suma + self.b[i]
        if not suma == 1:
            print "El metodo no es de orden 1"
            return None
        if N == 1:
            print "El metodo es de orden 1"
            return None
        suma = 0
        for i in range(N):
            suma = suma + self.b[i]*self.c[i]
        if not suma == 0.5:
            print "El metodo es de orden 1"
            return None
        if N == 2:
            print "El metodo es de orden 2"
            return None
        suma = 0
        for i in range(N):
            suma = suma + self.b[i]*(self.c[i]**2)
        if not suma == 1/3:
            print "El metodo es de orden 2"
            return None
        suma = vector(self.b)*self.A*vector(self.c)
        if not suma == 1/6:
            print "El metodo es de orden 2"
            return None
        if N == 3:
            print "El metodo es de orden 3"
            return None
        suma = 0
        suma2 = 0
        suma3 = 0
        suma4 = 0
        for i in range(N):
            suma = suma + self.b[i]*(self.c[i]**3)
            for j in range(N):
                suma2 = suma2 + self.b[i]*self.c[i]*self.A[i,j]*self.c[j]
                suma3 = suma3 + self.b[i]*self.A[i,j]*(self.c[j]**2)
                for k in range(N):
                    suma4 = suma4 + self.b[i]*self.A[i,j]*self.A[j,k]*self.c[k]
        if not suma == 0.25 or not suma2 == 0.125 or not suma3 == 1/12 or not suma4 == 1/24 :
            print "El metodo es de orden 3"
            return None
        print "El metodo es de orden 4"

    def calcular_k(self, f, x, y, h):
	k = []
	h = n(h, digits=16)
	for s in range(len(self.c)):
	    valory = 0
	    for s2 in range(s):
		valory = valory + h*n(self.A[s, s2])*k[s2]
	    xx = x + h*self.c[s]
	    yy = y + valory
	    try:
		dev = f(xx, yy).n(digits=16)
	    except:
		dev = []
		tmp = f(xx, yy)
		try:
		    tmp = len(tmp)
		except:
		    tmp = 1
		for i in range(tmp):
		    try:
		        dev.append(f(xx, yy)[i].n())
		    except:
			dev.append(f(xx,yy).n())
		dev = vector(dev)			
	    k.append(f(xx, yy).n(digits=16))
	return k
                    
    def iterar(self, f, x0, y0, h, pasos = 1):
        i = 0
        x = x0
	try:
            y = vector(y0)
	except:
	    y = y0
        resultados = [ [x,y] ]
        while i < pasos:
            k = self.calcular_k(f, x, y, h)
            resultado = 0
            for s in range(len(self.b)):
                resultado = resultado + h*self.b[s]*k[s]
	    x = x+h
	    y = y+resultado
            resultados.append([x, y])
            i = i+1
        return resultados
        
    def calcular(self, f, x0, y0, xf, N):
       h = n((xf-x0)/N)
       return self.iterar(f, x0, y0, h, N)
       
RKEuler = RungeKutta(matrix([[0]]), [1], [0])
RKModificat = RungeKutta(matrix([[0,0], [0.5, 0]]), [0,1], [0,0.5])
RKMillorat = RungeKutta(matrix([[0,0], [1,0]]), [0.5, 0.5], [0, 1])
RKHeun = RungeKutta(matrix([[0,0,0], [Rational('1/3'), 0, 0], [0, Rational('2/3'), 0]]), [0.25, 0, 0.75], [0, 1/3, 2/3])
RK3 = RungeKutta(matrix([[0,0,0], [0.5, 0,0], [-1, 2, 0]]), [Rational('1/6'), Rational('2/3'), Rational('1/6')], [0, 0.5, 1])

class RungeKuttaFehlberg(RungeKutta):
    def __init__(self, A, b1, b2, c, p):
	self.A = A
	self.b1 = b1
	self.b = b1
	self.b2 = b2
	self.c = c
	self.p = p

    def iterar(f, x0, y0, h, pasos = 1):
	print "En RKF hay que usar el metodo calcular()"

    def calcular(self, f, x0, xf, y0, h, hmax = 0.5, hmin = 1e-2, delta = 1e-2, Niter = False):
	i = 0
	x = x0
	if xf - x0 > 0:
	    signo = 1
	else:
	    signo = -1
	try:
	    y = vector(y0)
	except:
	    y = y0
        NMAX = int(abs(xf-x0)/hmin) + 10
        #print "Maximo: ", NMAX
	resultados = [ [x, y] ]
	iteraciones = 1
	while not x == xf and iteraciones < NMAX:
	    iteraciones = iteraciones + 1
	    k = self.calcular_k(f, x, y, h)
	    error = 0
	    for s in range(len(self.b1)):
		error += (self.b1[s] - self.b2[s])*h*k[s]
	    try:
		error = error.norm()
	    except:
		error = abs(error)
	    h = abs(h)
	    if not error == 0:
	        h = min(hmax, h*0.84*sqrt(delta/error, self.p), abs(xf-x))
            else:
		h = min(hmax, h, abs(xf-x))
            if not(h == abs(xf-x)):
                h = max(hmin, h)
	    h = h*signo
            #print iteraciones, h, x
	    k = self.calcular_k(f, x, y, h)
	    resultado = y
	    for s in range(len(self.b1)):
		resultado += h*self.b1[s]*k[s]
	    x = x+h
	    y = resultado
	    i += 1
	    resultados.append([x,y])
        if iteraciones == NMAX:
            return False
        if Niter:
	    print "[RKF] Han hecho falta " + str(i) + " iteraciones."
	return resultados

RKF23 = RungeKuttaFehlberg(matrix([[0,0,0], [1,0,0], [Rational('1/4'), Rational('1/4'), 0]]), [Rational('1/2'), Rational('1/2'), 0], [Rational('1/6'), Rational('1/6'), Rational('4/6')], [0, 1, Rational('1/2')], 2)
RKF45 = RungeKuttaFehlberg(matrix([ [0, 0, 0, 0, 0, 0], [Rational('1/4'), 0, 0, 0, 0, 0], [Rational('3/32'), Rational('9/32'), 0, 0, 0, 0], [Rational('1932/2197'), -Rational('7200/2197'), Rational('7296/2197'), 0, 0, 0], [Rational('439/216'), -8, Rational('3680/513'), -Rational('845/4104'), 0, 0], [-Rational('8/27'), 2, -Rational('3544/2565'), Rational('1859/4104'), -Rational('11/40'), 0]]), [Rational('25/216'), 0, Rational('1408/2565'), Rational('2197/4104'), -Rational('1/5'), 0], [Rational('16/135'), 0, Rational('6656/12825'), Rational('28561/56430'), -Rational('9/50'), Rational('2/55')], [0, Rational('1/4'), Rational('3/8'), Rational('12/13'), 1, 0.5], 4)

class Trapezi:
    @staticmethod
    def evaluar(f, x, y, h, variables):
        val1 = f.subs(variables[0] == x)
        val2 = f.subs(variables[0] == x+h)
        variables = variables[1:len(variables)]
        try:
            val1 = val1.subs(variables[0] == y)
	    y = y - variables[0]
        except:
            for j in range(len(variables)):
                val1 = val1.subs(variables[j] == y[j])
            y = y - vector(variables)
	try:
	    return (y + h*0.5*(val1() + val2()))
	except:
            return (y[0] + h*0.5*(val1() + val2()))
    
    @staticmethod
    def iterar(f, variables, x0, y0, y1, h, pasos = 1, NMAX = Infinity):
        i = 0
        x = x0
        try:
            y = vector(y0)
        except:
            y = y0
        longitud = len(variables)-1
        if longitud == 1:
	    longitud = 0
        resultados = [ [x, y] ]
        while i < pasos:
            f2 = Trapezi.evaluar(f, x, y, h, variables)
            funciones = []
            if longitud == 0:
                funciones.append(f2)
            else:
                for j in range(longitud):
                    funciones.append(f2[j])
            NN = NNewton(funciones, variables[1:len(variables)])
            y = NN.cercar_zero(y1, NMAX = NMAX)
            x = x + h
            resultados.append( [x, y] )
            i = i+1
        return resultados

    @staticmethod
    def calcular(f, variables, x0, y0, x1, y1, N, NMAX = Infinity):
	return Trapezi.iterar(f, variables, x0, y0, y1, n(x1-x0)/n(N), N, NMAX)

def taylor2(f, variables, x0, y0, x1, N):
    fx = f.diff(variables[0])
    fy = f.diff(variables[1])
    h = n(abs(x1-x0)/N)
    x = x0
    y = y0
    resultados = [ [x, y] ]
    k = 0
    while k < N:
	y = y + h*f(x,y) + ((h**2)/2.0)*(fx(x,y) + fy(x,y)*f(x,y))
        x = x+h
        resultados.append( [x, y] )
	k = k+1
    return resultados
