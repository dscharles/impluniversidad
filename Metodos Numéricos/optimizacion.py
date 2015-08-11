"""
Prueba ayuda
"""

from sage.all import *
from Metodos.interpolacion import *
from utiles import *

def forward_backward(f, alpha0, h0, t = 2, NMAX = Infinity, Niter=False):
    """
    Busca un intervalo inicial donde la funcion f tenga un posible minimo.

    Recibe como parametros un valor inicial alpha, y un salto inicial h.
    t indica cuanto se amplia el paso en cada iteracion.
    """
    alpha = [alpha0]
    h = [h0]
    phi = [f(alpha[0])]
    k = 0
    i = 0
    while i < NMAX:
        i = i+1
        try:
            alpha[k+1] = alpha[k]+h[k]
            phi[k+1] = f(alpha[k+1])
        except:
            alpha.append(alpha[k]+h[k])
            phi.append(f(alpha[k+1]))
        if(phi[k+1] < phi[k]):
            try:
                h[k+1] = t*h[k]
            except:
                h.append(t*h[k])
            k = k+1
        else:
            if k == 0:
                h[k] = -h[k]
                alpha[k] = alpha[k+1]
                phi[k] = phi[k+1]
            else:
                a = min(alpha[k-1], alpha[k+1])
                b = max(alpha[k-1], alpha[k+1])
		if Niter:
		   print("Se ha resuelto en " +str(i) + " iteraciones")
                return([a,b])
    print("[Forward-backward] Superado maximo iteraciones")
    a = min(alpha[k-1], alpha[k])
    b = max(alpha[k-1], alpha[k])
    return ([a,b])

def golden_section(f, a0, b0, delta = 1e-6, NMAX = Infinity, digitos = 16,Niter=False):
    """
    Aplica el metodo de la seccion aurea a la funcion f en el intervalo [a, b]
    """
    a0 = obtener_precision(a0, digitos)
    b0 = obtener_precision(b0, digitos)
    razon = obtener_precision((sqrt(5)-1)/2, digitos)
    a = [a0]
    b = [b0]
    lam = []
    mu = []
    k = 0
    while k < NMAX:
        if len(lam) == k:
            lam.append(a[k] + (1-razon)*(b[k]-a[k]))
        if len(mu) == k:
            mu.append(a[k] + razon*(b[k]-a[k]))
        if(f(lam[k]) < f(mu[k])):
            if(mu[k] - a[k] < delta):
		if Niter:
		   print("Se ha resuelto en " + str(k) + " iteraciones")
                return lam[k]
            else:
                a.append(a[k])
                b.append(mu[k])
                mu.append(lam[k])
        elif(f(lam[k]) == f(mu[k])):
            print("Puede no ser preciso")
	    if Niter:
		   print("Se ha resuelto en " + str(k) + " iteraciones")
            return (lam[k]+mu[k])/2
        else:
            if(b[k] - lam[k] < delta):
		if Niter:
		   print("Se ha resuelto en " + str(k) + " iteraciones")
                return mu[k]
            else:
                a.append(lam[k])
                b.append(b[k])
                lam.append(mu[k])
        k = k+1
    print("[Golden section] Superado maximo iteraciones")
    return (lam[k]+mu[k])/2

def obtener_n_fibonacci(a, b, epsilon):
    n = 1
    while abs(a-b)/epsilon > fibonacci(n):
        n = n+1
    return n

def metodo_fibonacci(f, a0, b0, delta = 1e-6, NMAX = Infinity, digitos = 16, Niter=False):
    """
    Aplica el metodo de fibonacci a la funcion f en el intervalo [a,b]
    """
    a = [obtener_precision(a0, digitos)]
    b = [obtener_precision(b0, digitos)]
    lam = []
    mu = []
    k = 0
    n = obtener_n_fibonacci(a[k], b[k], delta)
    sol = []
    while k < NMAX:
        sol.append((a[k]+b[k])/2)
        if len(lam) == k:
            lam.append(a[k] + obtener_precision(fibonacci(n-k-1)/fibonacci(n-k+1), digitos)*(b[k]-a[k]))
        if len(mu) == k:
            mu.append(a[k] + obtener_precision(fibonacci(n-k)/fibonacci(n-k+1), digitos)*(b[k]-a[k]))
        if(f(lam[k]) - f(mu[k]) <0):
            if(mu[k]-a[k] < delta):
		if Niter:
		   print("Se ha resuelto en " + str(k) + " iteraciones")
                return (lam[k], sol)
            else:
                a.append(a[k])
                b.append(mu[k])
                mu.append(lam[k])
        else:
            if(f(lam[k]) - f(mu[k]) ==0):
		print("Puede haber error de precision")
		if Niter:
		   print("Se ha resuelto en " + str(k) + " iteraciones")
                return (((lam[k]+mu[k])/2), sol)
            if(b[k] - lam[k] < delta):
		if Niter:
		   print("Se ha resuelto en " +str(k)+" iteraciones")
                return (mu[k], sol)
            else:
                a.append(lam[k])
                b.append(b[k])
                lam.append(mu[k])
        k = k+1
    print("[Fibonacci] Maximo iteraciones alcanzado")
    return ((lam[k]+mu[k])/2, sol)

def metodo_interpolador(f, x1, x2, x3, delta = 1e-6, digitos = 16, NMAX = Infinity,Niter=False):
    k = 0
    x1 = obtener_precision(x1, digitos)
    x2 = obtener_precision(x2, digitos)
    x3 = obtener_precision(x3, digitos)
    if f(x1) < f(x2) or f(x2) > f(x3):
        print("[Interpolador] Falla condicion inicial")
        return None
    while k < NMAX:
        if abs(x1-x3)<delta:
	    if Niter:
		   print("Se ha resuelto en "+str(k)+" iteraciones")
            return x2
        P = PolinomioLagrange([[x1, f(x1)], [x2, f(x2)], [x3, f(x3)] ])
        P.generar_polinomio()
        if(len(P.poli.coeffs())==3):
            v = -obtener_precision(P.poli.coeffs()[1]/(2*P.poli.coeffs()[2]), digitos)
        else:
	    print("Se ha dado el caso en que f(x1) = f(x3)")
            return x1
        if v < x2:
            if f(x2) >= f(v):
                # x1 = x1
		x3 = x2
		x2 = v
	    else:
		x1 = v
		# x2 = x2
		# x3 = x3
	elif v > x2:
    	    if f(v) <= f(x2):
		x1 = x2
		x2 = v
		# x3 = x3
	    else:
		# x1 = x1; x2 = x2
		x3 = v
	else:
            print("[Interpolador] El vertice cae en x2, cambiar condicion inicial.")
	    return None
	k = k+1
    print("[Metodo interpolador] Maximo iteraciones superado")
    return x3

def metodo_interpolador_hermite(f, df, x1, h, delta=1e-6, NMAX = Infinity, Niter = False, digitos=16):
    print "##", f(0), f(0.0005), f(-0.0005)
    print "###", df(0)
    h = n(h, digits = digitos)
    k = 0
    i = 0
    mu = obtener_precision(x1, digitos)
    fi = f(mu)
    fip = df(mu)
    sol = []
    while k < NMAX:
	i = i+1
        sol.append(mu)
 	print "MINT", i, mu, h
        if abs(fip) < delta:
            if Niter:
                print("Ha habido " + str(k) + " iteraciones.")
            return mu
        if fip > 0:
            h = -abs(h)
        else:
            h = abs(h)
        mu2 = mu + h
        fi2 = f(mu2)
        fip2 = df(mu2)
        while fip2*fip >= 0:
	    print "@", mu, fip, mu2, fip2, h
            h = 2*h
            mu = mu2
            fi = fi2
            fip = fip2
            mu2 = mu + h
            fi2 = f(mu2)
            fip2 = df(mu2)
        if True:
            if h > 0:
                a = mu
                b = mu2
            else:
                a = mu2
                b = mu
                tfi = fi
                fi = fi2
                fi2 = tfi
                tfi = fip
                fip = fip2
                fip2 = tfi
	    if abs(b-a)<delta:
		print("petaaaa")
		return mu
	    s = 3*(f(mu2)-f(mu))/(b-a)
            z = s - df(mu2) - df(mu)
            w = sqrt(z**2 - df(mu)*df(mu2))
            q = 1 - (df(mu2)+w+z)/(df(mu2)-df(mu)+2*w)
            mu = mu + (b-a)*q
            k = k+1
            if abs(df(mu)) < delta:
                if Niter:
                    print("Ha habido " + str(k) + " iteraciones.")
                return mu
            else:
                h = n(h/10)
    print("[Interpolador 2] Fuera iteraciones")
    return False, mu

def metodo_armijo(f, df, a, r, NMAX = Infinity,Niter=False):
    """
    Aplica el metodo de armijo para obtener un valor de alpha > 0 donde la funcion ha decrecido
    """
    k = 0
    while k < NMAX:
        k = k+1
        if f(a) < f(0):
            if f(a) <= f(0) + r*a*df(0):
                if Niter:
		   print("Se ha resuelto en "+str(k)+" iteraciones")
                return n(a)
            else:
                a = a/2
        else:
            a = a/2
    print("[Armijo] Maximo iteraciones superado")
    return a
def metodo_goldstein(f, df, a, r, NMAX = Infinity, Niter=False):
    """
    Aplica el metodo de golstein para obtener un valor de alpha > 0 donde la funcion ha decrecido
    """
    t=2
    ap=0
    ag=0
    phi0=f(0)
    phi10=df(0)
    phi=f(a)
    k = 0
    while k < NMAX:
        if phi < phi0+r*a*phi10:
            if phi >= phi0 + (1-r)*a*phi10:
		return n(a)
	    else:
		ap=a
	        if ag != 0:
		    a=(ap+ag)/2
		else:
		    a=t*a
        else:
	     ag=a
	     a=(ap+ag)/2
        k = k+1
    print("[Goldstein] Maximo iteraciones superado")
    return a
   
def metodo_goldstein2(f, df, a, r, NMAX = Infinity, Niter=False):
    """
    Aplica el metodo de golstein para obtener un valor de alpha > 0 donde la funcion ha decrecido
    """
    k = 0
    while k < NMAX:
        k = k+1
        if f(a) < f(0):
            if f(a) <= f(0) + r*a*df(0):
                if f(a) >= f(0) + (1-r)*a*df(0):
		    if Niter:
		   	print("Se ha resuelto en "+str(k)+" iteraciones")
                    return n(a)
                else:
                    a = n(a)*3/2
            else:
                a = n(a)/2
        else:
            a = n(a)/2
    print("[Goldstein] Maximo iteraciones superado")
    return a

def metodo_wolfe_powell(f, df, a, r, ro, NMAX = Infinity,Niter=False):
    """
    Aplica el metodo de Wolfe-Powell para obtener un valor de alpha > 0 donde la funcion ha crecido
    """
    k = 0
    while k < NMAX:
        k = k+1
        if f(a) < f(0):
            if f(a) <= f(0) + r*a*df(0):
                if df(a) >= ro*df(0):
		    if Niter:
		   	print("Se ha resuelto en "+str(k)+" iteraciones")
                    return n(a)
                else:
                    a = a*3/2
            else:
                a = a/2
        else:
            a = a/2
    print("[Wolfe-Powell] Maximo iteraciones superado")
    return a

class Funcion:
    def __init__(self, f, Df = None, digitos=16):
        self.funcion = f
        if Df == None:
            self.derivada = self.funcion.diff()
        else:
            self.derivada = Df
        self.digitos = digitos
        
    def __call__(self, punto):
        return self.funcion(punto)
        
    def __float__(self, punto):
        return self(punto)
                
    def forward_backward(self, alpha0, h0, t=2, NMAX = Infinity,Niter=False):
        return forward_backward(self, alpha0, h0, t, NMAX,Niter)
        
    def golden_section(self, a0, b0, delta = 1e-6, NMAX = Infinity,Niter=False):
        return golden_section(self, a0, b0, delta, NMAX, self.digitos, Niter)
        
    def metodo_fibonacci(self, a0, b0, delta = 1e-6, NMAX = Infinity, Niter=False):
        return metodo_fibonacci(self, a0, b0, delta, NMAX, self.digitos,Niter)
        
    def metodo_interpolador(self, x1, x2, x3, delta = 1e-6, NMAX = Infinity,Niter=False):
        return metodo_interpolador(self, x1, x2, x3, delta,NMAX, self.digitos,Niter)

    def metodo_interpolador_hermite(self, x1, h, delta=1e-6, NMAX = Infinity, Niter = False, digitos = 16):
        return metodo_interpolador_hermite(self, self.derivada, x1, h, delta, NMAX, Niter, digitos)

        
    def metodo_armijo(self, a = 1, r = 0.3, NMAX = Infinity,Niter=False):
        return metodo_armijo(self, self.derivada, a, r, NMAX,Niter)
        
    def metodo_goldstein(self, a = 1, r = 0.3, NMAX = Infinity,Niter=False):
        return metodo_goldstein(self, self.derivada, a, r, NMAX,Niter)
        
    def metodo_wolfe_powell(self, a = 1, r = 0.3, ro = 0.8, NMAX = Infinity,Niter=False):
        return metodo_wolfe_powell(self, self.derivada, a, r, ro, NMAX,Niter)

class NFuncion:
    def __init__(self, f, Df = None, digitos = 16):
        self.funcion = f
        self.variables = f.arguments()
        self.digitos = digitos
	self.funcionesdireccion = []
	self.funcionesderivada = []
        if Df != None:
            self.derivada = Df
        else:
            try:
                self.derivada = f.diff()
            except:
                self.derivada = None
             
    def __call__(self, punto):
        f = self.funcion
        for i in range(len(self.variables)):
            f = f.subs(self.variables[i] == punto[i])
        return f()
        
    def diff(self):
        return self.derivada
        
    def gradiente(self, punto):
        grad = self.derivada
        for i in range(len(self.variables)):
            grad = grad.subs(self.variables[i] == punto[i])
        return vector(grad())
        
    def funcion_direccion(self, punto, direccion):
        return lambda e : self(vector(punto) + e*vector(direccion))
        
    def funcion_derivada(self, punto, direccion):
        return lambda e : self.gradiente(vector(punto) + e*vector(direccion))*vector(direccion)
        
    def _limpiar_(self, punto):
        for i in range(len(punto)):
            punto[i] = obtener_precision(punto[i], self.digitos)
        return punto
        
    def metodo_gradiente_con_reinicio(self, punto, a=1, r=0.3, delta = 1e-6, NMAX = Infinity, Niter = False, cuadratica = False):
        punto = self._limpiar_(punto)
        k = 0
        i = 0
        aprox = vector(punto)
	aproximaciones = [aprox]
	todas = []
        while k < NMAX:
	    print k, i, n(self.gradiente(aprox).norm()), aproximaciones[k]
            if i == 0:
                direcciones = [-self.gradiente(aprox)]
            if self.gradiente(aprox).norm() < delta:
		todas.append(direcciones)
                if Niter:
                    print "Han hecho falta " + str(k+1) + " iteraciones."
                return aprox, todas, aproximaciones
	    print vector(self.gradiente(aprox))*direcciones[i] <= 0
            F = Funcion(self.funcion_direccion(aprox, direcciones[i]), self.funcion_derivada(aprox, direcciones[i]))
	    self.funcionesdireccion.append(self.funcion_direccion(aprox, direcciones[i]))
	    self.funcionesderivada.append(self.funcion_derivada(aprox, direcciones[i]))
	    if not cuadratica:
#                alpha = n(F.metodo_goldstein(a = a, r = r, NMAX = Infinity, Niter = False))
		#alpha=minimize(self.funcion_direccion(aprox,direcciones[i]),[0],disp=0)[0]
		alpha = n(F.metodo_interpolador_hermite(x1=a,h=2,NMAX=NMAX))
	    else:
		alpha = self.gradiente(aprox).norm()**2/(direcciones[i]*self.funcion.diff(2)()*direcciones[i])
            aproxanterior = aprox
            aprox = aprox + alpha*direcciones[i]
	    aproximaciones.append(aprox)
            direcciones.append(-self.gradiente(aprox))
            i = i+1
            k = k+1
            if self.gradiente(aprox).norm() < delta or (aproximaciones[k-1]-aproximaciones[k]).norm()<delta :
		todas.append(direcciones)
                if Niter:
                    print "Han hecho falta " + str(k+1) + " iteraciones."
                return aprox, todas, aproximaciones
            if i == len(punto):
                i = 0
		todas.append(direcciones)
                direcciones = None
            else:
                beta = (self.gradiente(aprox)*self.gradiente(aprox))/(self.gradiente(aproxanterior)*self.gradiente(aproxanterior))
                direcciones[i] = direcciones[i] + beta*direcciones[i-1]
                if direcciones[i]*self.gradiente(aprox) > 0:
                    i = 0
		    todas.append(direcciones)
                    direcciones = None
        print "[Metodo gradiente con reinicio] Me he pasado de las iteraciones permitidas"
        return aprox, todas, aproximaciones
        
    def _matriz_precondicionado_(self, punto):
        M = matrix(RR, len(punto), len(punto))
        for i in range(len(punto)):
            f = self.funcion.diff(self.variables[i]).diff(self.variables[i])
            for j in range(len(punto)):
                f = f.subs(self.variables[j] == punto[j])
            if f() == 0:
                M[i,i] = 1
            else:
                M[i,i] = n(1/f)
        return M
            

    def metodo_gradiente_precondicionado(self, punto, a = 1, r = 0.3, delta = 1e-6, NMAX = Infinity, Niter = False, cuadratica = False):
        punto = self._limpiar_(punto)
        aprox = vector(punto)
        k = 0
        i = 0
        d = []
        while k < NMAX:
	    print k, i, self.gradiente(aprox).norm()
            if i == 0:
                g = [self.gradiente(aprox)]
                v = [vector(self._matriz_precondicionado_(aprox)*g[0])]
                d = [-v[0]]
            if g[i].norm() < delta :
                if Niter == True:
                    print "Han hecho falta " + str(k+1) + " iteraciones."
                return aprox
            F = Funcion(self.funcion_direccion(aprox,d[i]), self.funcion_derivada(aprox, d[i]))
            if not cuadratica:
#                alpha = F.metodo_goldstein(a = a, r = r, NMAX = Infinity, Niter = False)
		alpha = n(F.metodo_interpolador_hermite(x1=a,h=2,NMAX=NMAX))
	    else:
		alpha = g[i]*v[i]/(d[i]*self.funcion.diff(2)()*d[i])
            aproxanterior = aprox
            aprox = aprox + alpha*d[i]
            k = k+1
            i = i+1
            if i == len(punto):
                i = 0
                g = None
                v = None
                d = None
            else:
                g.append(self.gradiente(aprox))
                v.append(vector(self._matriz_precondicionado_(aprox)*g[i]))
                if g[i].norm() < delta:
                    if Niter == True:
                        print "Han hecho falta " + str(k+1) + " iteraciones."
                    return aprox
                beta = (g[i]*v[i])/(g[i-1]*v[i-1])
                d.append(-v[i]+beta*d[i-1])
                if d[i]*g[i] > 0:
                    i = 0
                    g = None
                    v = None
                    d = None
        print "[Metodo gradiente precondicionado] Maximo iteraciones superado."
        return aprox
