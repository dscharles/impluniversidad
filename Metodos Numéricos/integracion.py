from sage.all import *

def simpsons_rule(f,a,b):
    """
    Devuelve la aproximacion de la integral de f en [a,b] mediante la regla de simpsons

    EJEMPLO:
    >>> f(x) = x
    >>> simpsons_rule(f,0,1)
    0.5"""
    c = (a+b)/2.0
    h3 = abs(b-a)/6.0
    return n(h3*(f(a) + 4.0*f(c) + f(b)))

def simpsons_rule_M(f, a, b, N):
	"""
	Devuelve la aproximacion de la integral de f en [a,b] dividiendo el intervalo en N trozos y aplicando la regla de simpsons en cada intervalo"""
	h = (b-a)/N
	valor = 0
	for i in range(N):
		valor = valor + simpsons_rule(f,a + i*h, a + (i+1)*h)
	return valor

def adaptative_simpsons_rule(f,a,b,eps):
    """
    Devuelve la aproximacion de la integral de f en [a,b] con un error mas pequenyo que eps aplicando la regla de simpsons adaptativa"""
    return recursive_asr(f,a,b,eps, simpsons_rule(f,a,b))

def recursive_asr(f,a,b,eps,sum):
    """
    Usar adaptative_simpsons_rule()"""
    c = (a+b)/2.0
    left = simpsons_rule(f,a,c)
    right = simpsons_rule(f,c,b)
    if abs(left + right - sum) <= 15*eps :
        return left + right + (left + right - sum)/15.0
    return recursive_asr(f,a,c,eps/2.0,left) + recursive_asr(f,c,b,eps/2.0,right)
