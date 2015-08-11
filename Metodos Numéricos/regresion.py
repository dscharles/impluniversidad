from sage.all import *

def orden(lista, minimo, p):
    ord = []
    for i in range(1, len(lista)):
        ord.append(abs(lista[i] - minimo)/abs(lista[i-1] - minimo)^p)
    return ord
    
def puntitos(lista, minimo):
    ord = []
    for i in range(1, len(lista)):
	if not (lista[i] - minimo == 0 or lista[i-1] - minimo == 0):
            ord.append( [log(abs(lista[i-1] - minimo), e),log(abs(lista[i] - minimo), e)])
    return ord

def valores_recta_regresion(puntos):
    sumax = 0
    sumay = 0
    sumax2 = 0
    sumay2 = 0
    sumaxy = 0
    for punto in puntos:
        sumax = sumax + punto[0]
        sumax2 = sumax2 + punto[0]**2
        sumay = sumay + punto[1]
        sumay2 = sumay2 + punto[1]**2
        sumaxy = sumaxy + punto[0]*punto[1]
    mediax = sumax/len(puntos)
    mediay = sumay/len(puntos)
    varx = sumax2/len(puntos) - mediax**2
    varxy = sumaxy/len(puntos) - mediax*mediay
    return [mediax, mediay, varxy/varx]

def recta_regresion(valores, minimo):
    puntos = puntitos(valores, minimo)
    x = var('x')
    p = point(puntos[0])
    mi = 0
    ma = 0
    for punto in puntos:
	if mi > punto[0]:
	    mi = punto[0]
	if ma < punto[0]:
	    ma = punto[0]
        p = p + point(punto)
    mx, my, m = valores_recta_regresion(puntos)
    p = p + plot(my+m*(x-mx),mi, ma, color = "red")
    return [p, m, exp(my-m*mx)]
