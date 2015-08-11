from sage.all import *

def dibujar(matriz, fila, x0, x1, M):
    vector = matriz[fila]
    h = n(abs(x1-x0))/n(M)
    puntos = []
    for i in range(len(vector)):
	puntos.append( [x0 + i*h, vector[i] ] )
    return line(puntos)

#f(x, t)
def upwindb(f, a, x0, x1, T, N, M):
    h = n(abs(x1-x0))/n(M)
    k = n(T)/n(N)
    matriz = matrix(RR, N+1, M+1)
    for i in range(N+1):
	matriz[i, 0] = 0
    for i in range(M+1):
	#print x0 + i*h,  f(x0 + i*h, 0)
	matriz[0, i] = f(x0 + i*h)
    for nn in range(1, N+1):
	for m in range(1, M+1):
	    matriz[nn, m] = matriz[nn-1, m] - a*n(k)/n(h)*(matriz[nn-1, m] - matriz[nn-1, m-1])
    return matriz

##############################################################################################

#f(x)
def upwindf(f, a, x0, x1, T, N, M):
    h = n(abs(x1-x0))/n(M)
    k = n(T)/n(N)
    matriz = matrix(RR, N+1, M+1)
    for i in range(N+1):
	matriz[i, M] = 0
    for i in range(M+1):
	#print x0 + i*h,  f(x0 + i*h, 0)
	matriz[0, i] = f(x0 + i*h)
    for nn in range(1, N+1):
	for m in range(0, M):
	    matriz[nn, m] = matriz[nn-1, m] - a*n(k)/n(h)*(matriz[nn-1, m+1] - matriz[nn-1, m])
    return matriz

def upwindfm(f, a, x0, x1, T, N, M):
    h = n(abs(x1-x0))/n(M)
    k = n(T)/n(N)
    vector = []
    vector2 = []
    for i in range(M+1):
	vector.append( f(x0 + i*h) )
        vector2.append( 0 )
    it = 0
    while it < N:
        vector[M] = 0
	for i in range(M):
            vector2[i] = vector[i] - a*k/h*(vector[i+1] - vector[i])
        it += 1
        vector = vector2
    return vector

def dibujar2(lista, x0, x1, M, col="blue"):
    h = n(abs(x1-x0))/n(M)
    puntos = []
    for i in range(M+1):
	puntos.append( [x0 + i*h, lista[i] ])
    return line(puntos, color = col)

def errores(f, a, x0, x1, T, M, N):
    err2 = []
    h = n(abs(x1-x0))/n(M)
    for i in N:
        k = n(T)/n(i)
	sol = upwindfm(f, a, x0, x1, T, i, M)
        errtmp = []
        for j in range(len(sol)):
	    errtmp.append( abs( f(x0 + j*h - a*T)  - sol[j] ) )
        err2.append( [ k, max(errtmp) ] )
    return err2

def errores2(f, a, x0, x1, T, M, N):
    err2 = []
    k = n(T)/n(N)
    for i in M:
        h = n(abs(x1-x0))/n(i)
	sol = upwindfm(f, a, x0, x1, T, N, i)
        errtmp = []
        for j in range(len(sol)):
	    errtmp.append( abs( f(x0 + j*h - a*T)  - sol[j] ) )
        err2.append( [ h, max(errtmp) ] )
    return err2
