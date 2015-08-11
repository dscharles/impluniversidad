from sage.all import *

#y^2 = x^3 + 9x + 13, 101
def CurvaEliptica(poli, n):
    return EllipticCurve( poli ).change_ring(Zmod(n))

def DiffieHellman(CE, P, a, b):
    print "Queremos compartir: " + str(P) + " (orden: " + str(P.order()) + ")"
    enviaA = a*P
    print "El usuario A envia: a*P = " + str(enviaA)
    enviaB = b*P
    print "El usuario B envia: b*P = " + str(enviaB)
    calA = a*enviaB
    calB = b*enviaA
    print "El usuario A hace a*Pb = " + str(calA)
    print "El usuario B hace b*Pa = " + str(calB)

def ElGamalClave(CE, P):
    r = CE.order()
    d = randint(2, r-2)
    print "El usuario b ha cogido d = " + str(d)
    clave = d*P
    print "Su clave publica es " + str(clave)
    return ((CE, P, clave), d)

def ElGamalCod(clave, E):
    print "El usuario A quiere enviar a B el punto " + str(E)
    k = randint(2, clave[0][0].order()-2)
    C1 = k*clave[0][1]
    C2 = E + k*clave[0][2]
    print "El usuario A envia"
    return (C1, C2)

def ElGamalDes(clave, C):
    print "El usuario B recibe C1 = " + str(C[0]) + " y C2 = " + str(C[1])
    M = C[1] - clave[1]*C[0]
    print "Ha recibido "
    return M

def firmar(clave, M):
    k = randint(2, clave[0][0].order()-2)
    x,y = (k*clave[0][1]).xy()
    modulo = Zmod(clave[0][0].order())
    R = modulo(x)
    S = modulo(M + R*clave[1])/modulo(k)
    return (R, S)

def verificar(clave, firma, M):
    r = clave[0][0].order()
    modulo = Zmod(r)
    w = modulo(1)/modulo(firma[1])
    u1 = modulo(M*w)
    u2 = modulo(firma[0]*w)
    x,y = (int(u1)*clave[0][1] + int(u2)*clave[0][2]).xy()
    v = modulo(x)
    print v, firma[0]
    return v == firma[0]
