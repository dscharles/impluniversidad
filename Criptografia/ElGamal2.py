from sage.all import *

#y^2 = x^3 + 9x + 13, 101

def DiffieHellman(p, alpha, a, b):
    print "Queremos compartir: " + str(alpha) 
    modulo=Zmod(p)
    enviaA = modulo(a*alpha)
    print "El usuario A envia: a*P = " + str(enviaA)
    enviaB = modulo(b*P)
    print "El usuario B envia: b*P = " + str(enviaB)
    calA = modulo(a*enviaB)
    calB = modulo(b*enviaA)
    print "El usuario A hace a*Pb = " + str(calA)
    print "El usuario B hace b*Pa = " + str(calB)

def ElGamalClave(p, alpha):
    d = randint(2, p-2)
    modulo=Zmod(p)
    print "El usuario b ha cogido d = " + str(d)
    clave = modulo(d*alpha)
    print "Su clave publica es " + str(clave)
    return ((p, alpha, clave), d)

def ElGamalCod(clave, E):
    print "El usuario A quiere enviar a B " + str(E)
    modulo=Zmod(clave[0][0])
    k = randint(2, clave[0][0]-2)
    C1 =modulo(k*clave[0][1])
    C2 = modulo(E + k*clave[0][2])
    print "El usuario A envia"
    return (C1, C2)

def ElGamalDes(clave, C):
    print "El usuario B recibe C1 = " + str(C[0]) + " y C2 = " + str(C[1])
    M = C[1] - clave[1]*C[0]
    print "Ha recibido "
    return M

def firmar(clave, M):
    k = randint(2, clave[0][0]-2)
    x,y = (k*clave[0][1]).xy()
    modulo = Zmod(clave[0][0].order())
    R = modulo(x)
    S = modulo(M + R*clave[1])/modulo(k)
    return (R, S)

def verificar(clave, firma, M):
    r = clave[0][0]
    modulo = Zmod(r)
    w = modulo(1)/modulo(firma[1])
    u1 = modulo(M*w)
    u2 = modulo(firma[0]*w)
    x,y = (int(u1)*clave[0][1] + int(u2)*clave[0][2]).xy()
    v = modulo(x)
    print v, firma[0]
    return v == firma[0]

def EnviarMensaje(claves, mensaje):
	envio = []
	for letra in mensaje:
		# Toca codificar 'letra'. 
		# Pasamos el caracter a codigo ASCII
		m = ord(letra)
		print(letra + " -> " + str(m))
		envio.append(ElGamalCod(claves,m))
	return envio

def RecibirMensaje(claves, mensaje):
	recibido = []
	for letra in mensaje:
		# Primero descodificamos
		m = ElGamalDes(claves, letra)
		recibido.append("%c" % m) #Para pasar de ASCII a caracter
	final = "" # Para que el resultado sea bonito
	for i in recibido:
		final = final + i
	return final

def firmar2(clave, mensaje):
	envio = []
	for letra in mensaje:
		# Toca codificar 'letra'. 
		# Pasamos el caracter a codigo ASCII
		m = ord(letra)
		envio.append(firmar(clave,m))
	return envio

def verificar2(clave,firma,mensaje):
	i=0
	for letra in mensaje:
		# Toca codificar 'letra'. 
		# Pasamos el caracter a codigo ASCII
		m = ord(letra)
		ver=verificar(clave,firma[i],m)
		if not ver:
			return false
		i=i+1
	return true
