from sage.all import *

def ElGamalClave(p, alpha, r = None):
    MOD = Zmod(p)
    if r == None:
        r = MOD(randint(2, p-2))
    y = MOD(alpha)**r
    return [ [MOD, y, MOD(alpha), p], r ]

def ElGamalCod(clave, m, coprimo = False):
    k = randint(2, clave[0][3]-2)
    if coprimo:
        while gcd(k, int(clave[0][3]) - 1) != 1:
            k = randint(2, clave[0][3]-2)
    K = clave[0][2]**k
    c = clave[0][0](m)*clave[0][1]**k
    return ( (K, c), k )
def ElGamalSend(clave,mensaje,coprimo=False):
    envio = []
    for letra in mensaje:
	# Toca codificar 'letra'. 
	# Pasamos el caracter a codigo ASCII
	m = ord(letra)
	print(letra + " -> " + str(m))
	envio.append(ElGamalCod(clave,m))
    return envio

def ElGamalDes(clave, R):
    b = R[0][0]**clave[1]
    return clave[0][0](R[0][1])/clave[0][0](b)

def ElGamalRec(clave,mensaje):
    recibido = []
    for letra in mensaje:
	# Primero descodificamos
	m = ElGamalDes(clave, letra)
	recibido.append("%c" % m) #Para pasar de ASCII a caracter
    final = "" # Para que el resultado sea bonito
    for i in recibido:
	final = final + i
    return final


# EnviarFirmado(Clave de quien envia, Clave de quien recibe, mensaje)
def EnviarFirmado(claves, clave, m):
    cod = ElGamalCod(clave, m, coprimo = True)
    MOD = Zmod(claves[0][3]-1)
    s = MOD(m - MOD(claves[1])*MOD(cod[0][0]))*MOD(cod[1])**(-1)
    return ( (cod[0][0], cod[0][1], s), cod[1])
    
def SendFirm(claves,clave,mensaje):
    envio = []
    for letra in mensaje:
	# Toca codificar 'letra'. 
	# Pasamos el caracter a codigo ASCII
	m = ord(letra)
	print(letra + " -> " + str(m))
	envio.append(EnviarFirmado(claves,clave,m))
    return envio

# ComprobarFirma(Clave de quien ha enviado, Clave de quien ha recibido, mensaje)
def ComprobarFirma(claves, clave, R):
    m = ElGamalDes(clave, R)
    print "El mensaje enviado es ", m
    MOD = Zmod(claves[0][3]-1)
    return [claves[0][2]**m == (claves[0][1]**R[0][0])*(R[0][0]**R[0][2]), m]

def checkFirm(claves,clave,mensaje):
    envio = True
    for letra in mensaje:
	# Toca codificar 'letra'. 
	# Pasamos el caracter a codigo ASCII
	envio=ComprobarFirma(claves,clave,letra)
	print envio
    return envio

