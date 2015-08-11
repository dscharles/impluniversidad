from sage.all import *
from Tema1 import *

# Devuelve la lista
# [p, q, n, fi(n), e, d]
# Estos valores son (x = privados, o = publicos)
# [x, x, o, x,     o, x]
def Claves(p, q, e):
	lista = []
	lista.append(p)
	lista.append(q)
	lista.append(p*q)
	fi = (p-1)*(q-1)
	lista.append(fi)
	lista.append(e)
	lista.append(Inverso(int(e), int(fi)))
	return lista
	
def Cifrar(m, A):
	return ExpRapida(m, A[4], A[2])

def Descifrar(m, A):
	return ExpRapida(m, A[5], A[2])

def Enviar(mensaje, A, B):
	s = Descifrar(mensaje, A)
	c = Cifrar(s, B)
	print(str(mensaje) + " --sA--> " + str(s) + " --cB--> " + str(c))
	return c

def Recibir(mensaje, A, B):
	c = Descifrar(mensaje, B)
	s = Cifrar(c, A)
	print(str(mensaje) + " --dB--> " + str(c) + " --sA--> " + str(s))
	return s

def EnviarMensaje(mensaje, A):
	envio = []
	for letra in mensaje:
		# Toca codificar 'letra'. 
		# Pasamos el caracter a codigo ASCII
		m = ord(letra)
		#print(letra + " -> " + str(m))
		envio.append(Cifrar(m, A))
	return envio

def RecibirMensaje(mensaje, A):
	recibido = []
	for letra in mensaje:
		# Primero descodificamos
		m = Descifrar(letra, A)
		recibido.append("%c" % m) #Para pasar de ASCII a caracter
	final = "" # Para que el resultado sea bonito
	for i in recibido:
		final = final + i
	return final

def EnviarMensaje2(mensaje, A, B):
	envio = []
	for letra in mensaje:
		# Pasamos la letra a un numero
		m = ord(letra)		
		#Firmamos la letra
		s = Descifrar(m, A)
		# ciframos con la clave de B
		c = Cifrar(s, B)
		print(letra + " -> " + str(m) + " -> " + str(s) + " -> " + str(c))
		envio.append(c)
	return envio

def RecibirMensaje2(mensaje, A, B):
	recibido = []
	for letra in mensaje:
		s = Descifrar(letra, B)
		m = Cifrar(s, A)
		recibido.append("%c" % m)
		l = recibido[-1]
		print(str(letra) + " -> " + str(s) + " -> " + str(m) + " -> " + str(l))
	final = ""
	for i in recibido:
		final = final + i
	return final
