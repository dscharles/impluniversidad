from sage.all import *

class Lectura:
    def __init__(self, texto):
        self.texto = texto
        self.text2 = texto
        self.i = 0
        
    def read(self, num):
        st = self.text2[0:num]
        self.text2 = self.text2[num: len(self.text2)]
        return st
        
    def seek(self, i):
        self.text2 = self.texto

    def close(self):
	self.i = 0

class Elemento:
    def __init__(self, caracter, k):
        self.caracter = caracter
        self.probabilidad = 0
        self.k = k
        
    def __cmp__(self, x):
        if self.k < x.k:
            return -1
        if self.k > x.k:
            return 1
        return 0
        
    def otro_mas(self, valor = 1):
        self.k += valor
        
    def getCaracter(self):
        return self.caracter
        
    def getFrecuencia(self):
        return self.k
        
    def getProbabilidad(self):
        return self.probabilidad
        
    def __repr__(self):
        return str( (self.caracter, self.k, self.probabilidad) )


def abrir(nombre):
    return codecs.open(nombre, "r", "utf-8")

def array_to_str(lista):
    ret = ""
    for p in lista:
	ret = ret + " " + p
    return ret[1:len(ret)]

def incluir(diccionario, palabra):
    if palabra in diccionario:
	diccionario[palabra].otro_mas()
    else:
	diccionario[palabra] = Elemento(palabra, 1)

def leer2(texto):
    texto.seek(0)
    texto1 = {}
    texto2 = {}
    texto3 = {}
    texto4 = {}
    textoword = {}
    textoword2 = {}
    textoword3 = {}
    tmp2 = []
    tmp3 = []
    tmp4 = []
    tmpword = ""
    tmpword2 = []
    tmpword3 = []
    car = texto.read(1)
    while not car == "":
	if car == " " and not tmpword == "":
           incluir(textoword, tmpword)
	   tmpword2.append(tmpword)
	   if len(tmpword2) == 2:
		incluir(textoword2, array_to_str(tmpword2))
	   	tmpword2.remove(tmpword2[0])
	   tmpword3.append(tmpword)
	   if len(tmpword3) == 3:
		incluir(textoword3, array_to_str(tmpword3))
		tmpword3.remove(tmpword3[0])
	   tmpword = ""
	if not car == " ":
	   tmpword += car
	incluir(texto1, car)
	for i in range(len(tmp2)):
	    tmp2[i] += car
	tmp2.append(car)
	if len(tmp2[0]) == 2:
	    incluir(texto2, tmp2[0])
	    tmp2.remove(tmp2[0])
	for i in range(len(tmp3)):
	    tmp3[i] += car
	tmp3.append(car)
	if len(tmp3[0]) == 3:
	    incluir(texto3, tmp3[0])
	    tmp3.remove(tmp3[0])
	for i in range(len(tmp4)):
	    tmp4[i] += car
	tmp4.append(car)
	if len(tmp4[0]) == 3:
	    incluir(texto4, tmp4[0])
	    tmp4.remove(tmp4[0])
	car = texto.read(1)
    if not tmpword == "":
	incluir(textoword, tmpword)
	tmpword2.append(tmpword)
	if len(tmpword2) == 2:
	    incluir(textoword2, array_to_str(tmpword2))
	tmpword3.append(tmpword)
	if len(tmpword3) == 3:
	    incluir(textoword3, array_to_str(tmpword3))
    return [texto1, texto2, texto3, texto4, textoword, textoword2, textoword3]
	   
class Markov:
    def __init__(self, texto):
	self.texto = texto
	self.dic = {}

    def crear_aleatorio(self, palabra, k):
	#k = 1 --> 1
	#k = 2 --> 2
	#k = 3 --> 3
	#k = "word" --> 5
	#k = "word2" --> 6
        if k == "word":
	    k = 5
	if k == "word2":
	    k = 6
	listado = []
	probabilidades = []
	suma = 0
	for p in self.texto[k]:
	    if k == 5 or k == 6:
		if len(p) > len(palabra):
		    if p[0:len(palabra)] == palabra and p[len(palabra)] == " ":
			listado.append(p)
			probabilidades.append(self.texto[k][p].getFrecuencia())
			suma += self.texto[k][p].getFrecuencia()
	    else:
		if p[0:len(palabra)] == palabra:
		    listado.append(p)
		    probabilidades.append(self.texto[k][p].getFrecuencia())
	            suma += self.texto[k][p].getFrecuencia()
        for i in range(len(probabilidades)):
	    probabilidades[i] = n(probabilidades[i]/suma)
        if len(probabilidades) == 0:
	    return False
	else:
	    self.dic[palabra] = [GeneralDiscreteDistribution(probabilidades), listado]

    def sortear2(self, palabra, k):
	if not palabra in self.dic:
	    if self.crear_aleatorio(palabra, k) == False:
		return False
        return self.dic[palabra][1][self.dic[palabra][0].get_random_element()]

    def sortear(self, palabra, k):
	res = self.sortear2(palabra, k)
	if not res:
	    return False
	return res[len(palabra):len(res)]

    def sorteo_inicial(self, k):
	palabras = []
	prob = []
	suma = 0
	if k == "word":
	    k = 5
	if k == "word2":
	    k = 6
	for p in self.texto[k]:
	    palabras.append(p)
	    prob.append(self.texto[k][p].getFrecuencia())
	    suma += self.texto[k][p].getFrecuencia()
	for i in range(len(prob)):
	    prob[i] = n(prob[i]/suma)
	if len(prob) == 0:
	    return False
	X = GeneralDiscreteDistribution(prob)
	return palabras[X.get_random_element()]

    def generar_texto(self, longitud, k):
	#self.dic = {}
	anterior = self.sorteo_inicial(k)
	if anterior == False:
	    print "No hay palabras suficientes"
	    return 
	texto = anterior
	if k == "word" or k == "word2":
	    anterior = anterior.split(" ")[1]
	else:
	    anterior = anterior[1:len(anterior)]
	for i in range(longitud):
	    anterior = self.sortear(anterior, k)
	    if anterior == False:
		print "No hay palabras suficientes"
		return texto
	    texto += anterior
	    if k == "word" or k == "word2":
		anterior = anterior[1:len(anterior)]
	return texto

class Markov2:
    def __init__(self, texto):
        self.texto = texto
        self.dic = {}
        
    def generar_texto_orden1(self, longitud):
        prob = []
        palabras = []
        suma_total = 0
        for i in self.texto[0]:
            suma_total += self.texto[0][i].getFrecuencia()
            palabras.append(i)
        for i in self.texto[0]:
            prob.append(n(self.texto[0][i].getFrecuencia()/suma_total))
        X = GeneralDiscreteDistribution(prob)
        st = ""
        for i in range(longitud):
            st += palabras[X.get_random_element()]
        return st
        
    def generar_texto_orden0(self, longitud):
        prob = []
        palabras = []
        suma_total = len(self.texto[0])
        for i in self.texto[0]:
            palabras.append(i)
            prob.append(n(1/suma_total))
        X = GeneralDiscreteDistribution(prob)
        st = ""
        for i in range(longitud):
            st += palabras[X.get_random_element()]
        return st
        
    def generar_texto_ordenword1(self, longitud):
        prob = []
        palabras = []
        suma_total = 0
        for i in self.texto[4]:
            suma_total += self.texto[4][i].getFrecuencia()
            palabras.append(i)
        for i in self.texto[4]:
            prob.append(n(self.texto[4][i].getFrecuencia()/suma_total))
        X = GeneralDiscreteDistribution(prob)
        st = ""
        for i in range(longitud):
            st += palabras[X.get_random_element()] + " "
        return st
