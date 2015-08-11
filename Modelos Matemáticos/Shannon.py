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
        
def leer(texto, k = 1):
    valor = []
    if not k == "word":
        st = texto.read(k)
        while not st == "":
            valor.append(st)
            st = texto.read(k)
    else:
        st = " "
        while not st == "":
            st = ""
            caracter = "a"
            while not(caracter == '' or caracter == ' '):
                caracter = texto.read(1)
                st = st + caracter
            if not st == "":
                if st == " ":
                    valor.append(" ")
                else:
                    if caracter == ' ':
                        valor.append(st[0:len(st)-1])
                        valor.append(' ')
                    else:
                        valor.append(st)
    return valor

def frecuencias(texto, k = 1):
    resultado = []
    diccionario = {}
    for car in texto:
	if car in diccionario:
	    diccionario[car].otro_mas()
	else:
	    a = Elemento(car, 1)
	    resultado.append(a)
	    diccionario[car] = a
    resultado.sort()
    resultado.reverse()
    return resultado

def suma_frecuencias(lista):
    l = [0]
    for i in range(1,len(lista)):
        l.append(l[i-1] + lista[i-1].getProbabilidad())
    return l

def probabilidades(lista):
    lista_frecuencias = [lista[i].getFrecuencia() for i in range(len(lista))]
    n = float(sum(lista_frecuencias))
    for i in range(len(lista)):
        lista[i].probabilidad = lista[i].getFrecuencia()/n

def to_binary(valor, l):
    i = 0
    st = ""
    while i < l:
        valor = 2*valor
        st += str(int(valor))
        valor = valor - int(valor)
        i = i+1
    return st

def codigo_shannon(lista):
    diccio = {}
    freq = suma_frecuencias(lista)
    for i in range(len(lista)):
        long = n(-log(lista[i].getProbabilidad(), 2))
        if not(long - int(long) == 0):
            long = int(long) + 1
        diccio[lista[i].getCaracter()] = to_binary(freq[i], long)
    return diccio        

class Shannon:
    def __init__(self, texto, k = 1):
        self.freq = frecuencias(texto, k)
        probabilidades(self.freq)
        self.codigo = codigo_shannon(self.freq)

    def getCodigo(self, car):
        return self.codigo[car]
        
    def entropia(self):
        ret = 0
        for i in self.freq:
            ret += i.getProbabilidad()*log(i.getProbabilidad(), 2)
        return -n(ret)

class MiHuffman:
    def __init__(self, texto, k = 1):
	self.freq = frecuencias(texto, k)
	listado = {}
	for i in self.freq:
	    listado[i.getCaracter()] = i.getFrecuencia()
	self.H = Huffman(table=listado)
	self.codigo = self.H.encoding_table()

    def getCodigo(self, car):
	return self.codigo[car]

class Escribir:
    def __init__(self, shannon, fitx = None):
        self.shannon = shannon
        self.fichero = fitx
        self.texto = ""
        self.numero = 0
        self.digitos = 0
        self.longitud = 0
        
    def codificar(self, texto):
        for b in texto:
            self.escribir(b)
        self.completarTexto()
        self.fichero.close()
        
        
    def escribir(self, caracter):
        cod = self.shannon.getCodigo(caracter)
        for letra in cod:
            self.longitud += 1
            self.numero *= 2
            if letra == '1':
                self.numero += 1
            self.digitos += 1
            if self.digitos == 8:
                if self.fichero == None:
                    self.texto += chr(self.numero)
                else:
                    self.fichero.write(chr(self.numero))
                self.numero = 0
                self.digitos = 0
                
    def completarTexto(self):
        if not (self.digitos == 0):
            for i in range(self.digitos, 8):
                self.numero *= 2
            if self.fichero == None:
                self.texto += chr(self.numero)
            else:
                self.fichero.write(chr(self.numero))
            self.numero = 0
                
    def returnTexto(self):
        return self.texto
        
def ocupacion_codigo(codigo):
    max = 0
    for palabra in codigo:
        i = len(palabra)*8 + len(codigo[palabra])
        if max < i:
            max = i
    return [max, len(codigo), max*len(codigo)]

def longitud_codigo(codigo):
    ocupacion = 1
    for palabra in codigo:
	ocupacion = ocupacion + (len(palabra)+2)*8 + len(codigo[palabra])
    ocupacion = ocupacion
    return ocupacion

def longitud_media_codigo_freq(palabra, freq):
    for p in freq:
	if p.getCaracter() == palabra:
	    return p.getProbabilidad()

def longitud_media_codigo(freq, codigo):
    media = 0
    for palabra in codigo:
	media += longitud_media_codigo_freq(palabra, freq)*len(codigo[palabra])
    return media

def crear_dos_compresiones(fichero, k = 1):
    fichero.seek(0)
    texto = leer(fichero, k)
    S = Shannon(texto, k)
    H = MiHuffman(texto, k)
    return texto, S, H

def abrir(nombre):
    return codecs.open(nombre, "r", "utf-8")

def toda_informacion(fichero, k = 1):
    texto, S, H = crear_dos_compresiones(fichero, k)
    print "\n\n"
    print "#" * 10
    print "Estamos en el caso k = " + str(k)
    print "#" * 10
    print "\n\n"
    print "Un total de " + str(len(S.freq)) + " grupos diferentes"
    a = S.freq[0].getFrecuencia()/S.freq[0].getProbabilidad()
    print "entre aproximadamente " + str(a) + " grupos."
    entropia = S.entropia()
    print "Entropia: " + str(entropia)
    print "Longitud media codigo SHANNON: " + str(longitud_media_codigo(S.freq, S.codigo))
    print "Longitud media codigo HUFFMAN: " + str(longitud_media_codigo(H.freq, H.codigo))
    longitud_esperada = n(entropia*a/8)
    print "Longitud esperada (texto codificado): " + str(longitud_esperada)
    print "10 grupos mas utilizados:"
    grupos = S.freq[0:10]
    print "Grupo / Shannon / Huffman / Diferencia tamanio"
    for p in grupos:
        print "'" + p.getCaracter() + "'", S.codigo[p.getCaracter()], H.codigo[p.getCaracter()], len(S.codigo[p.getCaracter()]) - len(H.codigo[p.getCaracter()])
    print "10 grupos menos utilizados:"
    for i in range(1, 10):
	print "'" + S.freq[-i].getCaracter() + "'", S.codigo[S.freq[-i].getCaracter()], H.codigo[S.freq[-i].getCaracter()], len(S.codigo[S.freq[-i].getCaracter()]) - len(H.codigo[S.freq[-i].getCaracter()])
    print "Longitud codigo Shannon (bits / bytes)"
    a = longitud_codigo(S.codigo)
    print a, n(a/8)
    print "Longitud codigo Huffman (bits / bytes)"
    a = longitud_codigo(H.codigo)
    print a, n(a/8)
    outputS = open("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 2/outputShannon" + str(k) + ".txt", "w")
    outputH = open("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 2/outputHuffman" + str(k) + ".txt", "w")
    E = Escribir(S, outputS)
    E.codificar(texto)
    E = Escribir(H, outputH)
    E.codificar(texto)
    print "Longitud texto codificado Shannon (bytes)"
    print os.path.getsize("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 2/outputShannon" + str(k) + ".txt")
    print "Longitud texto codificado Huffman (bytes)"
    print os.path.getsize("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 2/outputHuffman" + str(k) + ".txt")

def haz_todo(nombre):
    fichero = abrir(nombre)
    toda_informacion(fichero, 1)
    toda_informacion(fichero, 2)
    toda_informacion(fichero, 3)
    toda_informacion(fichero, "word")

print "No te olvides de"
print "import os"
print "import codecs"
print "from sage.coding.source_coding.huffman import Huffman"
