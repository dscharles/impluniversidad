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

class Codigo:
    def __init__(self):
        self.texto = {}
        self.texto2 = []
        self.codigo = {}
	self.longitud_extra = 0
        self.freq = []
        self.caracteres = "aáàÁÀAbBcCdDeéèÉÈEfFgGhHiíìÍÌIjJkKlLmMnNoóòÓÒOpPqQrRsStTuúùÚÙUvVwWxXyYzZñÑçÇ"
        
    def insertar(self, palabra):
        if palabra in self.texto:
            self.texto[palabra] += 1
        else:
            self.texto[palabra] = 1
            
    def escaracter(self, letra):
        return letra in self.caracteres
        
    def leer(self, texto, k = 1):
        texto.seek(0)
        if not k == "word":
            st = texto.read(k)
            while not st == "":
                self.insertar(st)
                self.texto2.append(st)
                st = texto.read(k)
        else:
            st = " "
            while not st == "":
                caracter = texto.read(1)
                st = ""
                while not caracter == "" and self.escaracter(caracter):
                    st += caracter
                    caracter = texto.read(1)
                if st != "":
                    self.insertar(st)
                    self.texto2.append(st)
                if caracter != "":
                    self.insertar(caracter)
                    self.texto2.append(caracter)
                    st = caracter
                else:
                    st = ""
                    
    def frecuencias(self):
        self.freq = [ Elemento(car, self.texto[car]) for car in self.texto ]
        self.freq.sort()
        self.freq.reverse()
        
    def probabilidades(self):
        lista_freq = [ self.freq[i].getFrecuencia() for i in range(len(self.freq)) ]
        n = float(sum(lista_freq))
        for i in range(len(self.freq)):
            self.freq[i].probabilidad = self.freq[i].getFrecuencia()/n
            
    def fin_leer(self):
        self.frecuencias()
        self.probabilidades()
        
    def getCodigo(self, caracter):
        if caracter in self.codigo:
            return self.codigo[caracter]
        else:
            st = ""
            for i in caracter:
                bin = Integer(ord(i)).binary()
                st += "0"*(8-len(bin)) + bin
		self.longitud_extra += 8
            self.longitud_extra += 8
            return st
            
    def entropia(self):
        ret = 0
        for i in self.freq:
            ret += i.getProbabilidad()*log(i.getProbabilidad(), 2)
        return -n(ret)
        
    def esperanza_codigo(self):
        ret = 0
        for i in self.freq:
            ret += i.getProbabilidad()*len(self.getCodigo(i.getCaracter()))
        return ret
        
    def longitud_codigo(self):
        ocupacion = 0
        for palabra in self.codigo:
            ocupacion += (len(palabra)+2)*8 + len(self.codigo[palabra])
	print "Han hecho falta " + str(self.longitud_extra) + " bits de más."
	ocupacion += self.longitud_extra
        return ocupacion

class Shannon(Codigo):
    def suma_frecuencias(self):
        l = [0]
        for i in range(1, len(self.freq)):
            l.append(l[i-1] + self.freq[i-1].getProbabilidad())
        return l
        
    def to_binary(self, valor, l):
        i = 0
        st = ""
        while i < l:
            valor = 2*valor
            st += str(int(valor))
            valor = valor - int(valor)
            i = i+1
        return st
        
    def generar_codigo(self):
        freq = self.suma_frecuencias()
        for i in range(len(self.freq)):
            long = n(-log(self.freq[i].getProbabilidad(), 2))
            if not(long - int(long) == 0):
                long = int(long) + 1
            self.codigo[self.freq[i].getCaracter()] = self.to_binary(freq[i], long)

class MiHuffman(Codigo):
    def generar_codigo(self):
        try:
            H = Huffman(table = self.texto)
            self.codigo = H.encoding_table()
        except:
            print "Tienes que cargar la siguiente libreria"
            print "from sage.coding.source_coding.huffman import Huffman"

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

def abrir(nombre):
    return codecs.open(nombre, "r", "utf-8")

def toda_informacion(fichero, k = 1):
    S = Shannon()
    S.leer(fichero, k)
    M = MiHuffman()
    M.texto = S.texto
    M.texto2 = S.texto2
    S.fin_leer()
    M.freq = S.freq
    S.generar_codigo()
    M.generar_codigo()
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
    E.codificar(S.texto2)
    E = Escribir(H, outputH)
    E.codificar(S.texto2)
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

def prueba1(fichero1, fichero2, k):
    Sc = Shannon()
    Mc = MiHuffman()
    Sq = Shannon()
    Mq = MiHuffman()
    Sc.leer(fichero1, k)
    Sq.leer(fichero2, k)
    Mc.texto = Sc.texto
    Mc.texto2 = Sc.texto2
    Sc.fin_leer()
    Mc.freq = Sc.freq
    Sc.generar_codigo()
    Mc.generar_codigo()
    Mq.texto = Sq.texto
    Mq.texto2 = Sq.texto2
    Sq.fin_leer()
    Mq.freq = Sq.freq
    Sq.generar_codigo()
    Mq.generar_codigo()
    print "Pruebas normales:"
    print "Obra 1:"
    output = open("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 3/output.txt", "w+")
    E = Escribir(Sc, output)
    E.codificar(Sc.texto2)
    print "Tamaño codificado Shannon: " + str(os.path.getsize("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 3/output.txt"))
    print "Tamaño codigo Shannon: " + str(n(Sc.longitud_codigo()/8))
    output = open("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 3/output.txt", "w+")
    E = Escribir(Mc, output)
    E.codificar(Mc.texto2)
    print "Tamaño codificado Huffman: " + str(os.path.getsize("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 3/output.txt"))
    print "Tamaño codigo Huffman: " + str(n(Mc.longitud_codigo()/8))
    print "Obra 2:"
    output = open("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 3/output.txt", "w+")
    E = Escribir(Sq, output)
    E.codificar(Sq.texto2)
    print "Tamaño codificado Shannon: " + str(os.path.getsize("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 3/output.txt"))
    print "Tamaño codigo Shannon: " + str(n(Sq.longitud_codigo()/8))
    output = open("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 3/output.txt", "w+")
    E = Escribir(Mq, output)
    E.codificar(Mq.texto2)
    print "Tamaño codificado Huffman: " + str(os.path.getsize("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 3/output.txt"))
    print "Tamaño codigo Huffman: " + str(n(Mq.longitud_codigo()/8))
    print "Prueba nueva:"
    output = open("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 3/output.txt", "w+")
    E = Escribir(Sc, output)
    E.codificar(Sq.texto2)
    print "Tamaño codificado SHANNON: " + str(os.path.getsize("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 3/output.txt"))
    print "Tamaño extra: " + str(Sc.longitud_extra)
    output = open("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 3/output.txt", "w+")
    E = Escribir(Mc, output)
    E.codificar(Sq.texto2)
    print "Tamaño codificado HUFFMAN: " + str(os.path.getsize("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 3/output.txt"))
    print "Tamaño extra: " + str(Mc.longitud_extra)

def prueba2(fichero1, fichero2, k)
    S = Shannon()
    M = MiHuffman()
    texto1 = Shannon()
    texto2 = Shannon()
    S.leer(fichero1, k)
    texto1.texto2 = deepcopy(S.texto2)
    S.leer(fichero2, k)
    texto2.leer(fichero2, k)
    M.texto = S.texto
    M.texto2 = S.texto2
    S.fin_leer()
    M.freq = S.freq
    S.generar_codigo()
    M.generar_codigo()
    output = open("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 3/output.txt", "w+")
    E = Escribir(S, output)
    E.codificar(texto1.texto2)
    print "SHANNON TEXTO 1: " + str(os.path.getsize("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 3/output.txt"))
    output = open("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 3/output.txt", "w+")
    E = Escribir(S, output)
    E.codificar(texto2.texto2)
    print "SHANNON TEXTO 2: " + str(os.path.getsize("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 3/output.txt"))
    print "LONGITUD CODIGO SHANNON: " + str(n(S.longitud_codigo()/8))
    output = open("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 3/output.txt", "w+")
    E = Escribir(M, output)
    E.codificar(texto1.texto2)
    print "HUFFMAN TEXTO 1: " + str(os.path.getsize("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 3/output.txt"))
    output = open("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 3/output.txt", "w+")
    E = Escribir(M, output)
    E.codificar(texto2.texto2)
    print "HUFFMAN TEXTO 2: " + str(os.path.getsize("/home/david/Dropbox/Universidad/Cuarto/Modelos/Programa 3/output.txt"))
    print "LONGITUD CODIGO HUFFMAN: " + str(n(M.longitud_codigo()/8))
