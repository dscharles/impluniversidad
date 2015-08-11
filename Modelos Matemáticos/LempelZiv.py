from sage.all import *

class LempelZiv:
    def __init__(self, texto, output):
        self.texto = texto
        self.texto.seek(0)
        self.fichero = output
        self.bloques = []
        self.bloques2 = {}
        self.N = 1
        self.actual = ""
        self.digitos = 0
        self.numero = 0
        
    def leer_bloque(self):
        tmp = ""
        continuar = True
        while continuar:
            self.actual = self.texto.read(1)
            if self.actual == "":
                return tmp
            tmp += self.actual
            self.actual = self.actual[1:len(self.actual)]
            if not (tmp in self.bloques2):
                self.bloques.append(tmp)
                self.bloques2[tmp] = self.N
                self.N += 1
                return tmp
                
    def leer_texto(self):
        o = self.leer_bloque()
        while o != "":
            o = self.leer_bloque()
                
    def escribir2(self, caracter):
        for letra in caracter:
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
                
    def escribir(self, indice, digito):
        binario = Integer(indice).binary()
        binario = '0'*(self.N-len(binario)) + binario + Integer(ord(digito)).binary()
        self.escribir2(binario)
                
    def codificar(self):
        self.leer_texto()
        self.N = int(log(len(self.bloques), 2))+1
        for j in range(len(self.bloques)):
            try:
                i = self.bloques.index(self.bloques[j][0:len(self.bloques[j])-1]) + 1
            except:
                i = 0
            self.escribir(i, self.bloques[j][-1])
            cod = self.leer_bloque()
        self.fichero.close()
