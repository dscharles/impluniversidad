from .BasicIdent import BasicIdent
import random
from sage.crypto.cryptosystem import PublicKeyCryptosystem
from sage.all import EllipticCurve
from sage.all import Hom
from sage.all import Zmod, FiniteField, Integer
from copy import deepcopy

class FullIdent(BasicIdent):
    def H4(self, text, length = 0):
        random.seed(hash("".join(map(str,text))))
        mask = [None]*length
        for i in xrange(length):
            mask[i] = random.choice([0, 1])
        return mask

    def H3(self, text):
        random.seed(hash("".join(map(str,text))))
        return random.randint(2, self.order-1)

    def encrypt(self, message, pubkey, seed=None, text=False):
        random.seed(seed)
        
        tmp = None
        if not text:
            tmp = Integer(message).digits(2)
        else:
            tmp = 0
            for let in message:
                tmp = tmp*256
                tmp = (tmp + ord(let))
                
            tmp = Integer(tmp).digits(2)
            
        tmp.reverse()

	l = len(tmp)
        sigma = [0]*l
        for i in xrange(l):
            sigma[i] = random.choice([0,1])

        r = self.H3("".join(map(str,sigma)) + message)
        if self.pairing == "tate":
            pair = self._ext(pubkey[0]).tate_pairing(self.distortion(pubkey[1]), self.order, self.k, self.ec2.base_ring().cardinality())
        else:
            pair = self._ext(pubkey[0]).weil_pairing(self.distortion(pubkey[1]), self.order)
        
        print "Sin cifrar", tmp
        c1 = self._mask(list(sigma), pair**r, self.H2)
        c2 = self._mask(tmp, sigma, self.H4)
        return r*self.P, c1, c2

    def decrypt(self, ciphertext, privatekey, text=False):
        if self.pairing == "tate":
            pair = self._ext(privatekey).tate_pairing(self.distortion(ciphertext[0]), self.order, self.k, self.ec2.base_ring().cardinality())
        else:
            pair = self._ext(privatekey).weil_pairing(self.distortion(ciphertext[0]), self.order)
            
        sigma = self._mask(ciphertext[1], pair)
        msg = self._mask(sigma, self.H4)
        msg = int(self._mask(map(int, list(msg)), pair), base=2)
        if text:
            msg = map(chr, Integer(msg).digits(256))
            msg.reverse()
            msg = "".join(msg)
        r = self.H3("".join(sigma) + msg)
        return r*self.P == ciphertext[0], msg