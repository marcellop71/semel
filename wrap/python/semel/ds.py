import random
import math
import hashlib

import numpy as np

from decimal import *

def modf(x, y):
    x_mod_y = x % y
    if (x_mod_y < 0):
        x_mod_y = y - Decimal.copy_abs(x_mod_y)
    return x_mod_y

def mod2pi(x):
    return modf(x, Decimal(2.0 * math.pi))

def mod1(x):
    return modf(x, Decimal(1.0))

def to_bytes(s, n):
    s_ = int.from_bytes(s, 'little')
    return (s_).to_bytes(n, byteorder='little')

def bytes_to_dec(s):
    s_ = int.from_bytes(s, 'big')
    s_ = Decimal(str('0.' + str(s_)))
    return s_

def dec_to_bytes(s, n):
    s_ = int(str(s)[2:])
    return (s_).to_bytes(n, byteorder='big')

def chunk(s, chunk_size):
    chunks = []
    a = len(s) // chunk_size
    for i in range(0, a):
        chunks.append(s[chunk_size*i:chunk_size*(i+1)])
    chunks.append(to_bytes(s[chunk_size*a:], chunk_size))
    return chunks

class DDS():
    def iterate(self, num_iter, track = False):
        if (track == False):
            for _ in range(0, num_iter):
                self.iterate_()
        else:
            self.orbit = []
            for _ in range(0, num_iter):
                self.orbit.append(self.iterate_())
        return 0

    def get_point_cloud(self):
        orbit_ = [[float(x) for x in tmp] for tmp in self.orbit]
        return np.array(orbit_)

class Logistic(DDS):
    def __init__(self, K):
        self.K = Decimal(K)
        self.x = Decimal(0)
        return None

    def set_start(self, x0):
        self.x = Decimal(x0)
        return 0

    def set_start_random(self):
        x0 = mod1(Decimal(random.random()))
        return self.set_start(x0)

    def iterate_(self):
        xtmp = self.K * self.x * (Decimal(1.0) - self.x)
        self.x = mod1(xtmp)
        return [self.x]

    def hashing(self, s, n, num_iter):
        q = Decimal(0)
        for s_chunk in chunk(s, n):
            p = bytes_to_dec(s_chunk)
            self.set_start(p)
            self.iterate(num_iter)
            q = mod1(q + self.x)
        digest = dec_to_bytes(q, n)
        return digest

class Logistic2D(DDS):
    def __init__(self, K):
        self.K = Decimal(K)
        self.x = Decimal(0)
        self.y = Decimal(0)
        return None

    def set_start(self, x0, y0):
        self.x = Decimal(x0)
        self.y = Decimal(y0)
        return 0

    def set_start_random(self):
        x0 = mod1(Decimal(random.random()))
        y0 = mod1(Decimal(random.random()))
        return self.set_start(x0, y0)

    def iterate_(self):
        self.x = mod1(self.K * ((Decimal(3) * self.y) - 1) * self.x * (Decimal(1.0) - self.x))
        self.y = mod1(self.K * ((Decimal(3) * self.x)) *self.y * (Decimal(1.0) - self.y))
        return [self.x, self.y]

class Henon(DDS):
    def __init__(self, a, b):
        self.a = Decimal(a)
        self.b = Decimal(b)
        self.x = Decimal(0)
        self.y = Decimal(0)
        return None

    def set_start(self, x0, y0):
        self.x = Decimal(2 * x0 - 1)
        self.x = Decimal(2 * y0 - 1)
        return 0

    def set_start_random(self):
        x0 = mod1(Decimal(random.random()))
        y0 = mod1(Decimal(random.random()))
        return self.set_start(x0, y0)

    def iterate_(self):
        xtmp = Decimal(1.0) - self.a * (self.x * self.x) + self.y
        ytmp = self.b * self.x
        self.x = mod1(xtmp)
        self.y = mod1(ytmp)
        return [self.x, self.y]

class CSM(DDS):
    def __init__(self, K):
        getcontext().prec = 80
        self.K = Decimal(K)
        self.x = Decimal(0)
        self.p = Decimal(0)
        return None

    def set_start(self, x0, p0):
        self.x = Decimal(x0)
        self.p = Decimal(p0)
        return 0

    def set_start_random(self):
        x0 = mod2pi(Decimal(random.random() * (2 * math.pi)))
        p0 = mod2pi(Decimal(random.random() * (2 * math.pi)))
        return self.set_start(x0, p0)

    def iterate_(self):
        p = self.p + (self.K * Decimal(math.sin(self.x)))
        x = self.x + p
        self.x = mod2pi(x)
        self.p = mod2pi(p)
        return [self.x, self.p]
