import numpy as np
import hashlib

def hashsequence(n, seed, salt):
    x = [bytearray(seed + salt, 'utf-8')]
    for k in range(0, n):
        x.append(hashlib.sha256(x[-1]).digest())
    return x[1:]

class PointCloud:

    def __init__(self):
        pass

    def __del__(self):
        pass

    def random(self, n, d):
        P = 2.0 * np.random.rand(n, d) - 1.0
        return P

    def grid(self, n, eps):
        P = []
        for x in np.arange(-1.0, 1.0, 2.0 / n):
            for y in np.arange(-1.0, 1.0, 2.0 / n):
                e = eps * (2.0 * np.random.rand(1, 2) - 1.0)
                P.append([x + e[0,0],y + e[0,1]])
        P = np.array(P)
        return P

    def grid_random_0(self, n):
        P = []
        x_ = 2.0 * np.random.rand(n, 1) - 1.0
        for x in x_:
            for y in np.arange(-1.0, 1.0, 2.0 / n):
                P.append([x,y])
        P = np.array(P)
        return P

    def grid_random_1(self, n):
        P = []
        x_ = 2.0 * np.random.rand(n, 1) - 1.0
        y_ = 2.0 * np.random.rand(n, 1) - 1.0
        for x in x_:
            for y in y_:
                P.append([x,y])
        P = np.array(P)
        return P

    def circle(self, n, cx, cy, rmin, rmax, method):
        P = []; i = 0
        if (method == 0):
            while (i < n):
                rtmp = rmin + np.random.random() * (rmax - rmin)
                thetatmp = np.random.random() * (2.0 * np.pi)
                ptmp = [cx + rtmp * np.cos(thetatmp), cy + rtmp * np.sin(thetatmp)]
                P.append(ptmp); i += 1
        if (method == 1):
            while (i < n):
                rtmp = rmin + np.random.random() * (rmax - rmin)
                thetatmp = (i / n) * (2.0 * np.pi)
                ptmp = [cx + rtmp * np.cos(thetatmp), cy + rtmp * np.sin(thetatmp)]
                P.append(ptmp); i += 1
        P = np.array(P)
        return P

    def circle_q(self, n, cx, cy, rmin, rmax, ravg, rstd):
        P = []; i = 0
        rtmp = ravg
        while (i < n):
            rtmp += np.random.randn() * rstd
            if (rtmp < rmin):
                rtmp = rmin
            if (rtmp > rmax):
                rtmp = rmax
            thetatmp = (i / n) * (2.0 * np.pi)
            ptmp = [cx + rtmp * np.cos(thetatmp), cy + rtmp * np.sin(thetatmp)]
            P.append(ptmp); i += 1
        P = np.array(P)
        return P

    def sphere(self, n, rmin, rmax):
        P = []; i = 0
        while (i < n):
            rtmp = rmin + np.random.random() * (rmax - rmin)
            thetatmp = 1.0 + (np.random.random() * 0.20 * np.pi)
            phitmp = np.random.random() * np.pi
            ptmp = [rtmp * np.sin(thetatmp) * np.cos(phitmp), rtmp * np.sin(thetatmp) * np.sin(phitmp), np.cos(thetatmp)]
            P.append(ptmp); i += 1
        P = np.array(P)
        return P

    def rnd_0(self, n):
        eps = np.random.random()
        if (eps > 0.50):
            d = 2
            P = pc_random(n, d)
        else:
            cx = 0.0; cy = 0.0; rmin = 1.00; rmax = 1.00; method = 0
            P = pc_circle(n, cx, cy, rmin, rmax, method)
        return P

    def hashsequence_embcube(self, nP, d, seed, salt):
        x = hashsequence(nP, seed, salt)
        bs = int(len(x[0]) / d)
        P = []
        for x_ in x:
            Ptmp = []
            for k in range (0, d):
                b = x_[k*bs:(k+1)*bs] # x_[k::bs]
                #tmp = float('0.' + str(int.from_bytes(b, byteorder='little', signed=False)))
                tmp = float.fromhex('0x0.' + b.hex())
                Ptmp.append(tmp)
            P.append(Ptmp)
        P = np.array(P)
        return P

    def fetch_from_pdb(self, folder, item):
        f = gzip.open('%s/%s' % (folder, item), 'rt')
        X = []
        for s in f:
            if (s[0:6] == 'ATOM  '):
                serial = int(s[6:12])
                resseq = int(s[22:26])
                x = float(s[30:38])
                y = float(s[38:46])
                z = float(s[46:54])
                X.append([serial, resseq, x, y, z])
        X = np.array(X)
        P = np.column_stack([X[:,2], X[:,3], X[:,4]])
        return [X, P]
