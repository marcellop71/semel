import numpy as np
import random

class TimeSeries:

    def __init__(self, points=None):
        self.data = points
        self.cloud = None

    def __del__(self):
        pass

    def embedding(self, dim, tau, method, mute = False):
        x=self.data

        if method == 0:
            #[x0, x1, x2, ..]
            # maps to
            #[[x0,x(tau)...x(tau*(dim-1))],
            # [x1, x(tau+1), ... x(tau*(dim-1))]
            # ...]

            n=len(x)-(dim-1)*tau
            assert (n>=0), 'Not enough data for this embedding'
            P=[[x[i+j*tau] for j in range(dim)] for i in range(n)]

        elif method == 1:
            #[x0, x1, x2, ..]
            # maps to
            #[[x0, x(tau), .. x(tau*(dim-1))],
            # [x(tau*dim), x(tau*(dim+1)), ... x(tau*(2*dim-1))]
            # ...]

            n=len(x)-dim*tau
            assert (n>=0), 'Not enough data for this embedding'
            n=len(x)//(dim*tau)
            P=[[x[(i+j*dim)*tau] for i in range(dim)] for j in range(n)]

        elif method == 2:
            #[x0, x1, x2, ..]
            # maps to
            #[[x0, x1, .. x(dim-1)],
            # [x(tau), x(tau+1), ... x(tau+dim-1)]
            # ...]

            n=len(x)-dim
            assert (n>=0), 'Not enough data for this embedding'
            n=max((n-tau)//tau, 0)
            P=[[x[i+j*tau] for i in range(dim)] for j in range(n)]

        else:
            raise Exception('Please insert a valid method (0-2).')

        if mute:
            self.cloud = np.array(P)
        else:
            return np.array(P)


class IET(TimeSeries):

    def __init__(self, permutation = None, weight = None, start = None):
        TimeSeries.__init__(self)
        if permutation is not None:
            self.permutation  = permutation
        else:
            permutation = list(range(5))
            random.shuffle(permutation)
            self.permutation = permutation

        if weight is not None:
            self.weight = weight
        else:
            w=[random.random() for i in range(len(self.permutation))]
            self.weight=[i/sum(w) for i in w]

        if start is not None:
            self.start = start
        else:
            self.start = random.random()

    def __del__(self):
        pass

    def evolution(self, n, mute=False):
        a=[sum(self.weight[:i]) for i in range(len(self.weight))]
        b=[sum([self.weight[self.permutation.index(j)] for j in range(self.permutation[i])]) \
           for i in range(len(self.permutation))]
        y=[self.start]
        for j in range(n):
            flag=True
            i=0
            while flag and i < len(a):
                if a[i] < y[-1]:
                    i+=1
                else:
                    flag=False
            i-=1
            y.append(y[-1] - a[i] + b[i])

        if mute:
            self.data = y
        else:
            return y


class Weierstrass(TimeSeries):
    def __init__(self, w=None, s=None, iteration=None):
        if w is not None:
            self.weight = w
        else:
            self.weight = random.random() + 1

        if s is not None:
            self.s = s
        else:
            self.s = random.random() + 1

        if iteration is None:
            iteration = 1000

        self.function = lambda t: sum( self.weight ** ( (self.s -2)*k) * np.sin(t*self.weight ** k) for k in range(iteration))

    def __del__(self):
        pass

    def evaluate(self, points):
        self.data = [self.function(i) for i in points]
        return


class SelfAffine(TimeSeries):

    def __init__(self, m=None, b=None, r=None):
        if m is None:
            m = random.random()
        if b is None:
            b = random.random()
        if r is None:
            r = [(1-m) * random.random() + m,\
                 (1-m) * random.random() + m]

        self.mat = [ np.array([[m, 0],[((-1)**i)*b, r[i]]])\
                            for i in range(2)]
        self.tra = [ i*np.array([m, b]) for i in range(2)]
        self.maps = [ lambda v: self.mat[0]@v + self.tra[0],\
                      lambda v: self.mat[1]@v + self.tra[1]]

    def __del__(self):
        pass

    def plot(self, iteration, store=True, show=True):
        last = [self.maps[1](np.array([0,0]))]
        result = last
        result = np.concatenate((result, np.array([[0,0],[1,0]])), axis = 0)
        for i in range(iteration):
            aux=[]
            for point in last:
                for f in self.maps:
                    aux.append(f(point))
            last = aux
            result = np.concatenate((result, aux), axis = 0)

        expansion = [(i[1][0], i[0]) for i in enumerate(result)]
        expansion.sort()
        index=[ i[1] for i in expansion]
        result = result[index]
        if show:
            f, ax = plt.subplots(1)
            ax.plot([i[0] for i in result], [i[1] for i in result])
            ax.set_ylim(bottom=0)
            plt.show(f)
        if store:
            self.data = np.array(result)
            return
        else:
            return np.array(result)
