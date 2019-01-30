import numpy as np

# merge persistent interval and cycles
def coalesce(z, pi, rmax):
    d = {}

    for key_ in pi.keys():
        if (key_ not in d):
            d[key_] = []
        for i, [x, y] in enumerate(pi[key_]):
            d[key_].append([x, y])

    for key_ in z.keys():
        if (key_ not in d):
            d[key_] = []
        for i, [x] in enumerate(z[key_]):
            d[key_].append([x, rmax])

    for key_ in d.keys():
        d[key_] = np.asarray(d[key_])

    return d

def entropy(d):
    se = 0.0
    for p in np.nditer(d):
        se += - p * np.log(p)
    return se

# weighted Betti numbers
# add gaussian_blur
def persistent_1d_distribution(d):
    m = 0.0
    for i, [x, y] in enumerate(d):
        m += y - x

    d_ = []
    for i, [x, y] in enumerate(d):
        d_.append((y - x) / m)
    d_ = np.asarray(d_)
    return [d_, m]

# persistent diagram (binning 2-dim persistent diagram)
def persistent_2d_distribution(d):
    pass

# dei 2 1d distributions -> color strip (1d spectrum)
# dei 2 2d distributions -> color image (2d spectrum)

# distanza tra 2 diagrammi di persistenza
# -> distanza tra le 2 2d-distributions associate (sia delta che gaussian blur o altro)
# -> ora nel transport si porta sempre A in B

def transport(A, B):
    import ot
    [n, d] = A.shape
    a, b = np.ones((n,)) / n, np.ones((n,)) / n

    M = ot.dist(A, B)
    M /= M.max()

    G = ot.emd(a, b, M)
    return G

def persistent_landscape(t, pi):
    a = []
    for t_ in t:
        tmp1 = t_ - pi[:,0]
        tmp2 = pi[:,1] - t_
        q = np.min([tmp1, tmp2], axis=0)
        q_ = np.sort(q)
        q_[q_<0] = 0
        a.append(q_[-1:-5:-1])
    a = np.array(a)
    return [t, a]
