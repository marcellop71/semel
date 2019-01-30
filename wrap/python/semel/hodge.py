import numpy as np
import torch
import torch.nn as nn
from torch.autograd import Variable

def modf(x, y):
    x_mod_y = x % y
    if (x_mod_y < 0):
        x_mod_y = y - x_mod_y
    return x_mod_y

def mod1(x):
    return modf(x, 1.0)

def hodge(q_complex, q_cocycles, q_coboundary, learning_rate = 1e-3, max_iter = 20000):
    q_cc = []; y_cc = []; b_cc = []
    for iz, z in enumerate(q_cocycles):
        if (z['dimension'] == 1):
            batch_size = 1
            x = torch.randn(q_coboundary['n'], batch_size, dtype=torch.double, requires_grad=True)

            b = torch.zeros((q_coboundary['m'], batch_size), dtype=torch.double)
            for tmp in z['expression']:
                idx = q_complex['complex'][tmp[1]]['index']
                b[idx,:] = tmp[0]
            b_cc.append(np.array(b.pow(2).sum()))

            A_ = np.array(q_coboundary['matrix_sparse'], dtype=np.int32)
            i = torch.LongTensor([A_[:,0], A_[:,1]])
            v = torch.DoubleTensor(A_[:,2])
            A = torch.sparse.DoubleTensor(i, v, torch.Size([q_coboundary['m'], q_coboundary['n']]))

            for k, t in enumerate(range(max_iter)):
                y = torch.sparse.mm(A, x) + b
                loss = y.pow(2).sum()

                loss.backward()
                with torch.no_grad():
                    x -= learning_rate * x.grad
                    x.grad.zero_()

            y_cc.append(np.array(y.data))

            tmp = []
            for idx in torch.nonzero(x)[:,0]:
                tmp.append(mod1(float(x[idx,0].data)))
            q_cc.append(tmp)

    return [q_cc, y_cc, b_cc]
