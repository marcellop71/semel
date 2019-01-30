import matplotlib
matplotlib.use('Agg')

import os
import numpy as np
import torch

from semel import Semel, PointCloud, Graph

if __name__ == '__main__':
    os.makedirs('output', exist_ok=True)

    n = 8
    P = PointCloud().grid(n, 0.00)
    Graph().point_cloud(P, None, labels=True, xlim=[-1.2, 1.2], base_folder='output', name='grid_coboundary_pc')

    semel = Semel()

    pci = 0
    manifold = {
        'type': semel.MANIFOLD_FLAT_UNBOUNDED,
        'distance': semel.DISTANCE_L2,
        'border': []
    }
    semel.add_point_cloud(pci, manifold, P)
    semel.add_complex(pci)

    S = [
        [9, 16],
        [16, 25],
        [25, 26],
        [26, 27],
        [9, 18],
        [18, 27],
    ]
    typeS = 0
    semel.add_simplices(pci, typeS, S)

    dim_max = 2
    r_max = 0.50
    semel.build_complex(pci, dim_max, r_max)

    q_complex = semel.export_complex(pci)
    Graph().complex(P, q_complex, labels=True, xlim=[-1.2, 1.2], base_folder='output', name='grid_coboundary_complex')
    print(q_complex)

    q_cohomology = semel.cohomology(pci, measures_only=False)

    Graph().bundle_0(dim_max, r_max, q_cohomology, base_folder='output', name='grid_coboundary_bundle')

    q_cocycles = semel.export_cocyles(pci)
    Graph().complex_and_cocycles(P, q_complex, q_cocycles, labels=False, xlim=[-1.2, 1.2],
                                 base_folder='output', name='grid_coboundary_cocycles')

    q = semel.export_coboundary(pci, 1)
    print(q)

    Graph().complex(P, q_complex, labels=True, xlim=[-1.2, 1.2], base_folder='output', name='grid_coboundary_complex2')
    print(q_complex)

    for dim in [0, 1]:
        q = semel.export_coboundary(pci, dim)
        total = q['m'] * q['n']
        nnz = len(q['matrix_sparse'])
        pct = 100 * nnz / total if total > 0 else 0.0
        print("d^%d: C^%d (R^%d) -> C^%d (R^%d), nnz=%d / %d = %.2f%%" % (
            dim, dim, q['n'], dim + 1, q['m'], nnz, total, pct))

    # Torch cocycle optimisation
    q_complex = semel.export_complex(pci)
    q_cocycles = semel.export_cocyles(pci)
    q_cocycles_ = []
    for iz, z in enumerate(q_cocycles):
        if iz > 1:
            break
        if z['dimension'] == 1:
            print(z)
            q_coboundary = semel.export_coboundary(pci, 0)

            batch_size = 1
            x = torch.randn(q_coboundary['n'], batch_size, dtype=torch.float, requires_grad=True)

            b = torch.zeros((q_coboundary['m'], batch_size), dtype=torch.float)
            for tmp in z['expression']:
                idx = q_complex['complex'][tmp[1]]['index']
                b[idx, :] = tmp[0]
            print(b.pow(2).sum())

            A_ = np.array(q_coboundary['matrix_sparse'], dtype=np.int32)
            i = torch.LongTensor([A_[:, 1], A_[:, 0]])
            v = torch.FloatTensor(A_[:, 2])
            A = torch.sparse_coo_tensor(i, v, torch.Size([q_coboundary['m'], q_coboundary['n']]))
            print(A.to_dense())

            learning_rate = 1e-3
            for k in range(10000):
                y = torch.sparse.mm(A, x) + b
                loss = y.pow(2).sum()
                if (k % 1000) == 0:
                    print(loss.item())

                loss.backward()
                with torch.no_grad():
                    x -= learning_rate * x.grad
                    x.grad.zero_()

            print(len(torch.nonzero(x)[:, 0]))
            print(y)

            tmp_expr = []
            for idx in torch.nonzero(x)[:, 0]:
                tmp_expr.append([float(x[idx, 0].data), int(idx.data)])
            tmp = {'dimension': 0, 'expression': tmp_expr}
            q_cocycles_.append(tmp)

    print(q_cocycles_)

    # Hodge-Laplace eigenvectors
    q = semel.HodgeLaplace(pci, 2)
    v = q['eigenvectors'][0][0]
    print(v)
    print(np.sum(v * v))

    semel.unload_complex(pci)
    semel.unload_point_cloud(pci)
