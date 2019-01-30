import matplotlib
matplotlib.use('Agg')

import os
import numpy as np
import matplotlib.pyplot as plt

from semel import Semel, PointCloud, Graph

if __name__ == '__main__':
    os.makedirs('output', exist_ok=True)

    semel = Semel()

    def pc_measures(P):
        pci = 0
        manifold = { 'type': semel.MANIFOLD_FLAT_UNBOUNDED, 'distance': semel.DISTANCE_L2, 'border': [] }
        semel.add_point_cloud(pci, manifold, P)
        semel.add_complex(pci)

        semel.generate_delaunay(pci)

        dim_max = 2
        r_max = 1.0
        semel.build_complex(pci, dim_max, r_max)
        q_cohomology = semel.cohomology(pci, measures_only=False)
        semel.unload_complex(pci)
        semel.unload_point_cloud(pci)
        return q_cohomology

    x = np.linspace(0.0, 10 * 2 * np.pi, 400)
    s = np.concatenate([np.sin(x), np.sin(2.0 * x), np.sin(x) + np.sin(2.0 * x), np.sin(x) + np.sin(4.0 * x)])

    r = s

    wlen = 2 * 60
    K_ = range(0, len(s) - wlen)

    d = 2
    mass = np.zeros((len(K_), d))
    entropy = np.zeros((len(K_), d))

    for k in range(0, len(s) - wlen):
        if (k % 10) == 0:
            print(k)
        r_ = r[k:(k+wlen)]
        tau = 1
        P = np.array([r_[0:-tau], r_[tau:]]).transpose()
        if np.max(P) > np.min(P):
            P = 2.0 * ((P - np.min(P)) / (np.max(P) - np.min(P))) - 1.0
            q_cohomology = pc_measures(P)

            for i in range(0, d):
                if i in q_cohomology['1d']['mass']:
                    mass[k, i] += q_cohomology['1d']['mass'][i]
                if i in q_cohomology['1d']['entropy']:
                    entropy[k, i] += q_cohomology['1d']['entropy'][i]

    fig = plt.figure(figsize=(20, 20))

    ax00 = fig.add_subplot(321)
    ax00.scatter(K_, s[wlen:], s=8)
    ax00.plot(K_, s[wlen:])
    ax00.grid()

    ax01 = fig.add_subplot(322)
    ax01.scatter(K_, s[wlen:], s=8)
    ax01.plot(K_, s[wlen:])
    ax01.grid()

    ax0 = fig.add_subplot(323)
    ax0.scatter(K_, mass[:,0], s=8)
    ax0.plot(K_, mass[:,0])
    ax0.grid()

    ax1 = fig.add_subplot(324)
    ax1.scatter(K_, entropy[:,0], s=8)
    ax1.plot(K_, entropy[:,0])
    ax1.grid()

    ax2 = fig.add_subplot(325)
    ax2.scatter(K_, mass[:,1], s=8)
    ax2.plot(K_, mass[:,1])
    ax2.grid()

    ax3 = fig.add_subplot(326)
    ax3.scatter(K_, entropy[:,1], s=8)
    ax3.plot(K_, entropy[:,1])
    ax3.grid()

    plt.savefig('output/timeseries_rolling.jpg', dpi=150, bbox_inches='tight')
    plt.close()
