import numpy as np

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.lines as lines
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

class Graph():

    def __init__(self):
        pass

    def __del__(self):
        pass

    def point_cloud(self, P, S = None, labels = True, xlim = [0, 1], base_folder = None, name = None, c = 0):

        fig = plt.figure(figsize=(10, 10))

        nP = P.shape[0]; dim = P.shape[1]

        if (dim == 2):

            eps = 0.01
            ax1 = fig.add_subplot(111)

            if (S is not None):
                for s_ in S:
                    if (s_.shape[0] == 2):
                        [x0, y0] = P[s_[0],:]
                        [x1, y1] = P[s_[1],:]
                        ax2.plot([x0, x1], [y0, y1], c='b', linewidth=0.4)
                    if (s_.shape[0] == 3):
                        idx = [s_[0], s_[1], s_[2]]
                        t1 = plt.Polygon(P[idx,:], alpha=0.3, color='gray')
                        plt.gca().add_patch(t1)

            ax1.scatter(P[:,0], P[:,1], s=10, c='b')
            if (labels == True):
                for i, txt in enumerate(range(0, nP)):
                    ax1.annotate(txt, (P[i,0] + eps, P[i,1] + eps))
            ax1.set_xlim([xlim[0], xlim[1]])
            ax1.set_ylim([xlim[0], xlim[1]])
            ax1.set_xlabel("x")
            ax1.set_ylabel("y")
            ax1.set_title("%d points" % (P.shape[0]))
            ax1.grid()

        if (base_folder is not None):
            fig.savefig('%s/%s_seq_%d.jpg' % (base_folder, name, c), bbox_inches='tight')
            plt.close(fig)
        else:
            plt.show()

        return None

    def covering(self, P, rmax, labels = True, xlim = [0, 1], base_folder = None, name = None, c = 0):

        fig = plt.figure(figsize=(10, 10))

        nP = P.shape[0]; dim = P.shape[1]

        if (dim == 2):

            ax1 = fig.add_subplot(111)
            for i in range(0, nP):
                circletmp = plt.Circle((P[i,0], P[i,1]), rmax, color='r', alpha=0.2)
                ax1.add_artist(circletmp)
            scat1 = ax1.scatter(P[:,0], P[:,1], s=10, c='b')
            #if (labels == True):
            #    for i, txt in enumerate(v):
            #        ax1.annotate(txt, (P[i,0],P[i,1]))
            ax1.set_xlim([xlim[0], xlim[1]])
            ax1.set_ylim([xlim[0], xlim[1]])

        if (dim == 3):
            ax1 = fig.add_subplot(321, projection='3d')
            ax1.scatter(P[:,0], P[:,1], P[:,2], s=8, c='b')

        if (base_folder is not None):
            fig.savefig('%s/%s_seq_%d.jpg' % (base_folder, name, c), bbox_inches='tight')
            plt.close(fig)
        else:
            plt.show()

        return None

    def complex(self, P, S, labels = True, xlim = [0, 1], base_folder = None, name = None, c = 0):

        fig = plt.figure(figsize=(10, 10))

        nP = P.shape[0]; dim = P.shape[1]

        eps = 0.01

        if (dim == 2):

            ax2 = fig.add_subplot(111)

            for s in S['complex']:
                s_ = s['simplex']
                if (s['dimension'] == 0):
                    pass
                if (s['dimension'] == 1):
                    [x0, y0] = P[s_[0],:]
                    [x1, y1] = P[s_[1],:]
                    ax2.plot([x0, x1], [y0, y1], c='b', linewidth=0.4)
                if (s['dimension'] == 2):
                    idx = [s_[0], s_[1], s_[2]]
                    t1 = plt.Polygon(P[idx,:], alpha=0.3, color='gray')
                    plt.gca().add_patch(t1)

            ax2.scatter(P[:,0], P[:,1], s=14, alpha=0.6, color='blue')
            if (labels == True):
                for i, txt in enumerate(range(0, nP)):
                    ax2.annotate(txt, (P[i,0] + eps, P[i,1] + eps))
            ax2.set_xlim([xlim[0], xlim[1]])
            ax2.set_ylim([xlim[0], xlim[1]])

        if (dim == 3):
            ax2 = fig.add_subplot(111, projection='3d')

            for s in S:
                s_ = s['simplex']
                if (s['dimension'] == 0):
                    pass
                if (s['dimension'] == 1):
                    [x0, y0, z0] = P[s[0],:]
                    [x1, y1, z1] = P[s[1],:]
                    ax2.plot([x0, x1], [y0, y1], [z0, z1], c='b', linewidth=0.4)
                if (s['dimension'] == 2):
                    [x0, y0, z0] = P[s[0],:]
                    [x1, y1, z1] = P[s[1],:]
                    [x2, y2, z2] = P[s[2],:]
                    ax2.plot_trisurf([x0, x1, x2], [y0, y1, y2], [z0, z1, z2], linewidth=0.2, antialiased=True, alpha=0.3, color='gray')
                if (s['dimension'] == 3):
                    [x0, y0, z0] = P[s[0],:]
                    [x1, y1, z1] = P[s[1],:]
                    [x2, y2, z2] = P[s[2],:]
                    [x3, y3, z3] = P[s[3],:]
                    ax2.plot_trisurf([x0, x1, x2], [y0, y1, y2], [z0, z1, z2], linewidth=0.2, antialiased=True, alpha=0.2, color='green')
                    ax2.plot_trisurf([x0, x1, x3], [y0, y1, y3], [z0, z1, z3], linewidth=0.2, antialiased=True, alpha=0.2, color='green')
                    ax2.plot_trisurf([x0, x2, x3], [y0, y2, y3], [z0, z2, z3], linewidth=0.2, antialiased=True, alpha=0.2, color='green')
                    ax2.plot_trisurf([x1, x2, x3], [y1, y2, y3], [z1, z2, z3], linewidth=0.2, antialiased=True, alpha=0.2, color='green')

            ax2.scatter(P[:,0], P[:,1], P[:,2], s=10, alpha=0.6, color='blue')
            #ax.view_init(80, 50)

        if (base_folder is not None):
            fig.savefig('%s/%s_seq_%d.jpg' % (base_folder, name, c), bbox_inches='tight')
            plt.close(fig)
        else:
            plt.show()

        return None

    def complex_and_cocycles(self, P, S, Z, labels = True, xlim = [0, 1], base_folder = None, name = None, c = 0):

        fig = plt.figure(figsize=(10, 10))

        nP = P.shape[0]; dim = P.shape[1]

        eps = 0.01

        if (dim == 2):

            ax2 = fig.add_subplot(111)

            for s in S['complex']:
                s_ = s['simplex']
                if (s['dimension'] == 0):
                    pass
                if (s['dimension'] == 1):
                    [x0, y0] = P[s_[0],:]
                    [x1, y1] = P[s_[1],:]
                    ax2.plot([x0, x1], [y0, y1], color='gray', linewidth=0.4)
                if (s['dimension'] == 2):
                    idx = [s_[0], s_[1], s_[2]]
                    t1 = plt.Polygon(P[idx,:], alpha=0.3, color='lightgray')
                    plt.gca().add_patch(t1)

            ax2.scatter(P[:,0], P[:,1], s=14, alpha=0.6, color='gray')
            if (labels == True):
                for i, txt in enumerate(range(0, nP)):
                    ax2.annotate(txt, (P[i,0] + eps, P[i,1] + eps))
            ax2.set_xlim([xlim[0], xlim[1]])
            ax2.set_ylim([xlim[0], xlim[1]])

            colors = ['red', 'blue', 'green', 'yellow', 'cyan', 'black']
            for i, z in enumerate(Z):
                if (z['dimension'] == 1):
                    q = z['expression']
                    for tmp in q:
                        s_ = S['complex'][tmp[1]]['simplex']
                        [x0, y0] = P[s_[0],:]
                        [x1, y1] = P[s_[1],:]
                        ax2.plot([x0, x1], [y0, y1], alpha=1.0, color=colors[i % len(colors)], linewidth=2.0)

        if (base_folder is not None):
            fig.savefig('%s/%s_seq_%d.jpg' % (base_folder, name, c), bbox_inches='tight')
            plt.close(fig)
        else:
            plt.show()

        return None

    def point_cloud_cc(self, P, cc, labels = True, xlim = [0, 1], base_folder = None, name = None, c = 0):

        fig = plt.figure(figsize=(10, 10))

        nP = P.shape[0]; dim = P.shape[1]

        if (dim == 2):

            eps = 0.01
            ax1 = fig.add_subplot(111)
            for i in range(0, nP):
                cc_ = ((100 * cc[i]) % 10) / 10.0
                c = (0, cc_, cc_)
                circletmp = plt.Circle((P[i,0], P[i,1]), 0.04, color=c, alpha=0.8)
                ax1.add_artist(circletmp)
            if (labels == True):
                for i, txt in enumerate(range(0, nP)):
                    ax1.annotate(txt, (P[i,0] + eps, P[i,1] + eps))
            ax1.set_xlim([xlim[0], xlim[1]])
            ax1.set_ylim([xlim[0], xlim[1]])
            ax1.set_xlabel("x")
            ax1.set_ylabel("y")
            ax1.set_title("%d points" % (P.shape[0]))
            ax1.grid()

        if (base_folder is not None):
            fig.savefig('%s/%s_seq_%d.jpg' % (base_folder, name, c), bbox_inches='tight')
            plt.close(fig)
        else:
            plt.show()

        return None

    def persistent_barcode(self, ax, dim, r_max, q):
        eps = 0.1

        k0 = 0
        if ('intervals' in q):
            pi = q['intervals']
            if (dim in pi):
                for k0, tmp in enumerate(pi[dim]):
                    ax.plot([tmp[0], tmp[1]], [k0, k0], c='b', linewidth=4.0)

        k1 = 0
        if ('cocycles' in q):
            z = q['cocycles']
            if (dim in z):
                for k1, tmp in enumerate(z[dim]):
                    ax.plot([tmp[0], r_max], [k1+k0, k1+k0], c='r', linewidth=4.0)

        ax.set_xlim([0.0, r_max * (1 + eps)])
        ax.set_ylim([0.0, (k0+k1) * (1 + eps)])
        ax.set_xlabel("age (birth to death)")
        ax.set_ylabel("idx")
        ax.set_yticks(range(0, (k0 + k1) + 1))
        ax.set_yticklabels([''] * ((k0 + k1) + 1))
        ax.set_title("H^%d barcode" % (dim))
        ax.grid()

        return None

    def persistent_diagram(self, ax, dim, r_max, q, annotate = False, horizontal = False, weighted = False):
        eps = 0.1; eps_ = 0.001
        pi_size = 40
        z_size = 40

        if ('intervals' in q):
            pi = q['intervals']
            if (dim in pi):
                if (weighted == False):
                    if (pi[dim].shape[0] > 0):
                        if (horizontal == False):
                            ax.scatter(pi[dim][:,0], pi[dim][:,1], s=pi_size, c='b')
                            if (annotate == True):
                                for i, [x, y] in enumerate(pi[dim]):
                                    ax.annotate(i, (x + eps_, y + eps_))
                        else:
                            ax.scatter(pi[dim][:,0], pi[dim][:,1] - pi[dim][:,0], s=pi_size, c='b')
                            if (annotate == True):
                                for i, [x, y] in enumerate(pi[dim]):
                                    ax.annotate(i, (x + eps_, y - x + eps_))
                else:
                    if (pi[dim].shape[0] > 0):
                        for i, [x, y] in enumerate(pi[dim]):
                            y_ = y - x
                            circletmp = plt.Circle((x, y_), 0.1 * y_, color='b', alpha=0.4)
                            ax.add_artist(circletmp)
                            if (annotate == True):
                                ax.annotate(i, (x + eps_, y_ + eps_))
        if ('cocycles' in q):
            z = q['cocycles']
            if (dim in z):
                if (weighted == False):
                    if (z[dim].shape[0] > 0):
                        if (horizontal == False):
                            ax.scatter(z[dim][:,0], [r_max] * z[dim].shape[0], s=z_size, c='r')
                        else:
                            ax.scatter(z[dim][:,0], ([r_max] * z[dim].shape[0]) - z[dim][:,0], s=z_size, c='r')
                else:
                    if (z[dim].shape[0] > 0):
                        for i, [x] in enumerate(z[dim]):
                            y = r_max; y_ = y - x
                            circletmp = plt.Circle((x, y_), 0.1 * y_, color='r', alpha=0.4)
                            ax.add_artist(circletmp)
                            if (annotate == True):
                                ax.annotate(i, (x + eps_, y - x + eps_))

        ax.set_xlim([0, r_max * (1 + eps)])
        ax.set_ylim([0, r_max * (1 + eps)])

        if (horizontal == False):
            ax.plot([0.0, r_max * (1 + eps)], [0.0, r_max * (1 + eps)], c='r', linewidth=0.6)
            ax.set_ylabel("age (death)")
        else:
            ax.set_ylabel("age (lifespan)")

        ax.set_xlabel("age (birth)")
        ax.set_title("H^%d persistent diagram" % (dim))
        ax.grid()

        return None

    def persistent_distribution(self, ax, dim, q):
        if ('intervals' in q):
            pi = q['intervals']
            if (dim in pi):
                ax.hist(pi[dim][:,1], bins=40)
                ax.set_title("H^%d distribution" % (dim))
                ax.set_xlabel("age (death)")
                ax.set_ylabel("count")
                ax.grid()
        return None

    def bundle_0(self, dim_max, r_max, q, base_folder = None, name = None, c = 0):
        fig = plt.figure(figsize=(20, 20))
        i = 1
        for dim in range(0, dim_max):
            code = (dim_max * 100) + (3 * 10)

            self.persistent_barcode(fig.add_subplot(code + i), dim, r_max, q); i += 1
            self.persistent_diagram(fig.add_subplot(code + i), dim, r_max, q); i += 1
            self.persistent_distribution(fig.add_subplot(code + i), dim, q); i += 1

        if (base_folder is not None):
            fig.savefig('%s/%s_seq_%d.jpg' % (base_folder, name, c), bbox_inches='tight')
            plt.close(fig)
        else:
            plt.show()

    def bundle_1(self, dim_max, r_max, q, base_folder = None, name = None, c = 0):
        fig = plt.figure(figsize=(20, 20))
        i = 1
        for dim in range(0, dim_max):
            code = (dim_max * 100) + (2 * 10)

            self.persistent_diagram(fig.add_subplot(code + i), dim, r_max, q,
                annotate = False, horizontal = True, weighted = False); i += 1
            self.persistent_diagram(fig.add_subplot(code + i), dim, r_max, q,
                annotate = False, horizontal = True, weighted = True); i += 1

        if (base_folder is not None):
            fig.savefig('%s/%s_seq_%d.jpg' % (base_folder, name, c), bbox_inches='tight')
            plt.close(fig)
        else:
            plt.show()
