import numpy as np

def dd(mu1, mu2, r):
    r1 = np.linalg.norm(r)
    return (mu1@mu2)/r1**3 - 3*(mu1@r)*(mu2@r)/r1**5

def dds(mu1, mu2s, rs):
    r1s = np.linalg.norm(rs, axis=1)
    print(np.sum(mu2s*rs, axis=1))
    return (mu1@mu2s.T)/r1s**3 - 3*(mu1@rs.T)*(np.sum(mu2s*rs, axis=1))/r1s**5


import sys
from lattices import diamond, icerule, bodycenter, facecenter, simplecubic
import networkx as nx
from pairlist import pairs_iter

lattice_type, lattice_size = sys.argv[1:]

N = int(lattice_size)
if lattice_type == "dia":
    rpos, cell = diamond(N)
    ddd = (2/3)**0.5
elif lattice_type == "fcc":
    rpos, cell = facecenter(N)
    ddd = 1/2
elif lattice_type == "bcc":
    rpos, cell = bodycenter(N)
    ddd = 1/3**0.5
elif lattice_type == "scl":
    rpos, cell = simplecubic(N)
    # dipole-dipole distance in the lattice of which node-node length is unity.
    ddd = 2**0.5/2


g = nx.Graph([[i,j] for i,j in pairs_iter(rpos, 1.1, cell, distance=False)])
dg = icerule(g)
print(dg.edges())
pairs = dg.edges()



# この時点では、格子点間距離を1としている。しかし、1にしたいのはスピン間距離。
# 最近接スピン間距離を求める必要がある。

cell /= ddd


mus = []
rs  = []
for i,j in pairs:
    d = rpos[j] - rpos[i]
    d -= np.floor(d+0.5)
    mus.append(d@cell)
    rs.append(rpos[i]+d/2)
mus = np.array(mus)
rs  = np.array(rs)


import random
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import matplotlib.cm as cm


fig = plt.figure()

X = np.linspace(0, cell[0,0]/2, 1000)

Ys = []

Nsample=500
for center in random.sample(range(len(rpos)), Nsample):
    muc = mus[center]
    rc  = rs[center]

    d = rs - rc
    d -= np.floor(d+0.5)
    d = d @ cell
    Es = dds(muc, mus, d)
    rr = np.linalg.norm(d, axis=1)
    order = np.argsort(rr)
    Es[center] = 0.0
    Esc = np.cumsum(Es[order])
    rrc = rr[order]
    print(rrc[0])

    Escp = interp1d(rrc, Esc, kind='previous')
    Ys.append(Escp(X))

Ys = np.array(Ys)

# 距離4でのエネルギーの値でソートする。

e25 = np.zeros(Nsample)
for i, e in enumerate(Ys):
    e25[i] = e[X<2.5][-1]
order = np.argsort(e25)
redro = np.zeros(Nsample)
for i in range(Nsample):
    redro[order[i]] = i
for i in range(Nsample):
    j = redro[i]
    plt.plot(X, Ys[i], color=cm.coolwarm(j/Nsample), lw=0.5)

plt.xlim(0,cell[0,0]/2)
plt.xlabel(r"$r$", fontsize=18)
plt.ylabel(r"$4\pi\epsilon_0 E(r)$", fontsize=18)
plt.tick_params(labelsize=14)

fig.savefig(f"{lattice_type}{N}.raw.pdf", bbox_inches="tight")
fig.savefig(f"{lattice_type}{N}.raw.png", bbox_inches="tight")

fig = plt.figure()

plt.plot(X, np.std(Ys, axis=0))
plt.xlim(1,cell[0,0]/2)
plt.ylim(5e-2,None)
plt.xscale("log", base=10)
plt.yscale("log", base=10)
plt.xlabel(r"$r$", fontsize=18)
plt.ylabel(r"$4\pi\epsilon_0\sigma_E(r)$", fontsize=18)
plt.tick_params(labelsize=14)

fig.savefig(f"{lattice_type}{N}.sd.pdf", bbox_inches="tight")
fig.savefig(f"{lattice_type}{N}.sd.png", bbox_inches="tight")
