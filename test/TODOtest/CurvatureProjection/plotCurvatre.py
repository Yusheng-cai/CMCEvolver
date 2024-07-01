import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def read_data(file:str):
    with open(file, "r") as f:
        lines = f.readlines()
        lines = [l.rstrip("\n").lstrip().split() for l in lines]
        lines = [[float(r) for r in l] for l in lines]
        lines = np.array(lines)
    return lines

if __name__ == "__main__":
    read_file = "c.out" 
    data = read_data(read_file)

    c = np.zeros((400,400))
    for d in data:
        index1, index2 = int(d[0]), int(d[1])
        c[index1,index2] = d[-1]

    fig = plt.figure(dpi=300)
    ax  = fig.add_subplot(111)
    p = cm.ScalarMappable(cmap=cm.jet)
    p.set_array([0,0.5])
    p.set_clim(0,0.5)
    colors = p.to_rgba(c.T)
    ax.imshow(colors, interpolation="nearest")
    fig.colorbar(p)
    plt.savefig("test.png")

    # plot normal map
    norm = read_data("normalMap.out")
    normMap = np.zeros((1000,1000,4))
    for d in norm:
        index1,index2 = int(d[0]), int(d[1])
        normMap[index2,index1] = np.array([d[2],d[3],d[4],1])

    fig = plt.figure(dpi=300)
    ax  = fig.add_subplot(111)
    p   = ax.imshow(normMap, interpolation="none")
    fig.colorbar(p)
    plt.savefig("a.png")
