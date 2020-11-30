import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

nx = 1001;
ny = 1001;

x = np.linspace(0, 2, nx)
y = np.linspace(0, 2, ny)

X, Y = np.meshgrid(x, y)


for i in range(0,6000,60):
	u = np.genfromtxt('output/u-%i.out'%i,delimiter=',')
	v = np.genfromtxt('output/v-%i.out'%i,delimiter=',')

	fig = plt.figure(figsize=(11, 7), dpi=100)
	ax = fig.gca(projection='3d')
	ax.plot_surface(X, Y, u, cmap=cm.viridis, rstride=1, cstride=1)
	ax.plot_surface(X, Y, v, cmap=cm.viridis, rstride=1, cstride=1)
	ax.set_xlabel('$x$')
	ax.set_ylabel('$y$');
	ax.set_zlim(0.9, 2.1)
	plt.savefig('output/Graph_time_%i.png' % i)

