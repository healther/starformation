# -*- coding: utf-8 -*-
#test commentary
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.mlab import griddata
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

f = open('oldout/__expected_number', 'r')
headers = f.readline().strip().split(',')[1:]
data = np.loadtxt(f)
f.close()


fig = plt.figure()
ax = fig.gca(projection='3d')
avs = np.linspace(10.625, 38.75, 6)
aperas = np.linspace(15000, 50000, 7)
ages = np.linspace(500000, 2000000, 5)
numbers = data[:,3].reshape(6, 7, 5)
numbers = np.roll(numbers,2,2)

X, Y = np.meshgrid(avs, ages)
Z = np.transpose(numbers[:,3,:])/559.
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm)
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#fig.colorbar(surf, shrink=0.5, aspect=5)
fig.colorbar(surf)


x = data[:,0].reshape(6, 7, 5)[:,3,:].reshape(30)
y = data[:,2].reshape(6, 7, 5)[:,3,:].reshape(30)
z = data[:,3].reshape(6, 7, 5)[:,3,:].reshape(30)/559.
ax.scatter(x,y,z)
ax.set_xlabel('av')
ax.set_ylabel('age')
ax.set_zlabel('starformation rate in 0.001 M_sun/year')

plt.savefig('av-age.svg')

fig2 = plt.figure()
ax2 = fig2.gca(projection='3d')


X2, Y2 = np.meshgrid(ages, aperas)
Z2 = numbers[0,:,:]/559.
surf2 = ax2.plot_surface(X2, Y2, Z2, rstride=1, cstride=1, cmap=cm.coolwarm)
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#fig.colorbar(surf, shrink=0.5, aspect=5)
fig2.colorbar(surf2)


x2 = data[:,2].reshape(6, 7, 5)[0,:,:].reshape(35)
y2 = data[:,1].reshape(6, 7, 5)[0,:,:].reshape(35)
z2 = data[:,3].reshape(6, 7, 5)[0,:,:].reshape(35)/559.
ax2.scatter(x2,y2,z2)
ax2.set_xlabel('age')
ax2.set_ylabel('apera')
ax2.set_zlabel('starformation rate in 0.001 M_sun/year')



plt.savefig('age-apera.svg')
