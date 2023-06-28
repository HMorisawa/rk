from re import X
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

data=np.genfromtxt('data_r0.csv', delimiter=',')

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')


x=np.zeros(len(data))
y=np.zeros(len(data))
z=np.zeros(len(data))
# x, y, z成分のデータの作成
for i in range(len(data)):
  x[i]=data[i][1]
  y[i]=data[i][2]
  z[i]=data[i][3]

ax.plot(x, y, z, color='blue')
#ax.scatter(x, y, color='blue')
ax.set_box_aspect((2,1,1))

plt.show()
