import pickle
import random
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from numpy import cross,eye,dot
from scipy.linalg import expm,norm
from pandas.core.frame import DataFrame
import pandas as pd
from matplotlib import gridspec
import matplotlib.colors

class Cell(object):

    def __init__(self, position, A_fi, A_mo,P_fi,P_mo):
        self.position=position
        self.A_polarity=[A_fi,A_mo]
        self.P_polarity = [P_fi,P_mo]


    def compute_polarity(self,neighbor):
        u=0.95
        len=10
        dAfi=0
        dPfi=0
        dAmo=0
        dPmo=0
        for cell in neighbor:
            Aa=(self.P_polarity[1]/2)*np.cos(self.A_polarity[0]-self.P_polarity[0])-u
            Pa=(self.A_polarity[1]/2)*np.cos(self.P_polarity[0]-self.A_polarity[0])-u
            dAfi += cell.P_polarity[1] * np.sin(np.pi + cell.P_polarity[0] - self.A_polarity[0])
            dPfi += cell.A_polarity[1] * np.sin(np.pi + cell.A_polarity[0] - self.P_polarity[0])
            dAmo += self.A_polarity[1] * (1 - self.A_polarity[1]) * (self.A_polarity[1] - 1 / 2 - Aa)
            dPmo += self.P_polarity[1] * (1 - self.P_polarity[1]) * (self.P_polarity[1] - 1 / 2 - Pa)

        A_fi = self.A_polarity[0] + dAfi/len
        P_fi = self.P_polarity[0] + dPfi/len
        A_mo = self.A_polarity[1] + dAmo/len
        P_mo = self.P_polarity[1] + dPmo/len

        if A_fi >= 2 * np.pi:
            A_fi = A_fi % (2 * np.pi)
        if A_fi < 0:
            A_fi = A_fi % (-2 * np.pi) + 2 * np.pi
        if P_fi >= 2 * np.pi:
            P_fi = P_fi % (2 * np.pi)
        if P_fi < 0:
            P_fi = P_fi % (-2 * np.pi) + 2 * np.pi
        if A_mo >= 1:
            A_mo = 1
        if A_mo < 0:
            A_mo = 0
        if P_mo >= 1:
            P_mo = 1
        if P_mo < 0:
            P_mo = 0

        return A_fi, A_mo,P_fi,P_mo

    def refresh(self,neighbor):
        A_fi, A_mo,P_fi,P_mo= self.compute_polarity(neighbor)
        self.A_polarity=[A_fi, A_mo]
        self.P_polarity=[P_fi, P_mo]

x=60
y=60
z=60
flat=59
aliable_cells= np.empty([x,y,z], dtype = Cell)
aliable_position=[]

#XYZ 点坐标用于绘图 UVW向量坐标用于绘图
X=[]
Y=[]
Z=[]
Ua=[]
Va=[]
Wa=[]
Up=[]
Vp=[]
Wp=[]
C=[]
Check=[]
windrose=[]

#切面上的三个方向的平面矢量，取c(正右)为0°
a=np.array([0,1,-1])
b=np.array([-1,0,1])
c=np.array([1,-1,0])
rand_scal=100
dis=1
#设定常量

def rotation(axis, theta):
    return expm(cross(eye(3), axis / norm(axis) * theta))

class Playground(object):

    def __init__(self,x,y,z):
        for m in range(x):
            for n in range(y):
                for l in range(z):
                    if m + n + l == flat:
                        #确定参与的点
                        #if (pow((m - int(x / 3)),2) + pow((n - int(y / 3)),2) >= 8):
                        #if (m > int(x/2-1) or n> int(y/2-1) ):
                        #    continue
                        X.append(m)
                        Y.append(n)
                        Z.append(l)
                        A_fi=random.uniform(0,2*np.pi)
                        A_mo=random.random()
                        P_fi=random.uniform(0,2*np.pi)
                        P_mo=random.random()

                        v, axis, theta = [1, -1, 0], [1, 1, 1], A_fi
                        M0 = rotation(axis, theta)
                        vctor = dot(M0, v)
                        vctor = vctor / np.linalg.norm(vctor) * A_mo
                        Ua.append(vctor[0])
                        Va.append(vctor[1])
                        Wa.append(vctor[2])

                        v, axis, theta = [1, -1, 0], [1, 1, 1], P_fi
                        M0 = rotation(axis, theta)
                        vctor = dot(M0, v)
                        vctor = vctor / np.linalg.norm(vctor) * P_mo
                        Up.append(vctor[0])
                        Vp.append(vctor[1])
                        Wp.append(vctor[2])

                        C.append(A_fi)

                        aliable_cells[m, n, l] = Cell([m, n, l], A_fi, A_mo,P_fi,P_mo)
                        aliable_position.append([m,n,l])
        #初始化细胞以及随机极性

        self.step = 0
        print("----Initiation done----")

    def distance(self,Cell1,Cell2):
        return ((abs(Cell1.position[0] - Cell2.position[0])+abs(Cell1.position[1] - Cell2.position[1])+abs(Cell1.position[2] - Cell2.position[2]))/2)

    def get_neighbor(self,cell):
        neighborcells=[]

        for m in range(max(cell.position[0]-dis,0),min(cell.position[0]+dis+1,x)):
            for n in range(max(cell.position[1]-dis, cell.position[1]+cell.position[0] -m - dis, 0),min(cell.position[1]+dis+1, cell.position[1]+cell.position[0] -m + dis+1, y, flat-m+1)):
                l=flat-m-n
                if aliable_cells[m,n,l] is None:
                    continue
                neighborcells.append(aliable_cells[m,n,l])
        '''
        for cells in aliable_cells.flat:
            if cells is None:
                continue
            else:
                distance=self.distance(cell, cells)
                if distance <= dis and distance !=0 :
                    neighborcells.append(cells)
        '''
        return neighborcells



    def update_state(self):
        """更新一次状态"""
        #X.clear()
        #Y.clear()
        #Z.clear()
        C.clear()
        #Ua.clear()
        #Va.clear()
        #Wa.clear()
        #Up.clear()
        #Vp.clear()
        #Wp.clear()

        for pos in aliable_position:
            #if cell is None:
            #    continue
            #else:
            cell=aliable_cells[pos[0],pos[1],pos[2]]
            cell.refresh(self.get_neighbor(cell))
            C.append(cell.A_polarity[0])
            # v, axis, theta = [1,-1,0], [1,1,1], cell.A_polarity[0]
            # M0 = rotation(axis, theta)
            # vctor=dot(M0, v)
            # vctor=vctor/np.linalg.norm(vctor)*cell.A_polarity[1]
            # Ua.append(vctor[0])
            # Va.append(vctor[1])
            # Wa.append(vctor[2])

            #v, axis, theta = [1, -1, 0], [1, 1, 1], cell.P_polarity[0]
            #M0 = rotation(axis, theta)
            #vctor = dot(M0, v)
            #vctor = vctor / np.linalg.norm(vctor) * cell.P_polarity[1]
            #Up.append(vctor[0])
            #Vp.append(vctor[1])
            #Wp.append(vctor[2])

        self.step += 1

    def hulahula(self, rounds):
        """更新状态并画图
        Parameters
        ----------
        rounds : 更新的轮数
        """

        plt.ion()
        fig = plt.figure(figsize=(12,7),dpi=150)
        gs = gridspec.GridSpec(2, 2, width_ratios=[5, 1])
        ax = plt.subplot(gs[:,0], projection='3d')
        ax.view_init(45, 45)
        #fig2 = plt.figure(figsize=(10, 6),dpi=150)
        #bx = fig2.gca(projection='3d')
        #bx.view_init(45, 45)

        surf = ax.scatter(X, Y, Z, s=20, c=C, cmap='hsv', alpha=1)
        #fig.colorbar(surf, shrink=1, aspect=20, ticks=[0,np.pi/2,np.pi,3*np.pi/2,2*np.pi])

        m=60
        cx = plt.subplot(gs[0,-1], projection='polar')
        dx = plt.subplot(gs[1,-1], projection='polar')
        #创建图例
        theta = np.linspace(0.0, 2 * np.pi, m, endpoint=False)
        radius = np.ones(m)
        width = np.pi * 2 / m
        cmap = plt.cm.hsv
        norm = matplotlib.colors.Normalize(vmin=0, vmax=2*np.pi)
        dx.bar(theta, radius, width=width, bottom=0.0, color=cmap(norm(theta)), alpha=1)
        dx.set_rticks([])
        #创建绘布

        for i in range(rounds):
            ax.set_zlabel("Step: " + str(self.step))
            #ax.set_xlabel("Change: " + str(mis)+' '+ str(change))
            ax.scatter(X, Y, Z, s=20, c=C, cmap='hsv', alpha=1,vmin=0, vmax=2*np.pi)
            #bx.scatter(X, Y, Z, s=20, c=C, cmap='gray')
            #ax.quiver(X, Y, Z, Ua, Va, Wa, length=0.7, arrow_length_ratio=.3, pivot='tail', normalize=False)
            #ax.quiver(X, Y, Z, Up, Vp, Wp, length=0.7, arrow_length_ratio=.3, pivot='tail', normalize=False)
            #bx.scatter(i,mis, s=20, c="black")

            #绘制风向玫瑰
            n = 12
            D=[c/(2*np.pi/n) for c in C]
            a = [int(d) for d in D]
            data = pd.value_counts(a)
            dict=data.to_dict()
            for k in range(0,n):
                if k in dict.keys():
                    windrose.append(dict[k])
                else:
                    windrose.append(0)


            theta = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
            radius = np.array(windrose)
            colors = cmap(norm(theta))
            width = np.pi*2 /n
            cx.bar(theta, radius,  width=width, bottom=0.0, color=colors, alpha=0.5)

            plt.pause(0.01)
            self.update_state()

            filename = "./CA_HEX_2/"+ "Default" + str(i) + ".png"
            plt.savefig(filename)
            
            ax.clear()
            # bx.clear()
            cx.clear()
            windrose.clear()

        plt.ioff()

if __name__ == '__main__':

    game = Playground(x,y,z)
    game.hulahula(5000)