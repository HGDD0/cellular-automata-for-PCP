import pickle
import random
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

class Cell(object):

    def __init__(self, position, polarity):
        self.position=position
        self.polarity=polarity


    def compute_polarity(self,neighbor):

        arrays = []
        for cell in neighbor:
            #sub_vct=cell.polarity

            sub=np.array(self.position)-np.array(cell.position)
            sub=sub/np.linalg.norm(sub)
            cos_fi=np.dot(sub,cell.polarity)/(np.linalg.norm(cell.polarity))
            #if cos_fi < 0:
            #    cos_fi = 0
            sub_vct=sub*np.linalg.norm(cell.polarity)*cos_fi

            arrays.append(sub_vct)
        #方法一：直接向量求和
        # 方法二：只限于距离一，求出垂直边界量之后求和

        polarity = self.polarity + 0.2*sum(arrays)
        mo=np.linalg.norm(polarity)
        if mo != 0 and mo > 1 :
            polarity = polarity / mo
        return polarity

    def refresh(self,neighbor):
        new_polarity= self.compute_polarity(neighbor)
        return Cell(self.position,new_polarity)

x=30
y=30
z=30
flat=29
aliable_cells= np.empty([x,y,z], dtype = Cell)
aliable_position=[]

X=[]
Y=[]
Z=[]
U=[]
V=[]
W=[]

a=np.array([0,1,-1])
b=np.array([-1,0,1])
c=np.array([1,-1,0])
rand_scal=100
dis=1
#设定常量

mis=0
change=True

class Playground(object):

    def __init__(self,x,y,z):
        for m in range(x):
            for n in range(y):
                for l in range(z):
                    if m + n + l == flat:
                        #if (pow((m - int(x / 3)),2) + pow((n - int(y / 3)),2) >= 8):
                        #if (m > int(x/2-1) or n> int(y/2-1) ):
                        #    continue
                        X.append(m)
                        Y.append(n)
                        Z.append(l)
                        vctor = random.randint(-rand_scal, rand_scal + 1) * a + random.randint(-rand_scal,rand_scal + 1) * b + random.randint(-rand_scal, rand_scal + 1) * c
                        mo = np.linalg.norm(vctor)
                        if mo != 0:
                            vctor = vctor / mo
                        U.append(vctor[0])
                        V.append(vctor[1])
                        W.append(vctor[2])
                        aliable_cells[m, n, l] = Cell([m, n, l], vctor)
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
                if cell is aliable_cells[m,n,l] or aliable_cells[m,n,l] is None:
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
        X.clear()
        Y.clear()
        Z.clear()
        U.clear()
        V.clear()
        W.clear()

        global change
        change=False
        global mis
        mis = 0

        for pos in aliable_position:
            #if cell is None:
            #    continue
            #else:
            cell=aliable_cells[pos[0],pos[1],pos[2]]
            newcell=cell.refresh(self.get_neighbor(cell))
            aliable_cells[cell.position[0],cell.position[1],cell.position[2]]=newcell
            X.append(newcell.position[0])
            Y.append(newcell.position[1])
            Z.append(newcell.position[2])
            vctor=newcell.polarity
            U.append(vctor[0])
            V.append(vctor[1])
            W.append(vctor[2])
            #if((change is False) and not((newcell.polarity==cell.polarity).all())):
            #    change=True
            newmis=newcell.polarity-cell.polarity
            mis+=np.linalg.norm(newmis)
        self.step += 1

    def hulahula(self, rounds):
        """更新状态并画图
        Parameters
        ----------
        rounds : 更新的轮数
        """

        plt.ion()
        fig = plt.figure(figsize=(12,7),dpi=140)
        ax = fig.gca(projection='3d')
        ax.view_init(45, 45)
        fig2 = plt.figure(figsize=(8, 5))
        bx = fig2.gca()

        #创建绘布
        global change
        global mis
        for i in range(rounds):
            ax.set_zlabel("Step: " + str(self.step))
            ax.set_xlabel("Change: " + str(mis)+' '+ str(change))
            ax.scatter(X, Y, Z, s=20, c="black")
            ax.quiver(X, Y, Z, U, V, W, length=0.7, arrow_length_ratio=.3, pivot='tail', normalize=False)
            bx.scatter(i,mis, s=20, c="black")
            plt.pause(0.01)
            self.update_state()

            #filename = "./Test parameter/"+ "Default" + str(i) + ".png"
            #plt.savefig(filename)
            ax.clear()

        plt.ioff()

def save_value (Playground,X,Y,Z,U,V,W,aliable_cells,aliable_position):
    f = open("./DD/"+ "Playground", 'wb')
    pickle.dump(Playground, f)
    f = open("./DD/" + "X", 'wb')
    pickle.dump(X, f)
    f = open("./DD/" + "Y", 'wb')
    pickle.dump(Y, f)
    f = open("./DD/" + "Z", 'wb')
    pickle.dump(Z, f)
    f = open("./DD/" + "U", 'wb')
    pickle.dump(U, f)
    f = open("./DD/" + "V", 'wb')
    pickle.dump(V, f)
    f = open("./DD/" + "W", 'wb')
    pickle.dump(W, f)
    f = open("./DD/" + "aliable_cells", 'wb')
    pickle.dump(aliable_cells, f)
    f = open("./DD/" + "aliable_position", 'wb')
    pickle.dump(aliable_position, f)
    f.close()

def load_value ():
    f = open("./DD/"+ "Playground", 'rb')
    game=pickle.load(f)
    f = open("./DD/" + "X", 'rb')
    X=pickle.load(f)
    f = open("./DD/" + "Y", 'rb')
    Y=pickle.load(f)
    f = open("./DD/" + "Z", 'rb')
    Z=pickle.load(f)
    f = open("./DD/" + "U", 'rb')
    U=pickle.load(f)
    f = open("./DD/" + "V", 'rb')
    V=pickle.load(f)
    f = open("./DD/" + "W", 'rb')
    W=pickle.load(f)
    f = open("./DD/" + "aliable_cells", 'rb')
    aliable_cells =pickle.load(f)
    f = open("./DD/" + "aliable_position", 'rb')
    aliable_position = pickle.load(f)
    return game,X,Y,Z,U,V,W,aliable_cells,aliable_position
    f.close()

def add_bias ():
    X.clear()
    Y.clear()
    Z.clear()
    U.clear()
    V.clear()
    W.clear()
    for cell in aliable_cells.flat:
        if cell is None:
            continue
        else:
            X.append(cell.position[0])
            Y.append(cell.position[1])
            Z.append(cell.position[2])
            vctor = cell.polarity + np.array([0, -1, 1])/1.414
            mo = np.sqrt(np.sum(vctor ** 2))
            if mo != 0:
                vctor = vctor / mo
            cell.polarity=vctor
            U.append(vctor[0])
            V.append(vctor[1])
            W.append(vctor[2])
    return
    # 添加全局偏转

def add_line ():
    X.clear()
    Y.clear()
    Z.clear()
    U.clear()
    V.clear()
    W.clear()
    for cell in aliable_cells.flat:
        if cell is None:
            continue
        else:
            if cell.position[2]>14:
                X.append(cell.position[0])
                Y.append(cell.position[1])
                Z.append(cell.position[2])
                vctor = np.array([1,1,-2])
                mo = np.sqrt(6)
                vctor = vctor / mo
                cell.polarity=vctor
                U.append(vctor[0])
                V.append(vctor[1])
                W.append(vctor[2])
            else:
                X.append(cell.position[0])
                Y.append(cell.position[1])
                Z.append(cell.position[2])
                vctor = cell.polarity
                U.append(vctor[0])
                V.append(vctor[1])
                W.append(vctor[2])
    return


if __name__ == '__main__':

    game = Playground(x,y,z)
    save_value(game,X,Y,Z,U,V,W,aliable_cells,aliable_position)
    #game,X,Y,Z,U,V,W,aliable_cells,aliable_position=load_value()
    game.hulahula(50000)