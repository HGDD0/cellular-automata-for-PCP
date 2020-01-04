import random
import numpy as np
import matplotlib.pyplot as plt

m=100
n=100
#m行n列

cellsmap = np.zeros((m,n))

aliable_cells=[]
new_cells=[]

class Playground(object):

    def __init__(self, cells_shape,ini_cell):
        """
        Parameters
        ----------
        cells_shape : 一个元组，表示画布的大小。
        """

        # 矩阵的四周不参与运算
        self.step = 0
        cellsmap[ini_cell.position[0],ini_cell.position[1]]=1

    def get_neighbor(self,cell):
        return cellsmap[cell.position[0]-1:cell.position[0]+2,cell.position[1]-1:cell.position[1]+2]

    def update_state(self):
        """更新一次状态"""
        print(str(self.step)+"...")

        if self.step<4:
            self.step += 1
            return
        #建立基础6细胞起源

        for cells in aliable_cells:
            a=cells.divide(self.get_neighbor(cells))
            if a == 0:
                new_cells.append(cells)
                continue
            del cells
            new_cells.append(a[0])
            new_cells.append(a[1])
        aliable_cells.clear()
        aliable_cells.extend(new_cells)
        new_cells.clear()

        self.step += 1

    def hulahula(self, rounds):
        """更新状态并画图
        Parameters
        ----------
        rounds : 更新的轮数
        """
        plt.ion()
        for _ in range(rounds):
            plt.title('Iter :{}'.format(self.step))
            plt.imshow(cellsmap)
            self.update_state()
            plt.pause(0.2)
        plt.ioff()



class Cell(object):

    def __init__(self, position, polarity,pop1,lineage):
        self.position=position
        self.polarity=polarity
        self.pop1 = pop1
        self.lineage= lineage
        #细胞的位置，极性方向，POP1水平

    def divide(self,neighbor):
        #细胞分裂，已知分裂新细胞的方向，根据原细胞的信息判断新细胞的类型以及位置

        def neighbor_analyse(neighbor):
            num=0
            empty=[]
            for cells in neighbor:
                for cell in cells:
                    num+=1
                    if cell == 0:
                        empty.append(num)
            if len(empty):
                a = random.sample(empty, 1)[0]
                print(a)
                return (a)
            else:
                return 0


        def getdirction(num):
            dirction = {
                1: [self.position[0]-1, self.position[1]-1],
                2: [self.position[0]-1, self.position[1]],
                3: [self.position[0]-1, self.position[1]+1],
                4: [self.position[0], self.position[1]-1],
                6: [self.position[0], self.position[1]+1],
                7: [self.position[0]+1, self.position[1]-1],
                8: [self.position[0]+1, self.position[1]],
                9: [self.position[0]+1, self.position[1]+1]
            }
            return (dirction.get(num))

        def usual_divide(dirction):
            if dirction == 0:
                return 0
            daughtercells=[]
            daughtercell1=Cell(self.position,self.polarity,self.pop1,self.lineage)
            daughtercell2=Cell(self.position,self.polarity,self.pop1,self.lineage)
            daughtercell2.position=getdirction(dirction)
            cellsmap[daughtercell2.position[0],daughtercell2.position[1]]=1
            print(daughtercell1.position)
            print(daughtercell2.position)
            daughtercells.append(daughtercell1)
            daughtercells.append(daughtercell2)
            return(daughtercells)

        lineage_tree = {
            "P0": ["AB", "P1"],

            "AB": ["ABa", "ABp"],
            "P1": ["EMS", "P2"],

            "ABa": ["al", "ar"],
            "ABp": ["pl", "Pr"],
            "EMS": ["MS", "E"],
            "P2": ["C", "P3"],

            "al": ["ala", "alp"],
            "ar": ["ara", "arp"],
            "Pl": ["pla", "plp"],
            "Pr": ["pra", "prp"],
            "MS": ["MSa", "MSp"],
            "E": ["Ea", "Ep"],
            "C": ["Ca", "Cp"],
            "P3": ["D", "P4"],
        }

        if self.lineage in lineage_tree.keys():
            return
        else:
            return (usual_divide(neighbor_analyse(neighbor)))


if __name__ == '__main__':

    ini_cell=Cell([int(m/2),int(n/2)], 5, -1, "0")
    aliable_cells.append(ini_cell)
    game = Playground((m, n),ini_cell)
    game.hulahula(80)

'''
        real_width = cells_shape[0] - 2
        real_height = cells_shape[1] - 2

        self.cells[1:-1, 1:-1] = np.random.randint(2, size=(real_width, real_height))
        self.timer = 0
        self.mask = np.ones(9)
        self.mask[4] = 0
'''
'''
    def update_state(self):
        """更新一次状态"""
        buf = np.zeros(self.cells.shape)
        cells = self.cells
        for i in range(1, cells.shape[0] - 1):
            for j in range(1, cells.shape[0] - 1):
                # 计算该细胞周围的存活细胞数
                neighbor = cells[i - 1:i + 2, j - 1:j + 2].reshape((-1,))
                neighbor_num = np.convolve(self.mask, neighbor, 'valid')[0]
                if neighbor_num == 3:
                    buf[i, j] = 1
                elif neighbor_num == 2:
                    buf[i, j] = cells[i, j]
                else:
                    buf[i, j] = 0
        self.cells = buf
        self.timer += 1

    def plot_state(self):
        """画出当前的状态"""
        plt.title('Iter :{}'.format(self.timer))
        plt.imshow(self.cells)
        plt.show()

    def update_and_plot(self, n_iter):
        """更新状态并画图
        Parameters
        ----------
        n_iter : 更新的轮数
        """
        plt.ion()
        for _ in range(n_iter):
            plt.title('Iter :{}'.format(self.timer))
            plt.imshow(self.cells)
            self.update_state()
            plt.pause(2)
        plt.ioff()
'''