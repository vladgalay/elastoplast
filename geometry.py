import decider
import math
import iterator
import numpy as np
import tkinter as tk
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

coh = 0.3 # c
fi = math.pi*20/180 # deg
color = [[int(255*np.random.random()) for i in range(3)] for j in range(4)]

def to_hex(li):
    return '#{0:02X}{1:02X}{2:02X}'.format(li[0],li[1],li[2])

# def make_canvas(elms, pnts, U, load, scale=500):
#     """Рисунок в канвасе"""
#     cnvs.create_text(300, 15, text='Исходная схема', font=("Arial", "14"))
#     cnvs.create_text(300, 325, text='Деформированная схема', font=("Arial", "14"))
#     down = 300 # уровень низа построения сетки КЭ
#     # ось У направлена вниз, поэтому ниже для У использую знак "-"
#     for i in elms: # рисует квадратики и номера элементов
#         cnvs.create_rectangle(i.pnt[0].x*25+10, -i.pnt[0].y*25+down, \
#                             i.pnt[2].x*25+10, -i.pnt[2].y*25+down)
#         cnvs.create_text(i.pnt[0].x*25+25,-i.pnt[0].y*25+down-15,
#                         text=str(i.num-1), fill="red", \
#                         font=("Arial", "6"))
#     for i in pnts: # выводит номера узлов
#         cnvs.create_text(i.x*25+15,-i.y*25+down-5, \
#                         text=str(i.num-1), fill="blue", \
#                         font=("Arial", "6"))
#
#     for i in elms: # рисует деформированную схему
#         cnvs.create_polygon([i.pnt[0].x*25+U[2*(i.pnt[0].num-1)]*25*scale+10, \
#                         -i.pnt[0].y*25-U[2*i.pnt[0].num-1]*25*scale+down*2], \
#                         [i.pnt[1].x*25+U[2*(i.pnt[1].num-1)]*25*scale+10,\
#                         -i.pnt[1].y*25-U[2*i.pnt[1].num-1]*25*scale+down*2],\
#                         [i.pnt[2].x*25+U[2*(i.pnt[2].num-1)]*25*scale+10,\
#                         -i.pnt[2].y*25-U[2*i.pnt[2].num-1]*25*scale+down*2],\
#                         outline="black", fill='white')
#
#     for i in load:
#         if i[1] == 2 and i[2]<0:
#             Vload = -40
#         elif i[1] == 2 and i[2]>0:
#             Vload = 40
#         else:Vload = 0
#         if i[1] == 1 and i[2]<0:
#             Hload = 40
#         elif i[1] == 1 and i[2]>0:
#             Hload = -40
#         else:Hload = 0
#         cnvs.create_line(pnts[i[0]].x*25+U[2*i[0]]*25*scale+10+Hload,\
#                         -pnts[i[0]].y*25-U[2*i[0]+1]*25*scale+down*2+Vload,\
#                         pnts[i[0]].x*25+U[2*i[0]]*25*scale+10,\
#                         -pnts[i[0]].y*25-U[2*i[0]+1]*25*scale+down*2,\
#                         width=1,arrow='last')

def make_plot(elms, pnts):
    polygons = []
    codes = []
    for i in elms:  # рисует исходную схему
        polygons += [(i.pnt[0].x, i.pnt[0].y),
                     (i.pnt[1].x, i.pnt[1].y),
                     (i.pnt[2].x, i.pnt[2].y), (0, 0)]
        plt.text((i.pnt[0].x + i.pnt[1].x + i.pnt[2].x)/3-0.1,
                 (i.pnt[0].y + i.pnt[1].y + i.pnt[2].y)/3,
                 i.num, color='red', fontsize=6)
        codes += [Path.MOVETO] + [Path.LINETO] * 2 + [Path.CLOSEPOLY]
    for i in pnts:
        plt.text(i.x + 0.1,
                 i.y + 0.1,
                 i.num, color='blue', fontsize=6)

    polygons = np.array(polygons, float)
    path = Path(polygons, codes)
    pathpatch = PathPatch(path, facecolor='#ffffff')
    return pathpatch

def make_plot_plast(elms, pnts, U, load):
    # shows stress
    polygons_plast = []
    codes_plast = []
    polygons_nonplast = []
    codes_nonplast = []
    polygons_quad = []
    codes_quad = []
    for i in elms:  # рисует исходную схему
        # xx = math.sqrt((i.sx-i.sy)**2+4*i.tau**2) - (i.sx + i.sy+2*coh/math.tan(fi))*math.sin(fi)
        color = lambda xx: 'red' if xx >= 0 else 'blue'
        if len(i.pnt) == 3 and i.pl == True:
            polygons_plast += [(i.pnt[0].x, i.pnt[0].y),
                     (i.pnt[1].x, i.pnt[1].y),
                     (i.pnt[2].x, i.pnt[2].y), (0, 0)]

            codes_plast += [Path.MOVETO] + [Path.LINETO] * 2 + [Path.CLOSEPOLY]

        elif len(i.pnt) == 3 and i.pl != True:
            polygons_nonplast += [(i.pnt[0].x, i.pnt[0].y),
                     (i.pnt[1].x, i.pnt[1].y),
                     (i.pnt[2].x, i.pnt[2].y), (0, 0)]

            codes_nonplast += [Path.MOVETO] + [Path.LINETO] * 2 + [Path.CLOSEPOLY]

        if len(i.pnt) == 4:
            polygons_quad += [(i.pnt[0].x, i.pnt[0].y),
                     (i.pnt[1].x, i.pnt[1].y),
                     (i.pnt[2].x, i.pnt[2].y),
                     (i.pnt[3].x, i.pnt[3].y), (0, 0)]

            codes_quad += [Path.MOVETO] + [Path.LINETO] * 3 + [Path.CLOSEPOLY]

    polygons_plast = np.array(polygons_plast, float)
    print(polygons_plast)
    path_plast = Path(polygons_plast, codes_plast)
    pathpatch_PLAST = PathPatch(path_plast, facecolor='#fff000')
    polygons_nonplast = np.array(polygons_nonplast, float)
    path_nonplast = Path(polygons_nonplast, codes_nonplast)
    pathpatch_nonPLAST = PathPatch(path_nonplast, facecolor='#ffffff')
    polygons_quad = np.array(polygons_quad, float)
    path_quad = Path(polygons_quad, codes_quad)
    pathpatch_quad = PathPatch(path_quad, facecolor='#ffffff')
    return pathpatch_PLAST, pathpatch_nonPLAST, pathpatch_quad

def make_plot_def(elms, pnts, U, load, scale=100):
    polygons = []
    codes = []
    for i in elms: # рисует деформированную схему
        if len(i.pnt) == 3:
            polygons += [(i.pnt[0].x+U[2*(i.pnt[0].num-1)]*scale,\
                    i.pnt[0].y+U[2*i.pnt[0].num-1]*scale),\
                    (i.pnt[1].x+U[2*(i.pnt[1].num-1)]*scale,\
                    i.pnt[1].y+U[2*i.pnt[1].num-1]*scale),\
                    (i.pnt[2].x+U[2*(i.pnt[2].num-1)]*scale,\
                    i.pnt[2].y+U[2*i.pnt[2].num-1]*scale), (0, 0)]
            plt.text(i.pnt[0].x+U[2*i.pnt[0].num-1]*scale+0.1,
                    i.pnt[0].y+U[2*i.pnt[0].num-1]*scale+0.1,
                    i.num-1, color='red', fontsize=6)

        codes += [Path.MOVETO] + [Path.LINETO]*2 + [Path.CLOSEPOLY]
    polygons = np.array(polygons, float)
    path = Path(polygons, codes)
    pathpatch = PathPatch(path, facecolor='#ffffff')
    return pathpatch

def make_plot_z(elms, pnts, U, load, scale=100):
    polygons = []
    codes = []
    for i in elms:  # рисует деформированную схему
        if len(i.pnt) == 3:
            polygons += [(i.pnt[0].x + U[2 * (i.pnt[0].num-1)] * scale,
                      i.pnt[0].y + U[2 * i.pnt[0].num - 1] * scale),
                     (i.pnt[1].x + U[2 * (i.pnt[1].num-1)] * scale,
                      i.pnt[1].y + U[2 * i.pnt[1].num - 1] * scale),
                     (i.pnt[2].x + U[2 * (i.pnt[2].num-1)] * scale,
                      i.pnt[2].y + U[2 * i.pnt[2].num - 1] * scale), (0, 0)]
            codes += [Path.MOVETO] + [Path.LINETO] * 2 + [Path.CLOSEPOLY]
    for i in pnts:
        plt.text(i.x + U[2 * (i.num-1)] * scale + 0.1,
                 i.y + U[2 * i.num - 1] * scale + 0.1,
                 '{:.2f}'.format(U[2 * i.num - 1] * 1000), color='blue', fontsize=6)

    polygons = np.array(polygons, float)
    path = Path(polygons, codes)
    pathpatch = PathPatch(path, facecolor='#ffffff')
    return pathpatch

U = decider.matrix_U()
decider.matrix_tensor(U, decider.le)
U_nonlin = iterator.U

fig, ax = plt.subplots()
mpp=make_plot_plast(decider.le, decider.lp, U, decider.loads)
ax.add_patch(mpp[0])
ax.add_patch(mpp[1])
ax.add_patch(mpp[2])
# ax.add_patch(PathPatch(Path(np.array([(0,0),(1,0),(1,-1),(0,0)], float), [Path.MOVETO] + [Path.LINETO] * 2 + [Path.CLOSEPOLY]), facecolor='#000000'))
ax.set_title('Возникновение пластических деформации')
ax.autoscale_view()

# fig, ax = plt.subplots()
# ax.add_patch(make_plot_tensor(iterator.le, iterator.lp, U, decider.loads))
# # ax.add_patch(PathPatch(Path(np.array([(0,0),(1,0),(1,-1),(0,0)], float), [Path.MOVETO] + [Path.LINETO] * 2 + [Path.CLOSEPOLY]), facecolor='#000000'))
# ax.set_title('Пласт деформ nonlin')
# ax.autoscale_view()

fig, ax = plt.subplots()
ax.add_patch(make_plot(decider.le, decider.lp))
ax.set_title('Деформированная схема')
ax.autoscale_view()

fig, ax = plt.subplots()
ax.add_patch(make_plot_z(decider.le, decider.lp, U, decider.loads))
ax.set_title('Перемещения по Z(мм)\nЛинейная задача')
ax.autoscale_view()

fig, ax = plt.subplots()
ax.add_patch(make_plot_z(decider.le, decider.lp, U_nonlin, decider.loads, scale=10))
ax.set_title('Перемещения по Z(мм)\nНелинейная задача')
ax.autoscale_view()
#
fig, ax = plt.subplots()
y = iterator.deform
x = np.linspace(0.01, 1, 100)
plt.plot(x, y, 'k')

plt.show()

# root = tk.Tk()
# cnvs = tk.Canvas(root, width=800, height=800)
# make_canvas(decider.le, decider.lp, U, decider.loads)
# cnvs.pack()
# root.mainloop()
if __name__ == '__main__':
    pass# main()
