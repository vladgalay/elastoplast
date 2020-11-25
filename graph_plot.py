import numpy as np
import math
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import matplotlib.pyplot as plt


def make_plot_plast(elms, scale=0):
    '''Рисует исходную схему'''
    polygons_plast = []
    codes_plast = []
    polygons_nonplast = []
    codes_nonplast = []
    polygons_quad = []
    codes_quad = []
    for i in elms:  # рисует исходную схему
        if len(i.pnt) == 3 and i.pl is True:
            polygons_plast += [
                (i.pnt[0].x + i.pnt[0].tot_displ_x * scale,
                 i.pnt[0].y + i.pnt[0].tot_displ_y * scale),
                (i.pnt[1].x + i.pnt[1].tot_displ_x * scale,
                 i.pnt[1].y + i.pnt[1].tot_displ_y * scale),
                (i.pnt[2].x + i.pnt[2].tot_displ_x * scale,
                 i.pnt[2].y + i.pnt[2].tot_displ_y * scale),
                (0, 0)]
            codes_plast += [Path.MOVETO] + [Path.LINETO]*2 + [Path.CLOSEPOLY]

        elif len(i.pnt) == 3 and i.pl is not True:
            polygons_nonplast += [
                (i.pnt[0].x + i.pnt[0].tot_displ_x * scale,
                 i.pnt[0].y + i.pnt[0].tot_displ_y * scale),
                (i.pnt[1].x + i.pnt[1].tot_displ_x * scale,
                 i.pnt[1].y + i.pnt[1].tot_displ_y * scale),
                (i.pnt[2].x + i.pnt[2].tot_displ_x * scale,
                 i.pnt[2].y + i.pnt[2].tot_displ_y * scale),
                (0, 0)]
            codes_nonplast += [Path.MOVETO] + [Path.LINETO]*2 + [Path.CLOSEPOLY]

        if len(i.pnt) == 4:
            polygons_quad += [
                (i.pnt[0].x + i.pnt[0].tot_displ_x * scale,
                 i.pnt[0].y + i.pnt[0].tot_displ_y * scale),
                (i.pnt[1].x + i.pnt[1].tot_displ_x * scale,
                 i.pnt[1].y + i.pnt[1].tot_displ_y * scale),
                (i.pnt[2].x + i.pnt[2].tot_displ_x * scale,
                 i.pnt[2].y + i.pnt[2].tot_displ_y * scale),
                (i.pnt[3].x + i.pnt[3].tot_displ_x * scale,
                 i.pnt[3].y + i.pnt[3].tot_displ_y * scale),
                (0, 0)]
            codes_quad += [Path.MOVETO] + [Path.LINETO] * 3 + [Path.CLOSEPOLY]

    polygons_plast = np.array(polygons_plast, float)
    try:
        path_plast = Path(polygons_plast, codes_plast)
        pathpatch_PLAST = PathPatch(path_plast, facecolor='#fff000')
    except ValueError:
        pathpatch_PLAST = False
        print('не получилось')
    else:
        path_plast = Path(polygons_plast, codes_plast)
        pathpatch_PLAST = PathPatch(path_plast, facecolor='#fff000')

    polygons_nonplast = np.array(polygons_nonplast, float)
    try:
        path_nonplast = Path(polygons_nonplast, codes_nonplast)
        pathpatch_nonPLAST = PathPatch(path_nonplast, facecolor='#ffffff')
    except ValueError:
        pathpatch_nonPLAST = False
        print('не получилось')
    else:
        path_nonplast = Path(polygons_nonplast, codes_nonplast)
        pathpatch_nonPLAST = PathPatch(path_nonplast, facecolor='#ffffff')

    polygons_quad = np.array(polygons_quad, float)
    try:
        path_quad = Path(polygons_quad, codes_quad)
        pathpatch_quad = PathPatch(path_quad, facecolor='#ffffff')
    except ValueError:
        pathpatch_quad = False
        print('не получилось')
    else:
        path_quad = Path(polygons_quad, codes_quad)
        athpatch_quad = PathPatch(path_quad, facecolor='#ffffff')

    return pathpatch_PLAST, pathpatch_nonPLAST, pathpatch_quad

def make_add_txt(elms, pnts, scale=0, num_elms=False, num_pnts=False,
                displacement_y=False, displacement_x=False, sigma_eq=False,
                sigma_x=False):
    if num_elms == True:
        for i in elms:  # рисует номер элемента
            plt.text((i.pnt[0].x + i.pnt[1].x + i.pnt[2].x)/3+0.001,
                     (i.pnt[0].y + i.pnt[1].y + i.pnt[2].y)/3,
                     i.num, color='red', fontsize=6)

    if num_pnts == True:
        for i in pnts: # рисует номер узла
            plt.text(i.x + 0.1,
                     i.y + 0.1,
                     i.num, color='blue', fontsize=6)
    
    if displacement_y == True:
        for i in pnts:
            plt.text(i.x + 0.1,
                     i.y + 0.1,
                     '{:.2f}'.format(i.tot_displ_y*1000), color='blue', fontsize=6)
    if displacement_x == True:
        for i in pnts:
            plt.text(i.x + 0.001,
                     i.y + 0.001,
                     '{:.2f}'.format(i.tot_displ_x*1000), color='blue', fontsize=6)
    if sigma_eq == True:
        for i in elms:
            plt.text((i.pnt[0].x + i.pnt[1].x + i.pnt[2].x)/3,
                     (i.pnt[0].y + i.pnt[1].y + i.pnt[2].y)/3,
                     '{:.0f}'.format(i.sigma_eq),
                     color='blue',
                     fontsize=6)
    if sigma_x == True:
        for i in elms:
            plt.text((i.pnt[0].x + i.pnt[1].x + i.pnt[2].x)/3,
                     (i.pnt[0].y + i.pnt[1].y + i.pnt[2].y)/3,
                     '{:.0f}'.format(i.sigma_x),
                     color='blue',
                     fontsize=6)

def show_plast(le, lp):
    fig, ax = plt.subplots()
    a = make_plot_plast(le, 1)
    make_add_txt(le, lp, scale=1, num_elms=True)
    for i in range(3):
        if a[i] is not False:
            ax.add_patch(a[i])
    #~ ax.set_title(
        #~ 'Возникновение пластических деформации при {:.2f}P'.format(count))
    ax.autoscale_view()
    # fig.savefig('png/png{}.png'.format(int(count*100)))
    plt.show()



def graph_sigma(l_s1, l_s2, sigma_pr):
    fig, ax = plt.subplots()
    x = np.arange(-sigma_pr/50,sigma_pr/50,1)
    ax.plot(x, (x-np.sqrt(-3*x**2+4*(sigma_pr/100)**2))/2, label='se')
    ax.plot(x, (x+np.sqrt(-3*x**2+4*(sigma_pr/100)**2))/2, label='se')
    ax.plot(l_s1, l_s2, 'ro', label='se')
    for i in range(len(l_s1)):
            plt.text(l_s1[i] + 0.1,
                     l_s2[i] + 0.1,
                     '{}'.format(i), color='blue', fontsize=6)
    plt.show()


def show_graph_P_u(P,u):
    fig, ax = plt.subplots()
    ax.plot(u, P, label='u(p)')

    plt.show()
    
    
def graph_show(l_se, l_scr):
    '''Show graphic of Se and Scr'''
    fig, ax = plt.subplots()
    x = np.arange(0, 100, 1)
    ax.plot(x, l_se, 'ro', label='se')
    ax.plot(x, l_scr, 'k', label='scr')
    plt.show()

if __name__=="__main__":
    graph_sigma([0],[0])
