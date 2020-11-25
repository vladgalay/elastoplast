import liraparser
import linear as linear
import graph_plot as gp
import numpy as np


def select_P():

    pass
def iteration(P_init, le, lp, steps):
    K = None
    U = np.zeros(len(lp)*2)
    delta=1/steps
    count = delta
    P = P_init*delta
    s1 = []  # Сигма критическое
    s2 = []  # Сигма эквивалентное
    u = list()
    for j in range(steps):
        
        K = linear.matrix_K(le, lp)
        # пересчет матрицы К с новой геометрией
        U = np.linalg.solve(K, P.transpose())
        for i in lp:
            i.displ_x = U[2 * (i.num-1)]
            i.displ_y = U[2 * i.num - 1]

        for elm in le:
            elm.add_stress(*elm.define_stress())

        for i in lp:
            i.summ_displ()
        # show_plast(count)
        ne = 0
        s1.append(le[ne].sigma_1/100)
        s2.append(le[ne].sigma_2/100)
        
        u.append(abs(np.sqrt(lp[1].tot_displ_y**2+lp[1].tot_displ_x**2)))
        #~ if j in range(87,91):
            #~ print(le[0].D,'-', le[0].D_ep)
        print('{:.0%}'.format(count/1))
        
        count += delta

    gp.show_graph_P_u(np.arange(steps),u)
    return s1, s2


if __name__ == '__main__':
    print('Начинаем')
    le, lp = liraparser.list_of_elem()
    print('Загрузили список элементов')
    loads = liraparser.list_of_loads()
    print('Загрузили список загружений')
    P_init = linear.matrix_P(len(lp), loads)
    print('Получаем матрицу P')
    l_scr, l_se = iteration(P_init, le, lp, 20)
    print('Проитерировали')
    gp.graph_sigma(l_se, l_scr)
    
    gp.show_plast(le, lp)
