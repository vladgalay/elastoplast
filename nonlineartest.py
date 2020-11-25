import liraparser
import linear as linear
import graph_plot as gp
import numpy as np
import classes

def select_P():

    pass


def iterate_alg():
    pass


def iteration(P_init, le, lp, steps, iters):
    K = None
    U = np.zeros(len(lp)*2)
    P = np.zeros(len(lp)*2)
    delta=1/steps
    count = delta
    Pi = P_init*0.95*delta
    P += Pi
    s1 = []  # Сигма 1
    s2 = []  # Сигма 2
    sig_x = list() # сигма Х
    def_x = list() # деформации Х
    u = list()
    for j in range(steps):
        
        K = linear.matrix_K(le, lp)
        # пересчет матрицы К с новой геометрией
        U = np.linalg.solve(K, Pi.transpose())
        for i in lp:
            i.displ_x = U[2 * (i.num-1)]
            i.displ_y = U[2 * i.num - 1]
        for k in range(iters):
            # номер итерации кратный 10
            #~ if k % 10 == 0: print('Итерация №', k)
            for elm in le:
                elm.define_stress()
                # ~ if elm.num == 1: print('sx=', elm.d_sigma_x)
            q = classes.make_q(le, lp)
            r = P - q.reshape(len(lp)*2)
            #~ if k % 10 == 0: print('q=', q)
            #~ if k % 10 == 0: print('p=', P)
            # накладка граничных условий на матрицу r(невязка внешних и внутренних сил)
            for point in lp:
                if point.bc == (1,): r[2*(point.num-1)] = 0
                if point.bc == (3,): r[2*point.num-1] = 0
                if point.bc == (1,3):
                    r[2*(point.num-1)] = 0
                    r[2*point.num-1] = 0
            #~ print(np.linalg.norm(r), ' <= ', 0.001*np.linalg.norm(P))
            if np.linalg.norm(r) <= 0.001*np.linalg.norm(P):
                #~ print('break from {} iteration'.format(k))
                break
            else:
                d_U = np.linalg.solve(K, r.transpose())
                for point in lp:
                    
                    point.displ_x += float(d_U[2 * (point.num-1)])
                    point.displ_y += float(d_U[2 * point.num - 1])

        for elm in le:
            #~ for iteration in range(iters):
            elm.add_stress(*elm.define_stress())


        for i in lp:
            i.summ_displ()
        # show_plast(count)
        ne = 9
        s1.append(le[ne].sigma_1/100)
        s2.append(le[ne].sigma_2/100)
        sig_x.append(le[ne].sigma_y)
        le[ne].define_strains()
        def_x.append(le[ne].strain[1,0])
        if j % 1 == 0:
            print('step:{}\n\nlist of sigma:\n{:10.0f} тс/м**2\nlist of strain:\n{:10.5f}'
                .format(j,le[ne].sigma_y,le[ne].strain[1,0]))
        u_num = 9
        u.append(abs(np.sqrt(lp[u_num].tot_displ_y**2+lp[u_num].tot_displ_x**2)))
        #~ if j in range(87,91):
            #~ print(le[0].D,'-', le[0].D_ep)
        print('{:.0%}'.format(count/1))
        P += Pi
        count += delta
    gp.show_graph_P_u(np.arange(steps),u)
    #~ print('list of sigma:\n{}\nlist of strain:\n{}'
        #~ .format(sig_x,def_x))
    return s1, s2


if __name__ == '__main__':
    le, lp = liraparser.list_of_elem()
    loads = liraparser.list_of_loads()
    P_init = linear.matrix_P(len(lp), loads, lp)
    l_s1, l_s2 = iteration(P_init, le, lp, 200, 1)
    gp.graph_sigma(l_s2, l_s1, le[0].sigma_pr)
    #~ qp.multigraph(l_se, l_scr)
    gp.show_plast(le, lp)
