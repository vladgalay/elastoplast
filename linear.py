#! /usr/bin/env python
# -*- coding: utf-8 -*-

import liraparser
import classes
import graph_plot as gp
import nonlinear
import numpy as np
import math
import graph
import matplotlib.pyplot as plt

#newElement - булевая переменная. Использовать ли новый КЭ пластины при расчёте
def matrix_K(le, lp, newElement):
    """make matrix of stiffness with boundary condition"""
    size = len(lp)
    K = np.zeros((size * 3, size * 3))
    # определение геометрических параметров

    for elm in le:
        """Для стержней"""
        if len(elm.pnt) == 2:
            #print('{type = 2, points = {', elm.pnt[0].num, ', ', elm.pnt[1].num, '}},')
            #  k11 | k12
            #  k21 | k22
			
            """Закреплены ли 0 и 1 узлы
            имеется ввиду, принадлежат ли эти узлы пластинам"""
            #fix = [0, 0]
            
            """Сразу же проверяем принадлежность узлов пластинам"""
            """for i in range(2):
                if K[3*(elm.pnt[i].num-1) + 2, 3*(elm.pnt[i].num-1) + 2] == 1:
                    fix[i] = 1"""
            
            """Матрица жёсткости элемента"""
            K_e = elm.matrix_K()
            
            #print("---СТЕРЖЕНЬ---")
            #print(K_e)
            
            # k11
            K[3*(elm.pnt[0].num-1):3*elm.pnt[0].num,
              3*(elm.pnt[0].num-1):3*elm.pnt[0].num] += \
                K_e[0:3,0:3]
            # k12
            K[3*(elm.pnt[0].num-1):3*elm.pnt[0].num,
              3*(elm.pnt[1].num-1):3*elm.pnt[1].num] += \
                K_e[0:3,3:6]
            #  k21
            K[3*(elm.pnt[1].num-1):3*elm.pnt[1].num,
              3*(elm.pnt[0].num-1):3*elm.pnt[0].num] += \
                K_e[3:6,0:3]
            #  k22
            K[3*(elm.pnt[1].num-1):3*elm.pnt[1].num,
              3*(elm.pnt[1].num-1):3*elm.pnt[1].num] += \
                K_e[3:6,3:6]
            
            """Если имеет место принадлежность узла к пластине"""
            """for i in range(2):
                if fix[i] == 1:
                    K[3*(elm.pnt[i].num-1) + 2] = 0
                    K[:, 3*(elm.pnt[i].num-1) + 2] = 0
                    K[3*(elm.pnt[i].num-1) + 2,
                    3*(elm.pnt[i].num-1) + 2] = 1"""
            
        if len(elm.pnt) == 3:
            #print('{type = 3, points = {', elm.pnt[0].num, ', ', elm.pnt[1].num, ', ', elm.pnt[2].num, '}},')
            
            #Если по старому КЭ
            if newElement == False:
                #  k11 | k12 | k13
                #  k21 | k22 | k23
                #  k31 | k32 | k33
                if elm.pl==False:
                    D = elm.D
                if elm.pl==True:
                    D = elm.D_ep
                #~ D = elm.D
                # матрица неизвестных
                B = elm.matrix_B()[0]
                a = elm.pnt[1].x+elm.pnt[1].tot_displ_x \
                    - elm.pnt[0].x-elm.pnt[0].tot_displ_x
                b = elm.pnt[0].y+elm.pnt[0].tot_displ_y \
                    - elm.pnt[1].y-elm.pnt[1].tot_displ_y
                c = elm.pnt[0].x+elm.pnt[0].tot_displ_x \
                    - elm.pnt[2].x-elm.pnt[2].tot_displ_x
                d = elm.pnt[2].y+elm.pnt[2].tot_displ_y \
                    - elm.pnt[0].y-elm.pnt[0].tot_displ_y
                s = -b - d
                t = -a - c
                area = abs(a*d - b*c)/2
                K_e = np.transpose(B) @ D @ B * elm.h * area
                #print('Толщина пластины #', elm.num, ' = ', elm.h, ' м')
                
                #print(K_e)
                #print("---ПЛАСТИНА---")
                #print(K_e)
                
                # k11
                K[3*(elm.pnt[0].num-1):3*elm.pnt[0].num - 1,
                  3*(elm.pnt[0].num-1):3*elm.pnt[0].num - 1] += \
                    K_e[0:2,0:2]
                # k12
                K[3*(elm.pnt[0].num-1):3*elm.pnt[0].num - 1,
                  3*(elm.pnt[1].num-1):3*elm.pnt[1].num - 1] += \
                    K_e[0:2,2:4]
                #  k13
                K[3*(elm.pnt[0].num-1):3*elm.pnt[0].num - 1,
                  3*(elm.pnt[2].num-1):3*elm.pnt[2].num - 1] += \
                    K_e[0:2,4:6]
                #  k21
                K[3*(elm.pnt[1].num-1):3*elm.pnt[1].num - 1,
                  3*(elm.pnt[0].num-1):3*elm.pnt[0].num - 1] += \
                    K_e[2:4,0:2]
                #  k22
                K[3*(elm.pnt[1].num-1):3*elm.pnt[1].num - 1,
                  3*(elm.pnt[1].num-1):3*elm.pnt[1].num - 1] += \
                    K_e[2:4,2:4]
                # k23
                K[3*(elm.pnt[1].num-1):3*elm.pnt[1].num - 1,
                  3*(elm.pnt[2].num-1):3*elm.pnt[2].num - 1] += \
                    K_e[2:4,4:6]
                # k31
                K[3*(elm.pnt[2].num-1):3*elm.pnt[2].num - 1,
                  3*(elm.pnt[0].num-1):3*elm.pnt[0].num - 1] += \
                    K_e[4:6,0:2]
                # k32
                K[3*(elm.pnt[2].num-1):3*elm.pnt[2].num - 1,
                  3*(elm.pnt[1].num-1):3*elm.pnt[1].num - 1] += \
                    K_e[4:6,2:4]
                # k33
                K[3*(elm.pnt[2].num-1):3*elm.pnt[2].num - 1,
                  3*(elm.pnt[2].num-1):3*elm.pnt[2].num - 1] += \
                    K_e[4:6,4:6]
    			
                """Присваивем единички - это, как оказалось, помогает"""
                """for p0 in range(3):
                    K[3*(elm.pnt[p0].num-1) + 2] = 0
                    K[:, 3*(elm.pnt[p0].num-1) + 2] = 0
                    K[3*(elm.pnt[p0].num-1) + 2,
                    3*(elm.pnt[p0].num-1) + 2] += 1"""
            #И если по новому КЭ
            if newElement == True:
                #Получаем матрицу жёсткости
                K_e = elm.matrix_K()
                #И вставляем её в глобальную матрицу жёсткости
                for iy in range(3):
                    for ix in range(3):
                        K[3*(elm.pnt[iy].num-1):3*(elm.pnt[iy].num-1) + 3, 3*(elm.pnt[ix].num-1):3*(elm.pnt[ix].num-1) + 3] += K_e[iy*3:iy*3 + 3, ix*3:ix*3 + 3]
        
        if len(elm.pnt) == 4:
            # прямоугольные КЭ
            #  k11 | k12 | k13 | k14
            #  k21 | k22 | k23 | k24
            #  k31 | k32 | k33 | k34
            #  k41 | k42 | k43 | k44
            a = elm.pnt[1].x - elm.pnt[0].x
            b = elm.pnt[2].y - elm.pnt[1].y
            k = elm.E*elm.h/(1-elm.V**2)
            mu = (1-elm.V)/2
            c = a/(6*b)
            cm = mu*c
            d = b/(6*a)
            dm = mu*d
            s = elm.V/4
            t = mu/4
            # k11
            K[3 * (elm.pnt[0].num - 1):3 * elm.pnt[0].num - 1,
              3 * (elm.pnt[0].num - 1):3 * elm.pnt[0].num - 1] += \
                k * np.array([[2*(d+cm), s+t],
                              [s+t, 2*(c+dm)]])
            # k12
            K[3 * (elm.pnt[0].num - 1):3 * elm.pnt[0].num - 1,
              3 * (elm.pnt[1].num - 1):3 * elm.pnt[1].num - 1] += \
                k * np.array([[cm-2*d, s-t],
                              [t-s, c-2*dm]])
            # k13
            K[3 * (elm.pnt[0].num - 1):3 * elm.pnt[0].num - 1,
              3 * (elm.pnt[2].num - 1):3 * elm.pnt[2].num - 1] += \
                k * np.array([[-d-cm, -s-t],
                              [-s-t, -c-dm]])
            # k14
            K[3 * (elm.pnt[0].num - 1):3 * elm.pnt[0].num - 1,
              3 * (elm.pnt[3].num - 1):3 * elm.pnt[3].num - 1] += \
                k * np.array([[d-2*cm, t-s],
                              [s-t, dm-2*c]])
            # k21
            K[3 * (elm.pnt[1].num - 1):3 * elm.pnt[1].num - 1,
              3 * (elm.pnt[0].num - 1):3 * elm.pnt[0].num - 1] += \
                k * np.array([[cm-2*d, t-s],
                              [s-t, c-2*dm]])
            # k22
            K[3 * (elm.pnt[1].num - 1):3 * elm.pnt[1].num - 1,
              3 * (elm.pnt[1].num - 1):3 * elm.pnt[1].num - 1] += \
                k * np.array([[2*(d+cm), -s-t],
                              [-s-t, 2*(c+dm)]])
            # k23
            K[3 * (elm.pnt[1].num - 1):3 * elm.pnt[1].num - 1,
              3 * (elm.pnt[2].num - 1):3 * elm.pnt[2].num - 1] += \
                k * np.array([[d-2*cm, s-t],
                              [t-s, dm-2*c]])
            # k24
            K[3 * (elm.pnt[1].num - 1):3 * elm.pnt[1].num - 1,
              3 * (elm.pnt[3].num - 1):3 * elm.pnt[3].num - 1] += \
                k * np.array([[-d-cm, s+t],
                              [s+t, -c-dm]])
            # k31
            K[3 * (elm.pnt[2].num - 1):3 * elm.pnt[2].num - 1,
              3 * (elm.pnt[0].num - 1):3 * elm.pnt[0].num - 1] += \
                k * np.array([[-d-cm, -s-t],
                              [-s-t, -c-dm]])
            # k32
            K[3 * (elm.pnt[2].num - 1):3 * elm.pnt[2].num - 1,
              3 * (elm.pnt[1].num - 1):3 * elm.pnt[1].num - 1] += \
                k * np.array([[d-2*cm, t-s],
                              [s-t, dm-2*c]])
            # k33
            K[3 * (elm.pnt[2].num - 1):3 * elm.pnt[2].num - 1,
              3 * (elm.pnt[2].num - 1):3 * elm.pnt[2].num - 1] += \
                k * np.array([[2*(d+cm), s+t],
                              [s+t, 2*(c+dm)]])
            # k34
            K[3 * (elm.pnt[2].num - 1):3 * elm.pnt[2].num - 1,
              3 * (elm.pnt[3].num - 1):3 * elm.pnt[3].num - 1] += \
                k * np.array([[cm-2*d, s-t],
                              [t-s, c-2*dm]])
            # k41
            K[3 * (elm.pnt[3].num - 1):3 * elm.pnt[3].num - 1,
              3 * (elm.pnt[0].num - 1):3 * elm.pnt[0].num - 1] += \
                k * np.array([[d-2*cm, s-t],
                              [t-s, dm-2*c]])
            # k42
            K[3 * (elm.pnt[3].num - 1):3 * elm.pnt[3].num - 1,
              3 * (elm.pnt[1].num - 1):3 * elm.pnt[1].num - 1] += \
                k * np.array([[-d-cm, s+t],
                              [s+t, -c-dm]])
            # k43
            K[3 * (elm.pnt[3].num - 1):3 * elm.pnt[3].num - 1,
              3 * (elm.pnt[2].num - 1):3 * elm.pnt[2].num - 1] += \
                k * np.array([[cm-2*d, t-s],
                              [s-t, c-2*dm]])
            # k44
            K[3 * (elm.pnt[3].num - 1):3 * elm.pnt[3].num - 1,
              3 * (elm.pnt[3].num - 1):3 * elm.pnt[3].num - 1] += \
                k * np.array([[2*(d+cm), -s-t],
                              [-s-t, 2*(c+dm)]])
            
            """Копипаст с предыдущего условия, только вместо 3 - 4
            Присваивем единички - это, как оказалось, помогает"""
            """for p0 in range(4):
                K[3*(elm.pnt[p0].num-1) + 2] = 0
                K[:, 3*(elm.pnt[p0].num-1) + 2] = 0
                K[3*(elm.pnt[p0].num-1) + 2,
                3*(elm.pnt[p0].num-1) + 2] = 1"""
        """Матрица до граничных условий"""
        #print("---ДО ГРАНИЧНЫХ УСЛОВИЙ---")
        #print(K)
    """Проходимся по главной диагонали и ставим единички где нужно
    ПРИШЛО ВРЕМЯ АДСКИХ ЕОСТЫЛЕЙ
    КОСТЫЛИ САМИ СЕБЯ НЕ НАПИШУТ"""
    #Это если старый КЭ
    if newElement == False:
        for i in range(int(K.shape[0]/3)):
            """Локальная переменная - принадлежит ли данному узлу стержень"""
            t = False
            """Проверка на принадлежность стержня данному узлу"""
            for elm in le:
                """Если принадлежит - сразу обрываем внутренний цикл и присваиваем значение переменной"""
                if len(elm.pnt) == 2 and (elm.pnt[0].num - 1 == i or elm.pnt[1].num - 1 == i):
                    #print("НЕ ПОПАЛСО СТЕРЖЕНЬ ", elm.num, "ЭЛЕМЕНТ ", i*3 + 2)
                    t = True
                    break
                #print("ПОПАЛСО ПЛАСТИНА ", elm.num, "ЭЛЕМЕНТ ", i*3 + 2)
            
            """Далее проверяем переменную. Если принадлежит - заполняем единицу"""
            if not t:
                #print("ЗАПОЛНИЛИ ", i*3 + 2)
                K[i*3 + 2, i*3 + 2] = 1
    #print(K)
    # присвоение граничных условий
    for pnt in lp:
        if pnt.bc != None:
            #print('BC: ', pnt.bc)
            for i in range(3):
                if i*2 + 1 in pnt.bc:
                    K[3*(pnt.num-1)+i] = 0
                    K[:, 3*(pnt.num-1)+i] = 0
                    K[3*(pnt.num-1)+i, 3*(pnt.num-1)+i] = 1
        """if pnt.bc == (1,):
            K[3*(pnt.num-1)] = 0
            K[:, 3*(pnt.num-1)] = 0
            K[3*(pnt.num-1), 3*(pnt.num-1)] = 1
        if pnt.bc == (3,):
            K[3*(pnt.num-1)+1] = 0
            K[:, 3*(pnt.num-1)+1] = 0
            K[3*(pnt.num-1)+1, 3*(pnt.num-1)+1] = 1
        if pnt.bc == (1, 3):
            K[3*(pnt.num-1)] = 0
            K[:, 3*(pnt.num-1)] = 0
            K[3*(pnt.num-1), 3*(pnt.num-1)] = 1
            K[3*(pnt.num-1)+1] = 0
            K[:, 3*(pnt.num-1)+1] = 0
            K[3*(pnt.num-1)+1, 3*(pnt.num-1)+1] = 1
        if pnt.bc == (1, 3, 5):
            K[3*(pnt.num-1)] = 0
            K[:, 3*(pnt.num-1)] = 0
            K[3*(pnt.num-1), 3*(pnt.num-1)] = 1
            K[3*(pnt.num-1)+1] = 0
            K[:, 3*(pnt.num-1)+1] = 0
            K[3*(pnt.num-1)+1, 3*(pnt.num-1)+1] = 1
            K[3*(pnt.num-1)+2] = 0
            K[:, 3*(pnt.num-1)+2] = 0
            K[3*(pnt.num-1)+2, 3*(pnt.num-1)+2] = 1"""
    return K


def matrix_P(size, loads, lp):
    P = np.zeros(size*3)
    for load in loads:
        #print('LOAD ', load)
        if load[1] == 1:
            P[3*(load[0]-1)] -= load[2]
            #print('НАГРУЗКА ПО X\nP = ', P)
        if load[1] == 3:
            P[3*(load[0]-1) + 1] -= load[2]
        if load[1] == 5:
            P[3*(load[0]-1) + 2] -= load[2]
    for point in lp:
        if point.bc == (1,): P[3*(point.num-1)] = 0
        if point.bc == (3,): P[3*point.num-1] = 0
        if point.bc == (1, 3):
            P[3*(point.num-1)] = 0
            P[3*point.num-1] = 0
        if point.bc == (1, 3, 5):
            P[3*(point.num-1)] = 0
            P[3*point.num-2] = 0
            P[3*point.num-1] = 0
    return P

#newElement - булевая переменная. Использовать ли новый КЭ пластины при расчёте
#l, a - значения констант
def matrix_U(le, lp, loads, newElement, l, a):
    classes.l = l
    classes.a = a
    
    K = matrix_K(le, lp, newElement)
    #print(K)
    P = matrix_P(len(lp), loads, lp)
    #print(P)
    U = np.linalg.solve(K, P.transpose())
    """for count, i in enumerate(U):
        print('узел: {}   перемещения {}'.format(count//3 + 1, i))"""
    return U

def outp(n_elm):
    print('\u03c3x = {}'.format(le[n_elm].sx))
    print('\u03c3y = {}'.format(le[n_elm].sy))
    print('\u03C4 = {}'.format(le[n_elm].tau))


if __name__ == '__main__':
    le, lp = liraparser.list_of_elem()
    loads = liraparser.list_of_loads()
    with classes.Timer() as p:
        U = matrix_U(le, lp, loads, False, 0.1, 0.5)
        #Далее заносим полученные перемещения в точки
        for pnt in lp:
            pnt.displ_x = U[(pnt.num - 1)*3]
            pnt.displ_y = U[(pnt.num - 1)*3 + 1]
        
        #И считаем напряжения
        for elm in le:
            #Если стержень - то вызываем функцию с аргументом U
            if len(elm.pnt) == 2:
                elm.define_stress(U)
            #Если пластина - то без U
            else:
                elm.define_stress()
        #print('USHAPE', U.shape)
        #print(np.reshape(U, (int(U.shape[0]/3), 3)))
        """Выводим перемещения узлов в формате, который будет воспринимать Lua"""
        """for i in range(int(U.shape[0]/3)):
            print('{', U[i*3], ', ', U[i*3 + 1], '},')"""
        """for i in range(121):
            print('{', lp[i].x, ', ', lp[i].y, '},')"""
        #Графический вывод
        #graph.PrepareGraph()
        #Подготавливаем вектор координат узлов в нужном формате
        Pc = np.zeros(U.shape[0])
        #Заполняем его
        for point in lp:
            Pc[(point.num - 1)*3 + 0] = point.x
            Pc[(point.num - 1)*3 + 1] = point.y
        #И рисуем
        #Коэффициент увеличения
        kk = 300
        #Оригинал
        graph.AddGraph(le, Pc, np.zeros_like(U), kk, 'b-')
        #С перемещениями пластин с двумя степенями свободы
        graph.AddGraph(le, Pc, U, kk, 'r-')
        #И с новым типом КЭ
        """nK = matrix_K(le, lp, True)
        P = matrix_P(len(lp), loads, lp)
        nU = np.linalg.solve(nK, P.transpose())"""
        #nU = matrix_U(le, lp, loads, True, 1.5, -0.93)
        #nU = matrix_U(le, lp, loads, True, 1.0, -1.048)
        nU = matrix_U(le, lp, loads, True, 0.1, 1.0)
        #graph.AddContour(le, Pc, nU, kk, nU.reshape((int(nU.shape[0]/3), 3)).transpose()[2])
        """graph.AddGraph(le, Pc, nU, kk, 'g-')"""
        #Делаем массив векторов координат и моментов (0 - X, 1 - Z, 2 - uY)
        
        print('U:\n', U.reshape((int(U.shape[0]/3), 3)))
        print('nU:\n', nU.reshape((int(nU.shape[0]/3), 3)))
        """Av = np.zeros((3, int(U.shape[0]/3)))
        #Заполняем по точкам
        for point in lp:
            Av[0, (point.num - 1)] = point.x
            Av[1, (point.num - 1)] = point.y
        #И моментам
        for j in range(int(U.shape[0]/3)):
            Av[2, j] = nU[j*3 + 2]
        #Попытка в изолинии моментов
        plt.contour(Av[0], Av[1], Av, 1, cmap = 'RdGy')"""
        graph.ShowGraph()
        for i in lp:
            i.tot_displ_x += U[3 * (i.num-1)]
            i.tot_displ_y += U[3 * i.num - 1]
        for elm in le:
            #~ for iteration in range(iters):
            if len(elm.pnt) == 2:
                elm.add_stress(*elm.define_stress(U))
            else:
                elm.add_stress(*elm.define_stress())
    #print(classes.make_q(le, lp))
    #gp.show_plast(le, lp)
