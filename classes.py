#! /usr/bin/env python
# -*- coding: utf-8 -*-

import liraparser
import numpy as np
import time
import math

l = 0.1
a = 0.5

#Функция изменения значений констант
"""def SetConst(nl, na):
    l = nl
    a = na"""

class Fe:
    """ Класс описывающий поведение элементов """
    def __init__(self, num, types, n_stiff, *pnt):
        """присвоение характеристик при работе скрипта liraparser.py"""
        self.num = num
        self.types = types
        self.n_stiff = n_stiff
        self.pnt = pnt
        self.pl = False
        self.sigma_x = 0
        self.sigma_y = 0
        self.tau_xy = 0
        self.d_sigma_x = 0
        self.d_sigma_y = 0
        self.d_tau_xy = 0
        self.sigma_eq = 0
        self.sigma_pr = 44300 # T/m**2
        #Коэффициент Пуассона
        self.nu = 0.25
        #Константы
        #self.l = 0.1
        #self.a = 0.5


    def add_proporties(self, E, V, h):
        """присвоение характеристик при работе скрипта liraparser.py"""
        self.V = V
        self.E = E
        #print('E пластины #', self.num, ' = ', self.E , 'т/м2')
        #Присваиваем модуль сдвига
        self.G = self.E/(2*(1 + self.nu))
        self.Et = 340000
        self.h = h
        if liraparser.TASK == 'plane stress':
                    self.D=self.E/(1-self.V**2)*np.array([[1,self.V,0],
                                                    [self.V,1,0],
                                                    [0,0,(1-self.V)/2]])
                    #print('D пластины #', self.num, ' =\n', self.D)
        if liraparser.TASK == 'plane strain':
                    self.D=self.E/(1+self.V)/(1-2*self.V)*\
                    np.array([[1-self.V,self.V,0],
                             [self.V,1-self.V,0],
                             [0,0,(1-2*self.V)/2]])
        # self.coh = 0.3  # T/m**2
        # self.fi = math.pi*20/180  # friction angle


    def matrix_B(self):
        """Построение матрицы функций формы и вычисление площади КЭ"""
        #for i in range(3):
            #print('Точка ', i, ':\n\tНомер: ', self.pnt[i].num, '\n\tX: ', self.pnt[i].x, '\n\tY: ', self.pnt[i].y)
        a = self.pnt[1].x+self.pnt[1].tot_displ_x \
            - self.pnt[0].x-self.pnt[0].tot_displ_x
        b = self.pnt[0].y+self.pnt[0].tot_displ_y \
            - self.pnt[1].y-self.pnt[1].tot_displ_y
        c = self.pnt[0].x+self.pnt[0].tot_displ_x \
            - self.pnt[2].x-self.pnt[2].tot_displ_x
        d = self.pnt[2].y+self.pnt[2].tot_displ_y \
            - self.pnt[0].y-self.pnt[0].tot_displ_y
        #print('a, b, c, d = ', a, ', ', b, ', ', c, ', ', d,' м')
        s = -b - d
        t = -a - c
        area = abs(a*d - b*c)/2
        #print('s, t = ', s, ', ', t, ' м\nПлощадь пастины = ', area, 'м2')
        B = np.array([[s, 0, d, 0, b, 0],
                      [0, t, 0, c, 0, a],
                      [t, s, c, d, a, b]])/2/area
        #print('B пластины #', self.num, ' =\n', B.astype(int))
        return B, area
    
    #Построение матрицы жёсткости элемента
    def matrix_K(self):
        #1. Определение длин сторон пластины
        #Номер стороны = номер противоположного узла пластины
        l1 = math.sqrt((self.pnt[1].x - self.pnt[2].x)**2 + (self.pnt[1].y - self.pnt[2].y)**2)
        l2 = math.sqrt((self.pnt[0].x - self.pnt[2].x)**2 + (self.pnt[0].y - self.pnt[2].y)**2)
        l3 = math.sqrt((self.pnt[1].x - self.pnt[0].x)**2 + (self.pnt[1].y - self.pnt[0].y)**2)
        #2. Определение дополнительных параметров, связанных с размером пластины
        l12 = l1**2 + l2**2 - l3**2
        l23 = l2**2 + l3**2 - l1**2
        l31 = l3**2 + l1**2 - l2**2
        #3. Определение метрического тензора (что это?)
        gab = np.array([[l2**2, l12/2],
                        [l12/2, l1**2]])
        #4. Вычисление определителя метрического тензора
        g = np.linalg.det(gab)
        #5. Вычисление площади пластины через определитель метрического тензора
        A = math.sqrt(g/4)
        #6. Задаём известные матрицы
        B = np.array([[1, 0, 0, 0, -1, 0],
                      [0, 0, 0, 1, 0, -1],
                      [0, 1, 0, 0, 0, -1],
                      [0, 0, 1, 0, -1, 0]])
        D = self.G*np.array([[(2*(1 - self.nu))/(1 - 2*self.nu), (2*self.nu)/(1 - 2*self.nu), 0, 0],
                        [(2*self.nu)/(1 - 2*self.nu), (2*(1 - self.nu))/(1 - 2*self.nu), 0, 0],
                        [0, 0, 1 + a, 1 - a],
                        [0, 0, 1 - a, 1 + a]])
        Bcap = np.array([[1, 0, -1],
                         [0, 1, -1]])
        Dcap = (4*self.G*l**2)*np.array([[1, 0],
                                    [0, 1]])
        #7. Задаём константу
        C = 2*a
        #8. Вычисляем подматрицы
        Kuu = A*(B.transpose() @ D @ B)
        Kufi = ((2*C*A**2)/3)*np.array([[0, 0, 0],
                                        [-1, -1, -1],
                                        [1, 1, 1],
                                        [0, 0, 0],
                                        [-1, -1, -1],
                                        [1, 1, 1]])
        Kfifi = ((4*C*A**3)/3)*np.array([[1, 0.5, 0.5],
                                         [0.5, 1, 0.5],
                                         [0.5, 0.5, 1]]) + A*(Bcap.transpose() @ Dcap @ Bcap)
        #9. Составляем матрицу жесткости элемента
        Ke = np.zeros((9, 9))
        Ke[0:6, 0:6] = Kuu
        Ke[0:6, 6:9] = Kufi
        Ke[6:9, 0:6] = Kufi.transpose()
        Ke[6:9, 6:9] = Kfifi
        #10. Переставляем столбцы и строки матрицы из формата
        #x1, z1, x2, z2, x3, z3, uy1, uy2, uy3
        #в
        #x1, z1, uy1, x2, z2, uy2, x3, z3, uy3
        #Для этого составляем новую матрицу жёсткости и заполняем её
        #Kr - K результирующая
        Kr = np.zeros((9, 9))
        #Заполняем матрицу по подматрицам 3х3
        for iy in range(3):
            for ix in range(3):
                #Кусок матрицы Kuu
                Kr[iy*3:iy*3 + 2, ix*3:ix*3 + 2] = Ke[iy*2:iy*2 + 2, ix*2:ix*2 + 2]
                #Далее матрица Kufi
                Kr[iy*3:iy*3 + 2, ix*3 + 2:ix*3 + 3] = Ke[iy*2:iy*2 + 2, 6 + ix:6 + ix + 1]
                #Транспонированная матрица Kufi
                Kr[iy*3 + 2:iy*3 + 3, ix*3:ix*3 + 2] = Ke[6 + iy:6 + iy + 1, ix*2:ix*2 + 2]
                #И матрица Kfifi - тут всего один компонент
                Kr[iy*3 + 2, ix*3 + 2] = Ke[6 + iy, 6 + ix]
        #
        #Печать данных
        """print('1:\nl1 =', l1, '\nl2 =', l2, '\nl3 =', l3)
        print('2:\nl12 =', l12, '\nl23 =', l23, '\nl31 =', l31)
        print('3:\ngab =\n', gab)
        print('4:\ng =', g)
        print('5:\nA =', A)
        print('6:\nB =\n', B, '\nD =\n', D, '\nBcap =\n', Bcap, '\nDcap =\n', Dcap)
        print('7:\nC =', C)
        print('8:\nKuu =\n', Kuu, '\nKufi =\n', Kufi, '\nKfifi =\n', Kfifi)
        print('9:\nKe =\n', Ke.astype(int))
        print('10:\nKr =\n', Kr.astype(int))"""
        #И возвращаем Kr
        return Kr


    def define_D_ep(self, sigma_x, sigma_y, tau_xy):
        """присвоение упругопластической матрицы жесткости
            КЭ работающему в пластике"""
        self.n_vect = n_vect(sigma_x, sigma_y, tau_xy)
        self.D_ep = self.D-(self.D @ self.n_vect.reshape(3,1) @
                    self.n_vect @ self.D)/\
                    (self.n_vect @ self.D @ self.n_vect.reshape(3,1)+
                    self.Et/(1-self.Et/self.E))


    def define_d_strains(self):
        """матрица деформаций на данном шаге итераций (3х1)"""
        e = self.matrix_B()[0] @ \
        np.array([[self.pnt[0].displ_x,
                self.pnt[0].displ_y,
                self.pnt[1].displ_x,
                self.pnt[1].displ_y,
                self.pnt[2].displ_x,
                self.pnt[2].displ_y]]).transpose()
        self.d_strain = e


    def define_strains(self):
        """матрица полных деформаций (3х1)"""
        e = self.matrix_B()[0] @ \
        np.array([[self.pnt[0].tot_displ_x,
                self.pnt[0].tot_displ_y,
                self.pnt[1].tot_displ_x,
                self.pnt[1].tot_displ_y,
                self.pnt[2].tot_displ_x,
                self.pnt[2].tot_displ_y]]).transpose()
        self.strain = e


    def define_stress(self):

        if len(self.pnt) == 3:
            '''Вычисление напряжений для упругой работы'''
            # определение геометрических параметров
            B = self.matrix_B()[0]
            # ~ вычисление приращения напряжений в линейной постановке
            if self.pl != True:
                D = self.D
                d_sigma = D @ B @ np.array(
                    [[self.pnt[0].displ_x,
                    self.pnt[0].displ_y,
                    self.pnt[1].displ_x,
                    self.pnt[1].displ_y,
                    self.pnt[2].displ_x,
                    self.pnt[2].displ_y]]).reshape(6,1)
                self.d_sigma_x = d_sigma[0,0]
                self.d_sigma_y = d_sigma[1,0]
                self.d_tau_xy = d_sigma[2,0]
                sx = self.d_sigma_x + self.sigma_x
                sy = self.d_sigma_y + self.sigma_y
                txy = self.d_tau_xy + self.tau_xy

                print('Напряжения элемента №', self.num, ':\nx = ', sx, '\ny = ', sy, '\ntxy = ', txy)
                
                # проверка условия пластичности
                if sigma_eq(sx, sy, txy) >= self.sigma_pr:
                    self.add_plast()
                    self.define_D_ep(sx, sy, txy)

            if self.pl == True:

                d_sigma = self.D @\
                        B @ np.array(
                        [[self.pnt[0].displ_x,
                        self.pnt[0].displ_y,
                        self.pnt[1].displ_x,
                        self.pnt[1].displ_y,
                        self.pnt[2].displ_x,
                        self.pnt[2].displ_y]]
                        ).reshape(6,1)
                n = 100
                for it in range(n):
                    N = n_vect(d_sigma[0,0]+self.sigma_x,
                                d_sigma[1,0]+self.sigma_y,
                                d_sigma[2,0]+self.tau_xy)
                    d_sigma -= (sigma_eq(d_sigma[0,0]+self.sigma_x,
                                        d_sigma[1,0]+self.sigma_y,
                                        d_sigma[2,0]+self.tau_xy)
                    -self.sigma_pr)/(N @ self.D @ N.reshape(3,1)) *\
                    self.D @ N.reshape(3,1)
                    if abs(sigma_eq(d_sigma[0,0]+self.sigma_x,
                                    d_sigma[1,0]+self.sigma_y,
                                    d_sigma[2,0]+self.tau_xy)
                    - self.sigma_pr) < 10:
                        #~ print('sigma break from: ', it)
                        break
                if it == 99: print('so much')

                sx = d_sigma[0,0] + self.sigma_x
                sy = d_sigma[1,0] + self.sigma_y
                txy = d_sigma[2,0] + self.tau_xy
                self.define_D_ep(sx, sy, txy)
    

            
                

        if len(self.pnt) == 4:
            pass
        # TODO: Дописать для 4-хузловых КЭ

        return sx, sy, txy

    #Определение напряжений для второго типа КЭ
    def define_stress_2(self):
        if len(self.pnt) == 3:
            '''Вычисление напряжений для упругой работы'''
            # определение геометрических параметров
            B = self.matrix_B()[0]
            # ~ вычисление приращения напряжений в линейной постановке
            if self.pl != True:
                D = self.D
                d_sigma = D @ B @ np.array(
                    [[self.pnt[0].displ_x,
                    self.pnt[0].displ_y,
                    self.pnt[1].displ_x,
                    self.pnt[1].displ_y,
                    self.pnt[2].displ_x,
                    self.pnt[2].displ_y]]).reshape(6,1)
                self.d_sigma_x = d_sigma[0,0]
                self.d_sigma_y = d_sigma[1,0]
                self.d_tau_xy = d_sigma[2,0]
                sx = self.d_sigma_x + self.sigma_x
                sy = self.d_sigma_y + self.sigma_y
                txy = self.d_tau_xy + self.tau_xy

                print('Напряжения элемента №', self.num, ':\nx = ', sx, '\ny = ', sy, '\ntxy = ', txy)
                
                # проверка условия пластичности
                if sigma_eq(sx, sy, txy) >= self.sigma_pr:
                    self.add_plast()
                    self.define_D_ep(sx, sy, txy)
    
    def add_stress(self, sigma_x, sigma_y, tau_xy):
        """ add stresses to FE """
        self.sigma_x = sigma_x
        self.sigma_y = sigma_y
        self.tau_xy = tau_xy

        self.sigma_1 = (self.sigma_x + self.sigma_y)/2 + np.sqrt(
                        ((self.sigma_x - self.sigma_y)/2)**2 +
                        self.tau_xy**2)
        self.sigma_2 = (self.sigma_x + self.sigma_y)/2 - np.sqrt(
                        ((self.sigma_x - self.sigma_y)/2)**2 +
                        self.tau_xy**2)

    def add_plast(self):
        self.pl = True

class Fe_2:
    """Класс, описывающий поведение стержневых элементов"""
    def __init__(self, num, types, n_stiff, *pnt):
        """присвоение характеристик при работе скрипта liraparser.py"""
        """ОБЩИЕ СВОЙСТВА"""
        """Номер КЭ"""
        self.num = num
        """Тип КЭ. В нашем случае - 2"""
        self.types = types
        """Тип жесткости КЭ"""
        self.n_stiff = n_stiff
        """Кортеж узлов стержня"""
        self.pnt = pnt
        """СВОЙСТВА СТЕРЖНЕЙ"""
        """Для учёта сдвига"""
        self.k = 1.2
        self.mu = 0.3
        """ВЫЧИСЛЯЕМЫЕ СВОЙСТВА"""
        """Длина стержня"""
        self.L = np.sqrt((self.pnt[1].x - self.pnt[0].x)**2\
                         + (self.pnt[1].y - self.pnt[0].y)**2)
        """Угол наклона стержня
        в будущем он пригодится только для вычисления матрицы поворота стержня"""
        self.phi = np.arctan2(self.pnt[1].y - self.pnt[0].y, self.pnt[1].x - self.pnt[0].x)
        """Временные переменные для вычисления матрицы поворота стержня
        Синус и косинус угла поворота"""
        Ts = np.sin(self.phi)
        Tc = np.cos(self.phi)
        """Матрица поворота стержня"""
        self.Q = np.array([[Tc, Ts, 0, 0, 0, 0],
                           [-Ts, Tc, 0, 0, 0, 0],
                           [0, 0, 1, 0, 0, 0],
                           [0, 0, 0, Tc, Ts, 0],
                           [0, 0, 0, -Ts, Tc, 0],
                           [0, 0, 0, 0, 0, 1]])
    
    def add_proporties(self, E, b, h):
        """присвоение характеристик при работе скрипта liraparser.py"""
        self.E = E
        #print('E стержня #', self.num, ' = ', self.E , 'т/м2')
        """/100 - костыли-костылики
        преобразование сантиметров в метры"""
        self.b = b/100
        self.h = h/100
        """А дальше вычисляем зависимые параметры"""
        """Площадь поперечного сечения стержня"""
        self.A = self.b*self.h
        """Момент инерции поперечного сечения стержня (по оси Y)"""
        self.Iy = (self.b*self.h**3)/12
        """Модуль сдвига"""
        self.G = self.E/(2*(1 + self.mu))
        
        #print('L = ', self.L, '\n\nQ = \n', self.Q, '\n\nE = ', self.E, '\nb = ', self.b, '\nh = ', self.h, '\nA = ', self.A, '\nIy = ', self.Iy, '\nG = ', self.G)
    
    def matrix_K(self):
        """Функция, возвращающая матрицу жесткости элемента
        с учётом сдвиговых деформаций"""
        """Дважды встречающийся делитель"""
        kk1 = self.A*self.G*self.k*self.L**3 + 12*self.E*self.L*self.Iy
        
        """Временные переменные для матрицы жесткости - для легкости заполнения
        Нумерация по порядку появления"""
        k1 = (self.A*self.E)/self.L
        k2 = (12*self.A*self.E*self.G*self.k*self.Iy)/kk1
        k3 = (6*self.A*self.E*self.G*self.k*self.Iy)\
            /(self.A*self.G*self.k*self.L**2 + 12*self.E*self.Iy)
        k4 = (4*self.E*self.Iy*(self.A*self.G*self.k*self.L**2 + 3*self.E*self.Iy))/kk1
        k5 = (2*self.E*self.Iy*(self.A*self.G*self.k*self.L**2 - 6*self.E*self.Iy))/kk1
        
        """Задание матрицы жёсткости"""
        """K = np.array([[k1, 0, 0, -k1, 0, 0],
                      [0, k2, k3, 0, -k2, k3],
                      [0, k3, k4, 0, -k3, k5],
                      [-k1, 0, 0, k1, 0, 0],
                      [0, -k2, -k3, 0, k2, -k3],
                      [0, k3, k5, 0, -k3, k4]])"""
        
        K = np.array([[k1, 0, 0, -k1, 0, 0],
                      [0, k2, k3, 0, -k2, k3],
                      [0, k3, k4, 0, -k3, k5],
                      [-k1, 0, 0, k1, 0, 0],
                      [0, -k2, -k3, 0, k2, -k3],
                      [0, k3, k5, 0, -k3, k4]])
        
        """Повёрнутая матрица жесткости"""
        TK = np.transpose(self.Q) @ K @ self.Q
        
        return TK
    
    """U - вектор узловых перемещений системы"""
    def define_stress(self, U):
        """Получаем матрицу жёсткости элемента"""
        K = self.matrix_K()
        
        """Вектор перемещения"""
        UE = np.zeros(6)
        
        """Средние значения - N, Q и M элементов"""
        FEa = np.zeros(3)
        
        """Получаем перемещения узлов - нулевого и первого"""
        UE[0:3] = U[3*(self.pnt[0].num-1):3*self.pnt[0].num]
        UE[3:6] = U[3*(self.pnt[1].num-1):3*self.pnt[1].num]
        
        """Переходим к локальным координатам
        Заполняем вектор внутренних сил элемента"""
        FE = K @ (self.Q @ UE)
        
        """И вычисляем средние значения для элементов"""
        for i in range(3):
            FEa[i] = (-FE[i] + FE[i + 3])/2
        
        return FEa[0], FEa[1], FEa[2]
    
    def define_stress_2(self, U):
        """Получаем матрицу жёсткости элемента"""
        K = self.matrix_K()
        
        """Вектор перемещения"""
        UE = np.zeros(6)
        
        """Средние значения - N, Q и M элементов"""
        FEa = np.zeros(3)
        
        """Получаем перемещения узлов - нулевого и первого"""
        UE[0:3] = U[3*(self.pnt[0].num-1):3*self.pnt[0].num]
        UE[3:6] = U[3*(self.pnt[1].num-1):3*self.pnt[1].num]
        
        """Переходим к локальным координатам
        Заполняем вектор внутренних сил элемента"""
        FE = K @ (self.Q @ UE)
        
        """И вычисляем средние значения для элементов"""
        for i in range(3):
            FEa[i] = (-FE[i] + FE[i + 3])/2
        
        return FEa[0], FEa[1], FEa[2]

    def add_stress(self, N, Q, M):
        """ add stresses to FE """
        self.sigma_x = N
        self.sigma_y = Q
        self.tau_xy = M


class Point:
    """Класс описывающий поведение узлов"""
    def __init__(self, num, x, y):
        self.num = num
        self.x = x
        self.y = y
        self.tot_displ_x = 0 # полное перемещение от внешних сил
        self.tot_displ_y = 0
        self.displ_x = 0 # перемещение на последнем шаге итерации
        self.displ_y = 0
        self.bc = None

    def add_bound_cond(self, *bc):
        self.bc = bc

    def summ_displ(self):
        self.tot_displ_x += self.displ_x
        self.tot_displ_y += self.displ_y


class Timer(object):
    """Таймер- считает время выполнения блока кода"""
    def __enter__(self):
        self._startTime = time.time()

    def __exit__(self, type, value, traceback):
        print("Время вычисления матрицы: {:.3f} сек.".format(
            time.time() - self._startTime))


def sigma_eq(sx, sy, tau):
    if liraparser.TASK == 'plane stress':
        return np.sqrt(sx**2-sx*sy+sy**2+3*tau**2)
    if liraparser.TASK == 'plane strain':
        return np.sqrt(3/4*(sx-sy)**2+3*tau**2)


def n_vect(dsx, dsy, dtau):
    if liraparser.TASK == 'plane stress':
        return np.array([[
            ((2*dsx - dsy)/(2*np.sqrt(dsx**2-
            dsx*dsy + dsy**2+3*dtau))),
            ((2*dsy - dsx)/(2*np.sqrt(dsx**2-
            dsx*dsy + dsy**2+3*dtau))),
            ((6*dtau)/(2*np.sqrt(dsx**2-
            dsx*dsy + dsy**2+3*dtau)))]]).reshape(1,3)
    if liraparser.TASK == 'plane strain':
        return np.array([[
            3*(dsx-dsy)/4/np.sqrt(3/4*(dsx-dsy)**2+3*dtau**2),
            -3*(dsx-dsy)/4/np.sqrt(3/4*(dsx-dsy)**2+3*dtau**2),
            6*dtau/4/np.sqrt(3/4*(dsx-dsy)**2+3*dtau**2)]])


def make_q(le, lp):
    length = len(lp)
    q = np.zeros(length*2).reshape(length*2,1)
    for elm in le:
        B, area = elm.matrix_B()
        q_elm = B.transpose() @\
            np.array(
                [elm.sigma_x + elm.d_sigma_x,
                elm.sigma_y + elm.d_sigma_y,
                elm.tau_xy + elm.d_tau_xy])*area*elm.h
        for k, point in enumerate(elm.pnt):
            assert np.shape(q[2*(point.num-1):2*point.num]) == (2,1) # нарушена размерность
            assert np.shape(q_elm.reshape(6,1)[k*2:2*(k+1)]) == (2,1) # нарушена размерность
            q[2*(point.num-1):2*point.num] += q_elm.reshape(6,1)[k*2:2*(k+1)]
    for point in lp:
        if point.bc == (1,): q[2*(point.num-1)] = 0
        if point.bc == (3,): q[2*point.num-1] = 0
        if point.bc == (1,3):
            q[2*(point.num-1)] = 0
            q[2*point.num-1] = 0 
    return q
