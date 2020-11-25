#! /usr/bin/env python
# -*- coding: utf-8 -*-

#Библиотека вывода графика
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

fig, ax = plt.subplots(figsize = (16, 9))

#Подготовка к рисованию
#def PrepareGraph():
    #Задаём размеры графика
    #fig = plt.figure(figsize = (16, 9))

#Добавление графика в вывод
#le - список элементов
#Ps - вектор координат узлов в формате x1, y1, 0, x2, y2, 0...
#0 - чтобы легко складывать с вектором узловых перемещений
#U - вектор узловых перемещений
#k - коэффициент умножения перемещений
#style - стиль вывода графика
def AddGraph(le, Ps, U, k, style):
    #Складываем с учётом коэффициента
    Pf = Ps + U*k
    #Далее проходимся по элементам
    for elm in le:
        #Если это стержень
        if len(elm.pnt) == 2:
            #Рисуем линию
            #Массив точек линии
            pa = np.zeros((2, 2))
            #Заполняем по точкам
            for j in range(2):
                for ii in range(2):
                    pa[ii, j] = Pf[(elm.pnt[j].num - 1)*3 + ii]
            #И рисуем
            ax.plot(pa[0], pa[1], style, linewidth = 5.0)
        
        #И если это стержень
        if len(elm.pnt) == 3:
            #Рисуем линию
            #Массив точек линии
            pa = np.zeros((2, 4))
            #Заполняем по точкам
            for j in range(3):
                for ii in range(2):
                    pa[ii, j] = Pf[(elm.pnt[j].num - 1)*3 + ii]
            #Дублируем последнюю точку
            pa[0:2, 3:4] = pa[0:2, 0:1]
            #И рисуем
            ax.plot(pa[0], pa[1], style, linewidth = 1.0)

#Добавление изополей в график
#le - список элементов
#Ps - вектор координат узлов в формате x1, y1, 0, x2, y2, 0...
#0 - чтобы легко складывать с вектором узловых перемещений
#U - вектор узловых перемещений
#Ct - одномерный массив значений изополей по точкам
def AddContour(le, Ps, U, k, Ct):
    #Складываем с учётом коэффициента
    Pf = Ps + U*k
    
    #Формируем списки координат узлов по осям
    #А точнее переформировываем вектор координат и транспонируем его
    Pf = Pf.reshape((int(Pf.shape[0]/3), 3)).transpose()
    print(Ct)
    
    #Количество треугольников
    #Лютый бред их так считать, нужно потом сделать ввод вне функции
    nt = 0
    
    for elm in le:
        if len(elm.pnt) == 3:
            nt += 1
    
    #Теперь немного рутинной работы
    #Получаем список точек треугольных пластин в виде массива [[p1, p2, p3], [p1, p2, p3], ...]
    #Не по кол-ву элементов - потому что среди элементов могут быть и стержни
    lt = np.zeros((nt, 3))
    num = 0
    
    for elm in le:
        if len(elm.pnt) == 3:
            lt[num] = np.array([[elm.pnt[0].num - 1, elm.pnt[1].num - 1, elm.pnt[2].num - 1]])
            num += 1
    
    np.delete(lt, 0, 0)
    
    #Проводим триангуляцию точек
    print('ДО ТРИАНГУЛЯЦИИ')
    print(lt)
    triang = mtri.Triangulation(Pf[0], Pf[1], lt)
    print('ПОСЛЕ ТРИАНГУЛЯЦИИ')
    
    #И выводим график
    #plt.triplot(triang, 'ko-')
    cntr = ax.tricontourf(triang, Ct, cmap='plasma')
    
    fig.colorbar(cntr, ax = ax)
    

#Вывод результатов рисования
def ShowGraph():
    ax.axis('equal')
    ax.set_title('Общий вид системы')
    plt.show()