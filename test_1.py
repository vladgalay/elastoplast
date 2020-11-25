#Тестовая программа по греческой статье
import math
import numpy as np
import matplotlib.pyplot as plt
#Исходные данные
#Модуль сдвига
G = 2039000/(2*(1 + 0.25))
#Коэффициент Пуассона
nu = 0.3
#Константы
l = 0.0
a = 0.25
#Массив координат точек пластины
pc = np.array([[0, 0],
               [1, 0],
               [0, 1],
               [1, 1]])
#Массив элементов - в нём номера точек, на которых он строится
elm = np.array([[0, 1, 2],
                [1, 2, 3]])
#Вектор узловых нагрузок и узловых моментов
"""fu = np.array([0, 0, 0, 0, 0, 0])
ffi = np.array([1e-12, 0, 0])"""
P = np.array([0, 0, 0, 0, -200, 0, 0, 0, 0, 0, 0, 0])

#Функция получения матрицы жёсткости элемента
#pn - кортеж с номерами точек пластины
def matrix_K(pn):
    #1. Определение длин сторон пластины
    #Номер стороны = номер противоположного узла пластины
    l1 = math.sqrt((pc[pn[1], 0] - pc[pn[2], 0])**2 + (pc[pn[1], 1] - pc[pn[2], 1])**2)
    l2 = math.sqrt((pc[pn[0], 0] - pc[pn[2], 0])**2 + (pc[pn[0], 1] - pc[pn[2], 1])**2)
    l3 = math.sqrt((pc[pn[1], 0] - pc[pn[0], 0])**2 + (pc[pn[1], 1] - pc[pn[0], 1])**2)
    print('1:\nl1 =', l1, '\nl2 =', l2, '\nl3 =', l3)
    #2. Определение дополнительных параметров, связанных с размером пластины
    l12 = l1**2 + l2**2 - l3**2
    l23 = l2**2 + l3**2 - l1**2
    l31 = l3**2 + l1**2 - l2**2
    print('2:\nl12 =', l12, '\nl23 =', l23, '\nl31 =', l31)
    #3. Определение метрического тензора (что это?)
    gab = np.array([[l2**2, l12/2],
                    [l12/2, l1**2]])
    print('3:\ngab =\n', gab)
    #4. Вычисление определителя метрического тензора
    g = np.linalg.det(gab)
    print('4:\ng =', g)
    #5. Вычисление площади пластины через определитель метрического тензора
    A = math.sqrt(g/4)
    print('5:\nA =', A)
    #6. Задаём известные матрицы
    B = np.array([[1, 0, 0, 0, -1, 0],
                  [0, 0, 0, 1, 0, -1],
                  [0, 1, 0, 0, 0, -1],
                  [0, 0, 1, 0, -1, 0]])
    D = G*np.array([[(2*(1 - nu))/(1 - 2*nu), (2*nu)/(1 - 2*nu), 0, 0],
                    [(2*nu)/(1 - 2*nu), (2*(1 - nu))/(1 - 2*nu), 0, 0],
                    [0, 0, 1 + a, 1 - a],
                    [0, 0, 1 - a, 1 + a]])
    Bcap = np.array([[1, 0, -1],
                     [0, 1, -1]])
    Dcap = (4*G*l**2)*np.array([[1, 0],
                                [0, 1]])
    print('6:\nB =\n', B, '\nD =\n', D, '\nBcap =\n', Bcap, '\nDcap =\n', Dcap)
    #7. Задаём константу
    C = 2*a
    print('7:\nC =', C)
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
    print('8:\nKuu =\n', Kuu, '\nKufi =\n', Kufi, '\nKfifi =\n', Kfifi)
    #9. Составляем матрицу жесткости элемента
    Ke = np.zeros((9, 9))
    Ke[0:6, 0:6] = Kuu
    Ke[0:6, 6:9] = Kufi
    Ke[6:9, 0:6] = Kufi.transpose()
    Ke[6:9, 6:9] = Kfifi
    print('9:\nKe =\n', Ke.astype(int))
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
    print('10:\nKr =\n', Kr.astype(int))
    #И возвращаем Kr
    return Kr

#Функция установки граничного условия
#KG - глобальная матрица жёсткости
#n - номер столбца и строки, который мы закрепляем (начинается с 0)
def set_BC(KG, n):
    #Размер матрицы
    size = KG.shape[0]
    KG[0:size, n:n + 1] = np.zeros((size, 1))
    KG[n:n + 1, 0:size] = np.zeros(size).transpose()
    KG[n, n] = 1
    
    return KG

#Функция добавления матрицы к глобальной матрице жёсткости
#KG - глобальная матрица жёсткости
#Kelm - матрица жёсткости элемента
#n - номер элемента
def matrix_Add(KG, Kelm, n):
    for iy in range(3):
        for ix in range(3):
            KG[3*elm[n, iy]:3*elm[n, iy] + 3, 3*elm[n, ix]:3*elm[n, ix] + 3] += Kelm[iy*3:iy*3 + 3, ix*3:ix*3 + 3]
    
    return KG
#Количество точек
pn = pc.shape[0]
#И количество элементов
en = elm.shape[0]
#Получаем глобальную матрицу жёсткости
K = np.zeros((pn*3, pn*3))
#Заносим в матрицу жёсткости элементы
for i in range(en):
    K = matrix_Add(K, matrix_K(elm[i]), i)
#И присваиваем граничное условие
K = set_BC(K, 0)
K = set_BC(K, 1)
K = set_BC(K, 3)
K = set_BC(K, 4)
print('10:\nBC =\n', K.astype(int))
#10. Составляем вектор нагрузок
"""P = np.zeros(9)
P[0:6] = fu
P[6:9] = ffi"""
print('10:\nP =\n', P)
#11. Получаем вектор узловых перемещений U, решив систему линейных уравнений
U = np.linalg.solve(K, P.transpose())
print('11:\nU =\n', U)

#Графический вывод
pcd = pc.transpose()

plt.figure(figsize = (16, 9))
#Вывод исходной схемы
#Проходимся по элементам
for i in range(elm.shape[0]):
    #Массив точек треугольника
    #Четвёртая точка - чтобы замкнуть
    #Сразу транспонирован
    tp = np.zeros((2, 4))
    #Заполняем по точкам
    for j in range(3):
        for ii in range(2):
            tp[ii, j] = pc[elm[i][j]][ii]
    #Дублируем последнюю точку
    tp[0:2, 3:4] = tp[0:2, 0:1]
    #И рисуем
    plt.plot(tp[0], tp[1], 'bo-')
#Вывод схемы с перемещениями
#Проходимся по элементам
for i in range(elm.shape[0]):
    #Массив точек треугольника
    #Четвёртая точка - чтобы замкнуть
    #Сразу транспонирован
    tp = np.zeros((2, 4))
    #Заполняем по точкам
    for j in range(3):
        for ii in range(2):
            tp[ii, j] = pc[elm[i][j]][ii] + U[elm[i][j]*3 + ii]
    #Дублируем последнюю точку
    tp[0:2, 3:4] = tp[0:2, 0:1]
    #И рисуем
    plt.plot(tp[0], tp[1], 'ro-')
plt.axis('equal')
plt.title('Общий вид системы')
plt.show()