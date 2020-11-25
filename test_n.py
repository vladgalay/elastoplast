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
a = 0.0
#Массив координат точек пластины
pc = np.array([[0, 0],
               [1, 0],
               [0, 1]])
pc2 = np.array([[0, 0],
                [0.5, 0],
                [1, 0],
                [0.5, 0.5],
                [0, 1],
                [0, 0.5]])
#Массив элементов - в нём номера точек, на которых он строится
elm = np.array([0, 1, 2])
#Вектор узловых нагрузок и узловых моментов
"""fu = np.array([0, 0, 0, 0, 0, 0])
ffi = np.array([1e-12, 0, 0])"""
P1 = np.array([0, 0, 0, 0, 1000, 0])
P2 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1000, 0, 0, 0])

#Функция получения матрицы жёсткости элемента
#pn - кортеж с номерами точек пластины
#te - тип матрицы:
#1 - 6x6
#2 - 12x12
def matrix_K(pn, te):
    #1. Определение длин сторон пластины
    #Номер стороны = номер противоположного узла пластины
    l1 = math.sqrt((pc[pn[1], 0] - pc[pn[2], 0])**2 + (pc[pn[1], 1] - pc[pn[2], 1])**2)
    l2 = math.sqrt((pc[pn[0], 0] - pc[pn[2], 0])**2 + (pc[pn[0], 1] - pc[pn[2], 1])**2)
    l3 = math.sqrt((pc[pn[1], 0] - pc[pn[0], 0])**2 + (pc[pn[1], 1] - pc[pn[0], 1])**2)
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
    D = G*np.array([[(2*(1 - nu))/(1 - 2*nu), (2*nu)/(1 - 2*nu), 0, 0],
                    [(2*nu)/(1 - 2*nu), (2*(1 - nu))/(1 - 2*nu), 0, 0],
                    [0, 0, 1 + a, 1 - a],
                    [0, 0, 1 - a, 1 + a]])
    #7. Задаём константу
    C = 2*a
    #8. Вычисляем подматрицы
    Kuu = A*(B.transpose() @ D @ B)
    #9. Составляем матрицу жесткости элемента
    Ke = Kuu
    if te == 2:
        E1 = D[0][0] + D[0][3]
        E2 = D[0][1] + D[1][3]
        E3 = D[0][2] + D[2][3]
        E4 = D[0][3] + D[3][3]
        
        M1 = D[0][1] + D[0][2]
        M2 = D[1][1] + D[1][2]
        M3 = D[1][2] + D[2][2]
        M4 = D[1][3] + D[2][3]
        
        R1 = D[0][1] + D[2][3]
        R2 = D[0][2] + D[1][3]
        
        T1 = 8*(D[0][0] + D[0][3] + D[3][3])
        T2 = 8*(D[1][1] + D[1][2] + D[2][2])
        
        #Заполняем одну половину, другую отзеркалим
        Ke = np.array([[3*D[0][0], 3*D[0][2], -4*D[0][3], -D[0][1], E1, M1, 0, 0, -4*E1, -4*M1, 4*D[0][3], 4*D[0][1]],
        [0, 3*D[2][2], -D[2][3], -D[1][2], E3, M3, 0, 0, -4*E3, -4*M3, 4*D[2][3], 4*D[1][2]],
        [0, 0, 3*D[3][3], 3*D[1][3], E4, M4, -4*E4, -4*M4, 0, 0, 4*D[0][3], 4*D[2][3]],
        [0, 0, 0, 3*D[1][1], E2, M2, -4*E2, -4*M2, 0, 0, 4*D[0][1], 4*D[1][2]],
        [0, 0, 0, 0, 3*(E1 + E4), 3*(M1 + M4), -4*E4, -4*E2, -4*E1, -4*E3, 0, 0],
        [0, 0, 0, 0, 0, 3*(M2 + M3), -4*M4, -4*M2, -4*M1, -4*M3, 0, 0],
        [0, 0, 0, 0, 0, 0, T1, 4*R1 + 8*R2, 8*D[0][3], 4*R1, -8*E1, -4*(E3 + M1)],
        [0, 0, 0, 0, 0, 0, 0, T2, 4*R1, 8*D[1][2], -4*(E3 + M1), -8*M3],
        [0, 0, 0, 0, 0, 0, 0, 0, T1, 4*R1 + 8*R2, -8*E4, -4*(E2 + M4)],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, T2, -4*(E2 + M4), -8*M2],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, T1, 4*R1 + 8*R2],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, T2]])
        
        #Зеркалим незаполненную половину матрицы
        for i in range(12):
            for j in range(i):
                Ke[i][j] = Ke[j][i]
    #Секция вывода данных
    """print('1:\nl1 =', l1, '\nl2 =', l2, '\nl3 =', l3)
    print('2:\nl12 =', l12, '\nl23 =', l23, '\nl31 =', l31)
    print('3:\ngab =\n', gab)
    print('4:\ng =', g)
    print('5:\nA =', A)
    print('6:\nB =\n', B, '\nD =\n', D, '\nBcap =\n', Bcap, '\nDcap =\n', Dcap)
    print('7:\nC =', C)
    print('8:\nKuu =\n', Kuu, '\nKufi =\n', Kufi, '\nKfifi =\n', Kfifi)
    print('9:\nKe =\n', Ke.astype(int))"""
    #И возвращаем Ke
    return Ke

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

#Получаем матрицы жёсткости элементов
K1 = matrix_K(elm, 1)
K2 = matrix_K(elm, 2)
#Заносим в матрицу жёсткости элементы
#И присваиваем граничное условие
K1 = set_BC(K1, 0)
K1 = set_BC(K1, 1)
K1 = set_BC(K1, 3)
#K1 = set_BC(K1, 4)

K2 = set_BC(K2, 0)
K2 = set_BC(K2, 1)
K2 = set_BC(K2, 5)
#K2 = set_BC(K2, 6)
#print('10:\nBC =\n', K.astype(int))
#10. Составляем вектор нагрузок
"""P = np.zeros(9)
P[0:6] = fu
P[6:9] = ffi"""
#print('10:\nP =\n', P)
#11. Получаем вектор узловых перемещений U, решив систему линейных уравнений
U1 = np.linalg.solve(K1, P1.transpose())
U2 = np.linalg.solve(K2, P2.transpose())
U1 = U1.reshape((int(U1.shape[0]/2), 2)).transpose()
U2 = U2.reshape((int(U2.shape[0]/2), 2)).transpose()
#print('11:\nU =\n', U1)
print('U1:\n', U1, '\n\nU2:\n', U2)

#Графический вывод
k = 1000
pcd = pc.transpose()
pc2d = pc2.transpose()
print('pcd:\n', pcd, '\n\npc2d:\n', pc2d)

plt.figure(figsize = (16, 9))

plt.plot(pc2d[0], pc2d[1], 'bo-')
plt.plot(pcd[0] + U1[0]*k, pcd[1] + U1[1]*k, 'ro-')
plt.plot(pc2d[0] + U2[0]*k, pc2d[1] + U2[1]*k, 'go-')
#Вывод исходной схемы
#Проходимся по элементам
"""for i in range(elm.shape[0]):
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
    plt.plot(tp[0], tp[1], 'ro-')"""
plt.axis('equal')
plt.title('Общий вид системы')
plt.show()