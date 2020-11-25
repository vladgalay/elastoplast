#Анализ влияния констант на квадратичное расхождение от классического расчёта
import liraparser
import classes
import graph_plot as gp
import nonlinear
import numpy as np
import math
import graph
import linear
import matplotlib.pyplot as plt

#Начальные и конечные значения констант
sl = 0.0*1
el = 1.0*1
sa = -1.1*1
ea = -1.0*1
#Шаг изменения значений констант
il = 0.1
ia = 0.01

#Функция, возвращающая вектор абсолютных величин перемещений узлов
#U - вектор перемещений в формате [X, Z, UY, ...]
def GetDisplacementVector(U):
    #Проводим необходимую трансформацию
    U = U.reshape((int(U.shape[0]/3), 3)).transpose()
    #И возвращаем вектор
    return np.sqrt(U[0]**2 + U[1]**2)

#Получение процента различия (превышения)
#D1 - то, что сравниваем
#D2 - то, с чем сравниваем
def GetDifference(D1, D2):
    #Проходимся по вектору D2 и удаляем нулевые элементы в обоих векторах
    #Потому что на ноль делить нельзя
    """for i in range(D2.shape[0] - 1, -1, -1):
        if D2[i] == 0:
            D1 = np.delete(D1, i)
            D2 = np.delete(D2, i)
    
    #Далее производим необходимую арифметику
    return np.mean((D1/D2)*100)"""
    return (np.sum(D2)/np.sum(D1))*100

if __name__ == '__main__':
    le, lp = liraparser.list_of_elem()
    loads = liraparser.list_of_loads()
    with classes.Timer() as p:
        #Начинаем
        #Образец перемещений схемы
        #Сначала считаем перемещения по классической теории эластичности
        Uc = linear.matrix_U(le, lp, loads, False, 0.1, 0.5)
        #Оттуда возвращаем вектор перемещений узлов
        UcD = GetDisplacementVector(Uc)
        
        #Далее считаем все остальные вариации
        x, y = np.meshgrid(np.arange(sa, ea, ia), np.arange(sl, el, il))
        #print(x)
        #print(y)
        #z = GetDifference(GetDisplacementVector(linear.matrix_U(le, lp, loads, True, x, y)), UcD)
        z = np.zeros_like(x)
        
        indx = 0
        indy = 0
        
        #Всего вычислений
        total = z.flatten().shape[0]
        #И сколько уже посчитано
        count = 0
        
        for ll in y.transpose()[0]:
            indx = 0
            
            for la in x[0]:
                count += 1
                print(count, 'из', total, ', l =', ll, ', a =', la)
                
                if ll == 0 and la == 0:
                    z[indy][indx] = GetDifference(GetDisplacementVector(linear.matrix_U(le, lp, loads, True, 1e-15, 1e-15)), UcD)
                else:
                    z[indy][indx] = GetDifference(GetDisplacementVector(linear.matrix_U(le, lp, loads, True, ll, la)), UcD)
                indx += 1
            
            indy += 1
        
        #print(z)
        
        fig, ax = plt.subplots()

        cntr = ax.contourf(x, y, z, cmap='plasma')
        
        fig.colorbar(cntr, ax = ax)
        
        plt.xlabel('a')
        plt.ylabel('l')
        
        ax.set_title('Система ' + liraparser.INPUT)
        
        plt.show()
        
        #U = linear.matrix_U(le, lp, loads, False)
        #dv = GetDisplacementVector(U)
        #print(dv)
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
        """Pc = np.zeros(U.shape[0])
        #Заполняем его
        for point in lp:
            Pc[(point.num - 1)*3 + 0] = point.x
            Pc[(point.num - 1)*3 + 1] = point.y
        #И рисуем
        #Коэффициент увеличения
        kk = 500
        #Оригинал
        graph.AddGraph(le, Pc, np.zeros_like(U), kk, 'b-')
        #С перемещениями пластин с двумя степенями свободы
        graph.AddGraph(le, Pc, U, kk, 'r-')
        #И с новым типом КЭ
        nU = linear.matrix_U(le, lp, loads, True)
        #graph.AddContour(le, Pc, nU, kk, nU.reshape((int(nU.shape[0]/3), 3)).transpose()[2])
        graph.AddGraph(le, Pc, nU, kk, 'g-')
        #Делаем массив векторов координат и моментов (0 - X, 1 - Z, 2 - uY)
        graph.ShowGraph()"""
