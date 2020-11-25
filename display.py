import liraparser
import classes
import graph_plot as gp
import nonlinear
import numpy as np
import graph

if __name__ == '__main__':
	le, lp = liraparser.list_of_elem()
	loads = liraparser.list_of_loads()

	#Графический вывод
	#graph.PrepareGraph()
	#Количество точек
	npo = 0
	for point in lp:
		npo += 1
	npo *= 3
	#Подготавливаем вектор координат узлов в нужном формате
	Pc = np.zeros(npo)
	#Заполняем его
	for point in lp:
		Pc[(point.num - 1)*3 + 0] = point.x
		Pc[(point.num - 1)*3 + 1] = point.y
	#И рисуем
	#Коэффициент увеличения
	kk = 1
	#Оригинал
	graph.AddGraph(le, Pc, np.zeros(npo), kk, 'bo-')
	graph.ShowGraph()