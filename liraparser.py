#! /usr/bin/env python
# -*- coding: utf-8 -*-
import classes
import re
import math
import numpy as np

#INPUT = '1_1.txt'
#INPUT = 'test1В.txt'
#INPUT = 'test1.txt'
#INPUT = 'З3.txt'
#INPUT = '5.txt'
#INPUT = 'Кирш 4Н8.txt'
INPUT = 'Test System.txt'
#TASK = 'plane stress' # choose 'plane strain' or 'plane stress'
TASK = 'plane strain' # choose 'plane strain' or 'plane stress'

def list_of_elem():
    """return list of elements from txt file in type:
    [<number of point>, [<coordinates of point>]]
    """
    points_list = []

    with open(INPUT, 'r', encoding='cp1251') as file:
        st = re.compile(r"^\( 4/", re.S)  # start of a block
        en = re.compile(r"(.*?) \)$", re.S)  # end of a block

        line = file.readline()
        while line:
            line = file.readline()
            if st.search(line):
                break
        if line == '':
            raise RuntimeError("Строчка с '( 4/' не найдена")

        del st
        count = 1
        while line:
            line = file.readline()
            if en.search(line):
                break
            else:
                for i in line.split('/'):
                    if i == '\n':
                        continue
                    k = i[:-1].split(' ')
                    k = [float(elem) for elem in k]
                    points_list.append(classes.Point(count, k[0], k[2]))
                    count += 1

        st = re.compile(r"^\( 5/", re.S)  # start of a block
        en = re.compile(r"(.*?) \)$", re.S)  # end of a block

        while line:
            line = file.readline()
            if st.search(line):
                break
        if line == '':
            raise RuntimeError("Строчка с '( 5/' не найдена")
        del st

        while line:
            line = file.readline()
            if en.search(line):
                break
            else:
                for i in line.split('/'):
                    if i == '\n':
                        continue
                    k = i[:-1].split(' ')
                    k = [int(elem) for elem in k]
                    points_list[k[0]-1].add_bound_cond(*k[1:])

    print('Список узлов сформирован.\nКоличество узлов - {}'.format(
        len(points_list)))

    """return list of elements from txt file in type:
        [<type of elements>, <number of stiffness>,*<number of points>]
        """
    list_el = []

    with open(INPUT, 'r', encoding='cp1251') as file:
        st = re.compile(r"^\( 1/", re.S)  # start of a block
        en = re.compile(r"(.*?) \)$", re.S)  # end of a block

        line = file.readline()
        while line:
            line = file.readline()
            if st.search(line):
                break
        if line == '':
            raise RuntimeError("Строчка с '( 1/' не найдена")

        del st
        count = 1
        while line:
            line = file.readline()
            if en.search(line):
                break
            else:
                for i in line.split('/'):
                    if i == '\n':
                        continue
                    k = i[:-1].split(' ')
                    k = [int(elem) for elem in k]
                    """Для стержней"""
                    """if len(k[2:]) == 2:
                        v = [k[2], k[3]]
                    if len(k[2:]) == 3:
                        v = [k[4], k[2], k[3]]
                    if len(k[2:]) == 4:
                        v = [k[2], k[3], k[5], k[4]]"""
                    if len(k[2:]) == 2:
                        v = [k[2], k[3]]
                    if len(k[2:]) == 3:
                        v = [k[2], k[3], k[4]]
                    if len(k[2:]) == 4:
                        v = [k[2], k[3], k[4], k[5]]
                    pnt = map(lambda x: points_list[x-1], v)
                    """Для стержней"""
                    if len(k[2:]) == 2:
                        list_el.append(classes.Fe_2(count, k[0], k[1], *pnt))
                        """Для пластин"""
                    else:
                        list_el.append(classes.Fe(count, k[0], k[1], *pnt))
                    
                    count += 1

        # ниже идет присвоение жесткости элементам

        st = re.compile(r"^\( 3/", re.S)  # start of a block
        en = re.compile(r"(.*?) \)$", re.S)  # end of a block
        list_of_stiff = []
        while line:
            line = file.readline()
            if st.search(line):
                break
        if line == '':
            raise RuntimeError("Строчка с '( 3/' не найдена")

        while line:
            line = file.readline()
            if en.search(line):
                break
            else:
                for i in line.split('/'):
                    if i == '\n':
                        continue
                    k = i[:-1].split()
                    
                    """Мой костыль - последняя единица не распознаётся"""
                    if len(k) == 4:
                        k = [int(k[0]), k[1], float(k[2]), float(k[3]), float(k[3])]
                    else:
                        k = [int(k[0]), k[1], float(k[2]), float(k[3]), float(k[4])]
                    list_of_stiff.append(k)

    for elm in list_el:
        prop = list_of_stiff[elm.n_stiff - 1]
        elm.add_proporties(prop[2], prop[3], prop[4])
        # print('{:3d} - E: {:6g}т/м2; V: {:4g}; h: {:4g}м'\
        #     .format(elm.num, elm.E, elm.V, elm.h))
    print('Список элементов сформирован.\nКоличество элементов - {}'
          .format(len(list_el)))
    return list_el, points_list


def list_of_loads():
    direct = []
    load = []
    loads = []
    with open(INPUT, 'r', encoding='cp1251') as file:
        st = re.compile(r"^\( 6/", re.S)  # start of a block
        en = re.compile(r"(.*?) \)$", re.S)  # end of a block
        line = file.readline()
        while line:
            line = file.readline()
            if st.search(line):
                break
        if line == '':
            raise RuntimeError("Строчка с '( 6/' не найдена")
        while line:
            line = file.readline()
            if en.search(line):
                break
            else:
                for i in line.split('/'):
                    if i == '\n':
                        continue
                    k = i[:-1].split()
                    k = [int(i) for i in k]
                    direct.append(k)

        st = re.compile(r"^\( 7/", re.S)  # start of a block
        en = re.compile(r"(.*?) \)$", re.S)  # end of a block
        while line:
            line = file.readline()
            if st.search(line):
                break
        if line == '':
            raise RuntimeError("Строчка с '( 7/' не найдена")
        while line:
            line = file.readline()
            if en.search(line):
                break
            else:
                for i in line.split('/'):
                    if i == '\n':
                        continue
                    k = i[:-1].split()
                    
                    k = [int(k[0]), float(k[1])]
                    load.append(k)
        for d in direct:
            if d[1] == 0:
                loads.append([d[0], d[2], load[d[3]-1][1]])
            """if d[1] == 19:
                num_of_pnt1 =elm[d[0]-1].pnt[]
                num_of_pnt2 =elm[d[0]-1].pnt
                drct=
                ld=
                loads.append(elm[d[0]-1].num)
                loads.append("""
                
    return loads


if __name__ == "__main__":
    le, lp = list_of_elem()
    loads = list_of_loads()
