#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Robótica Computacional - 
# Grado en Ingeniería Informática (Cuarto)
# Práctica: Resolución de la cinemática inversa mediante CCD
#           (Cyclic Coordinate Descent).

import sys
from math import *
import numpy as np
import matplotlib.pyplot as plt
import colorsys as cs

# ******************************************************************************
# Declaración de funciones

def muestra_origenes(O,final=0):
  # Muestra los orígenes de coordenadas para cada articulación
  print('Origenes de coordenadas:')
  for i in range(len(O)):
    print('(O'+str(i)+')0\t= '+str([round(j,3) for j in O[i]]))
  if final:
    print('E.Final = '+str([round(j,3) for j in final]))

def muestra_robot(O,obj):
  # Muestra el robot graficamente
  plt.figure()
  plt.xlim(-L,L)
  plt.ylim(-L,L)
  T = [np.array(o).T.tolist() for o in O]
  for i in range(len(T)):
    plt.plot(T[i][0], T[i][1], '-o', color=cs.hsv_to_rgb(i/float(len(T)),1,1))
  plt.plot(obj[0], obj[1], '*')
  plt.pause(0.0001)
  # plt.show()
  plt.draw()
  plt.waitforbuttonpress(0)
  
#  input()
  plt.close()

def matriz_T(d,th,a,al):
  # Los ángulos deben estar en Radianes
  return [[cos(th), -sin(th)*cos(al),  sin(th)*sin(al), a*cos(th)]
         ,[sin(th),  cos(th)*cos(al), -sin(al)*cos(th), a*sin(th)]
         ,[      0,          sin(al),          cos(al),         d]
         ,[      0,                0,                0,         1]
         ]

# Devuelve una lista de coordenadas de (x,y) de cada uno de los puntos o
# respecto al referencial inicial
def cin_dir(th,a):
  #Sea 'th' el vector de thetas
  #Sea 'a'  el vector de longitudes
  T = np.identity(4)
  o = [[0,0]]
  for i in range(len(th)):
    T = np.dot(T,matriz_T(0,th[i],a[i],0))
    tmp=np.dot(T,[0,0,0,1])
    o.append([tmp[0],tmp[1]])
  return o

# ******************************************************************************
# Cálculo de la cinemática inversa de forma iterativa por el método CCD

# valores articulares arbitrarios para la cinemática directa inicial
th=[0.,0.,0.]
a =[5.,0.,5.]
prismatica = [False, True, False]
qMin = np.array([-30 * pi / 180, 0, -30 * pi / 180])
qMax = np.array([ 30 * pi / 180, 5,  30 * pi / 180])

L = sum(a) # variable para representación gráfica
EPSILON = .01

#plt.ion() modo interactivo

# introducción del punto para la cinemática inversa
if len(sys.argv) != 3:
  sys.exit("El programa requiere dos parámetros: x y")
objetivo=[float(i) for i in sys.argv[1:]]

O=cin_dir(th,a)
 # Calculamos la posicion inicial
print ("- Posicion inicial:")
muestra_origenes(O)

dist = float("inf")
prev = 0.
iteracion = 1
while (dist > EPSILON and abs(prev-dist) > EPSILON/100.):
  prev = dist
  O=[cin_dir(th,a)]
  # Para cada combinación de articulaciones:
  for i in range(len(th)-1,-1,-1):
    if not prismatica[i]:
      # cálculo de la cinemática inversa:
      #print("Objetivo: ", objetivo, " pos: ", O[][i])

      v_obj = np.subtract(objetivo,O[-1][i])
      v_end = np.subtract(O[-1][-1],O[-1][i])
      angulo_desplazado = atan2((v_obj[0]*v_end[1]) - (v_obj[1]*v_end[0]),
                                (v_obj[0]*v_end[0]) + (v_obj[1]*v_end[1]))

      th[i] = th[i] - angulo_desplazado
      
      th[i] = (th[i] + pi) % (2 * pi) - pi
      th[i] = np.clip(th[i], qMin, qMax)[i]

    else:
      #Calcular omega = sumatorio de thetas anteriores
      omega = np.sum(th[:i+1])

      #Calcular d = unitario · (R - Ox)
      v_obj = np.subtract(objetivo,O[-1][i])
      v_end = np.subtract(O[-1][-1],O[-1][i])
      d = (cos(omega), sin(omega)) * np.subtract(v_obj, v_end)
      #L = L + d
      a[i] = a[i] + d
      #Limitar L según qMin y qMax
      a[i] = np.clip(a[i], qMin, qMax)[i]

    O.append(cin_dir(th,a))

  dist = np.linalg.norm(np.subtract(objetivo,O[-1][-1]))
  print ("\n- Iteracion " + str(iteracion) + ':')
  muestra_origenes(O[-1])
  muestra_robot(O,objetivo)
  print ("Distancia al objetivo = " + str(round(dist,5)))
  iteracion+=1
  O[0]=O[-1]

if dist <= EPSILON:
  print ("\n" + str(iteracion) + " iteraciones para converger.")
else:
  print ("\nNo hay convergencia tras " + str(iteracion) + " iteraciones.")
print ("- Umbral de convergencia epsilon: " + str(EPSILON))
print ("- Distancia al objetivo:          " + str(round(dist,5)))
print ("- Valores finales de las articulaciones:")
for i in range(len(th)):
  print ("  theta" + str(i+1) + " = " + str(round(th[i],3)))
for i in range(len(th)):
  print ("  L" + str(i+1) + "     = " + str(round(a[i],3)))
