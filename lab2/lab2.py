from sympy import *
from sympy.utilities.lambdify import lambdify
import vectormath as vmath
import math
import numpy as np
import matplotlib.pyplot as plt
from pandas import DataFrame

x = Symbol('x')
y = Symbol('y')

f = - (x * y * (1 - y - x)) / 8

lam_f= lambdify([x,y],f)

f_diffX = lambdify([x,y], diff(f,x))
f_diffY = lambdify([x,y], diff(f,y))

print(diff(f,x))
print(diff(f,y))


def draw_gradient(min, max, levels, g_point_X_list, g_point_Y_list, title) :
  x = np.arange(min, max, 0.01)
  y = np.arange(min, max, 0.01)

  X, Y = np.meshgrid(x, y)
  Z = lam_f(X, Y)

  marker_size = 4

  #plt.contour(X, Y, Z, levels = levels)

  #plt.plot(g_point_X_list, g_point_Y_list, linewidth = 2)
  plt.plot(g_point_X_list, g_point_Y_list, 'bo', markersize = marker_size)

  annotation_start_point_index = 0
  #annotation_end_point_index = 
  annotation_end_point_index = len(g_point_Y_list)

  for i in range(annotation_start_point_index, annotation_end_point_index) :
    plt.annotate(str(i),
                (g_point_X_list[i], g_point_Y_list[i]),
                xycoords = 'data',
                xytext = (0, marker_size),
                textcoords = 'offset pixels',
                size = 10)

  if(annotation_end_point_index != len(g_point_X_list)) :
    plt.annotate(str(annotation_end_point_index) + " - " + str(len(g_point_X_list)-1),
                (g_point_X_list[annotation_end_point_index+1], g_point_Y_list[annotation_end_point_index+1]),
                xycoords = 'data',
                xytext = (0, marker_size),
                textcoords = 'offset pixels',
                size = 10)
  

  #plt.colorbar()
  plt.title(title)
  plt.xlabel('X')
  plt.ylabel('Y')
  plt.xlim([min, max])
  plt.ylim([min, max])

  plt.show()

def draw_simplex(min, max, levels, triangles_X, triangles_Y, annotate_point_X, annotate_point_Y, title) :
  x = np.arange(min, max, 0.01)
  y = np.arange(min, max, 0.01)

  X, Y = np.meshgrid(x, y)
  Z = lam_f(X, Y)

  #plt.contour(X, Y, Z, levels = levels)

  for i in range(len(triangles_X)) :
    plt.plot(triangles_X[i], triangles_Y[i], linewidth = 2)

  annotation_start_point_index = 0
  #annotation_end_point_index = 6
  annotation_end_point_index = len(annotate_point_X)

  for i in range(annotation_start_point_index, annotation_end_point_index) :
    plt.annotate(str(i),
                (annotate_point_X[i], annotate_point_Y[i]),
                xycoords = 'data',
                xytext = (0, 5),
                textcoords = 'offset pixels',
                size = 10)

  if(annotation_end_point_index != len(annotate_point_X)) :
    plt.annotate(str(annotation_end_point_index) + " - " + str(len(annotate_point_X)-1),
                (annotate_point_X[annotation_end_point_index+1], annotate_point_Y[annotation_end_point_index+1]),
                xycoords = 'data',
                xytext = (0, 5),
                textcoords = 'offset pixels',
                size = 10)

  #plt.colorbar()
  plt.title(title)
  plt.xlabel('X')
  plt.ylabel('Y')
  plt.xlim([min, max])
  plt.ylim([min, max])

  plt.show()

def nusiliedimo(x0, gama, e):
  iterations = 0
  function_calls = 0
  g_point_X_list = []
  g_point_Y_list = []

  X = x0
  g_point_X_list.append(X.x)
  g_point_Y_list.append(X.y)

  while(true):
    iterations += 1
    grad = vmath.Vector2(f_diffX(X.x, X.y), f_diffY(X.x, X.y))

    function_calls += 2

    X = vmath.Vector2(X.x - gama * grad.x, X.y - gama * grad.y)
    g_point_X_list.append(X.x)
    g_point_Y_list.append(X.y)
    print("Iteration: ", iterations, ", grad: ", grad, ", X: ", X)
    if(abs(grad.length) < e):
      break
    
  print("Iterations = ", iterations, "Function calls = ", function_calls, ", gama = ", gama, ", minPoint = ", X, ", result = ", lam_f(X.x, X.y))
  draw_gradient(0, 1.2, 100, g_point_X_list, g_point_Y_list, "Nusileidimo metodas")
  #draw_gradient(0.2, 0.5, 30, g_point_X_list, g_point_Y_list)

def auksinio_pjuvio_met(f, l, r, e):
  function_calls = 0

  ratio = (-1 + math.sqrt(5))/2

  L = r - l
  x1 = r - ratio * L
  x2 = l + ratio * L

  function = lambdify(x, f)

  fx1 = function (x1)
  fx2 = function (x2)
  function_calls += 2

  while (L >= e):
    if (fx2 < fx1):
      l = x1
      L = r - l
      x1 = x2
      x2 = l + ratio * L
      fx1 = fx2
      fx2 = function(x2)
    else:
      r = x2
      L = r - l
      x2 = x1
      x1 = r - ratio * L
      fx2 = fx1
      fx1 = function(x1)
  function_calls += 1
  if (fx1 < fx2):
    return x1, function_calls
  else:
    return x2, function_calls

def greiciausio_nusileidimo(x0, e):
  iterations = 0
  function_calls = 0
  g_point_X_list = []
  g_point_Y_list = []

  gamaList = []
  gradList = []
  XList = []

  X = x0
  gama = 0

  g_point_X_list.append(X.x)
  g_point_Y_list.append(X.y)

  while(true):
    iterations += 1
    grad = vmath.Vector2(f_diffX(X.x, X.y), f_diffY(X.x, X.y))

    function_calls += 2
    g = Symbol("g")

    Xx = X.x - g * f_diffX(X.x, X.y)
    Xy = X.y - g * f_diffY(X.x, X.y)

    newf = f.subs({Symbol('x'):Xx, Symbol('y'): Xy})
    newf = newf.subs({Symbol('x'):X.x, Symbol('y'): X.y})

    newf = newf.subs({Symbol('g'):Symbol('x')})

    gama, inner_f_calls = auksinio_pjuvio_met(newf, 0, 7.2, e)

    function_calls += inner_f_calls

    X = vmath.Vector2(X.x - gama * grad.x, X.y - gama * grad.y)

    gamaList.append(gama)
    gradList.append(grad)
    XList.append(X)

    g_point_X_list.append(X.x)  
    g_point_Y_list.append(X.y)

    if(abs(grad.length) < e):
      break
  df = DataFrame({'Gama': gamaList, 'Grad': gradList, 'X': XList}, index = range(1,iterations+1))
  df.to_excel('output1.xlsx', sheet_name='sheet1')
  print("Iterations = ", iterations, "Function calls = ", function_calls, ", gama = ", gama, ", minPoint = ", X, ", result = ", lam_f(X.x, X.y))
  draw_gradient(0, 1.2, 0, g_point_X_list, g_point_Y_list, "Greiciausio nusileidimo metodas")

def calculate_new_X(Xh, Xc, teta) :
  return vmath.Vector2((Xh.x + (1+teta)*(Xc.x- Xh.x), Xh.y + (1+teta)*(Xc.y- Xh.y)))

def deformuojamo_simplekso(x0, e, alpha = 0.5, beta = 0.5, gamma = 2, teta = 1, eta = -0.5) :
  iterations = 0
  function_calls = 0
  triangles_X = []
  triangles_Y = []
  annotate_point_X = []
  annotate_point_Y = []

  delta1 = alpha * (math.sqrt(3) + 1) / ( 2 * math.sqrt(2))
  delta2 = alpha * (math.sqrt(3) - 1) / ( 2 * math.sqrt(2))

  x1 = vmath.Vector2(x0.x + delta2, x0.y + delta1)
  x2 = vmath.Vector2(x0.x + delta1, x0.y + delta2)

  annotate_point_X.append(x0.x)
  annotate_point_Y.append(x0.y)

  X = [x0, x1, x2]
  X.sort(key = lambda x: lam_f(x.x, x.y))

  Xl = X[0]
  Xg = X[1]
  Xh = X[2]

  fXl = lam_f(Xl.x, Xl.y)
  fXg = lam_f(Xg.x, Xg.y)
  fXh = lam_f(Xh.x, Xh.y)

  function_calls += 3

  triangles_X.append([Xl.x, Xg.x, Xh.x, Xl.x])
  triangles_Y.append([Xl.y, Xg.y, Xh.y, Xl.y])

  while(true): 
    iterations += 1
    Xc = vmath.Vector2( (Xl.x + Xg.x) / 2, (Xl.y + Xg.y) / 2)
    Xnew = calculate_new_X(Xh, Xc, teta)
    fXnew = lam_f(Xnew.x, Xnew.y)

    function_calls += 1

    if (Xnew.x <= 0 or Xnew.y <= 0) :
      Xnew = calculate_new_X(Xh, Xc, eta)
      fXnew = lam_f(Xnew.x, Xnew.y)

      function_calls += 1

    elif (fXnew < fXl) : #simplekas pleciamas
      Z = calculate_new_X(Xh, Xc, gamma)
      yz = lam_f(Z.x, Z.y)

      function_calls += 1

      if(yz < fXnew):
        Xnew = Z
        fXnew = yz 
    elif (fXnew > fXh) : #simpleksas spaudziamas, nuo -1 iki 0
      Xnew = calculate_new_X(Xh, Xc, eta)
      fXnew = lam_f(Xnew.x, Xnew.y)

      function_calls += 1

    elif (fXg < fXnew and fXnew < fXh) : #simpleksas suspaudziamas, nuo 0 iki 1
      Xnew = calculate_new_X(Xh, Xc, beta)
      fXnew = lam_f(Xnew.x, Xnew.y)

      function_calls += 1

    X = [Xl, Xg, Xnew]

    annotate_point_X.append(Xnew.x)
    annotate_point_Y.append(Xnew.y)

    X.sort(key = lambda x: lam_f(x.x, x.y))

    Xl = X[0]
    Xg = X[1]
    Xh = X[2]

    triangles_X.append([Xl.x, Xg.x, Xh.x, Xl.x])
    triangles_Y.append([Xl.y, Xg.y, Xh.y, Xl.y])

    fXl = lam_f(Xl.x, Xl.y)
    fXg = lam_f(Xg.x, Xg.y)
    fXh = lam_f(Xh.x, Xh.y)

    if (abs((Xl - Xg).length) < e and abs((Xl - Xh).length) < e and abs((Xh - Xg).length) < e):
      break

  print("Iterations = ", iterations, "Function calls = ", function_calls, ", minPoint = ", Xnew, ", result = ", lam_f(Xnew.x, Xnew.y))
  draw_simplex(0, 1.5, 0, triangles_X, triangles_Y, annotate_point_X, annotate_point_Y, "Deformuojamo simplekso metodas")



def main():
  x0 = vmath.Vector2(0.0, 0.0)
  x1 = vmath.Vector2(1.0, 1.0)
  xM = vmath.Vector2(5/10,6/10)

  #print(x0, " ", lam_f(x0.x, x0.y), " (", f_diffX(x0.x, x0.y), ", ", f_diffY(x0.x, x0.y), ")")
  #print(x1, " ", lam_f(x1.x, x1.y), " (", f_diffX(x1.x, x1.y), ", ", f_diffY(x1.x, x1.y), ")")
  #print(xM, " ", lam_f(xM.x, xM.y), " (", f_diffX(xM.x, xM.y), ", ", f_diffY(xM.x, xM.y), ")")

  #for i in np.arange(0.1, 10.0, 0.2):
  #  min = nusiliedimo(xM, i, 0.0001)

  #nusiliedimo(xM, 9.3, 0.0001)
  #greiciausio_nusileidimo(xM, 0.0001)
  #greiciausio_nusileidimo(xM, 0.0001)
  deformuojamo_simplekso(xM, 0.0001)
  #deformuojamo_simplekso(x1, 0.0001)
  #deformuojamo_simplekso(xM, 0.0001)

if __name__ == "__main__":
    main()
