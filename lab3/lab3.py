from sympy import *
from sympy.utilities.lambdify import lambdify
import vectormath as vmath
import math
import numpy as np
import matplotlib.pyplot as plt
from pandas import DataFrame

x, y, z, r = symbols('x y z r');

f = - x * y * z;

g1 = 2*x*y + 2*x*z + 2*y*z - 1;

h1 = - x;
h2 = - y;
h3 = - z;

B_with_r = f + 1/r * (g1**2 + Max(0, h1)**2 + Max(0, h2)**2 +  Max(0, h3)**2)

def calc_function(func, point):
    return N(func.subs(
        {
            Symbol('x'):point.x,
            Symbol('y'): point.y,
            Symbol('z'): point.z
        }));

def sub_r(B, r_replacement):
    return B.subs(Symbol('r'), r_replacement);

def calc_Xnew(Xh, Xc, teta):
    return vmath.Vector3(Xh.x + (1 + teta) * (Xc.x- Xh.x), Xh.y + (1 + teta) * (Xc.y- Xh.y), Xh.z + (1 + teta) * (Xc.z - Xh.z));

def check_e(xl, xg1, xg, xh, e):
    if((xl - xg).length < e):
        if((xl - xh).length < e):
            if((xl - xg1).length < e):
                if((xg - xh).length < e):
                    if((xg - xg1).length < e):
                        if((xg1 - xh).length < e):
                            return true
    return False

def deformuojamo_simplekso_algoritmas(func, X0, e = 0.0001, alpha = 0.5, teta = 1, beta = 0.5, gamma = 2, eta = -0.5):   # skaiciuojame delta1 ir delta2 pagal formule
    iterations = 0
    function_calls = 0
    
    delta1 = alpha * (math.sqrt(4) + 3) / (3 * math.sqrt(2))
    delta2 = alpha * (math.sqrt(4) - 1) / (3 * math.sqrt(2))
    
    X1 = vmath.Vector3(X0.x + delta2, X0.y + delta1, X0.z + delta1)
    X2 = vmath.Vector3(X0.x + delta1, X0.y + delta2, X0.z + delta1)
    X3 = vmath.Vector3(X0.x + delta1, X0.y + delta1, X0.z + delta2)

    X = [X0, X1, X2, X3]

    X.sort(key = lambda x: calc_function(func, x))

    Xl = X[0]
    Xg1 = X[1]
    Xg = X[2]
    Xh = X[3]

    fXl = calc_function(func, Xl)
    # fXg1 = calc_function(func, Xg1)
    fXg = calc_function(func, Xg)
    fXh = calc_function(func, Xh)

    function_calls += 4
    
    while(true):
        iterations += 1
        Xc = vmath.Vector3(
            (Xl.x + Xg.x + Xg1.x) / 3,
            (Xl.y + Xg.y + Xg1.y) / 3,
            (Xl.z + Xg.z + Xg1.z) / 3
            )

        Xnew = calc_Xnew(Xh, Xc, teta)
        fXnew = calc_function(func, Xnew)

        function_calls += 1
        
        if (Xnew.x <= 0 or Xnew.y <= 0 or Xnew.z <= 0):
            Xnew = calc_Xnew(Xh, Xc, eta)
            fXnew = calc_function(func, Xnew)

            function_calls += 1

        elif(fXnew < fXl):
            Z = calc_Xnew(Xh, Xc, gamma)
            fZ = calc_function(func, Z)

            function_calls += 1
            if(fZ < fXnew):
                Xnew = Z
                fXnew = fZ

        elif(fXnew > fXh): #simpleksas spaudziamas, nuo -1 iki 0
            Xnew = calc_Xnew(Xh, Xc, eta)
            fXnew = calc_function(func, Xnew)

            function_calls += 1
           
        elif(fXg < fXnew and fXnew < fXh): #simpleksas suspaudziamas, nuo 0 iki 1
            Xnew = calc_Xnew(Xh, Xc, beta)
            fXnew = calc_function(func, Xnew)

            function_calls += 1

        X=[Xl, Xg1, Xg, Xnew] 
        X.sort(key = lambda x: calc_function(func, x))

        Xl = X[0]
        Xg1 = X[1]
        Xg = X[2]
        Xh = X[3]

        fXl = calc_function(func, Xl)
        # fXg1 = calc_function(func, Xg1)
        fXg = calc_function(func, Xg)
        fXh = calc_function(func, Xh)

        if (check_e(Xl, Xg1, Xg, Xh, e)):
            break;
    return Xnew, function_calls, iterations

def baudos_metodas(B_with_r, x0, r = 10, C = 10, e = 0.0001):
    iterations = 0

    iterationList = []
    rList = []
    dsaIterationList = []
    dsaFunctionCallList = []
    xNewList = []

    B = sub_r(B_with_r, r)
    
    iterations += 1

    xnew, dsaFunctionCalls, dsaIterations = deformuojamo_simplekso_algoritmas(B, x0)

    iterationList.append(iterations)
    rList.append(r)
    dsaIterationList.append(dsaIterations)
    dsaFunctionCallList.append(dsaFunctionCalls)
    xNewList.append(xnew)

    while(true):
        iterations += 1
        r = r/C
        x0 = xnew
        B = sub_r(B_with_r, r)
        
        xnew, dsaFunctionCalls, dsaIterations = deformuojamo_simplekso_algoritmas(B, x0)
        
        iterationList.append(iterations)
        rList.append(r)
        dsaIterationList.append(dsaIterations)
        dsaFunctionCallList.append(dsaFunctionCalls)
        xNewList.append(xnew)

        if((x0 - xnew).length <e):
            break;
            
            
    df = DataFrame({'Iteration': iterationList, 'r': rList, 'dsaIterations': dsaIterationList, 'dsaFunctionCalls': dsaFunctionCallList, 'XList': xNewList})
    df.to_excel('output1.xlsx', sheet_name='sheet1')
    return xnew

def compute_min(X):
	print("Pradzios taskas: ", X);
	print("Rastas minimumas: ", baudos_metodas(B_with_r, X, 2, 10));

lamf = lambdify([x,y,z], f)
lamg1 = lambdify([x,y,z], g1)
lamh1 = lambdify([x,y,z], h1)
lamh2 = lambdify([x,y,z], h2)
lamh3 = lambdify([x,y,z], h3)
lamB = lambdify([x,y,z,r], B_with_r)

def calc(x):
  print(lamg1(x.x, x.y, x.z))
  print(lamh1(x.x, x.y, x.z))
  print(lamh2(x.x, x.y, x.z))
  print(lamh3(x.x, x.y, x.z))
  print(lamf(x.x, x.y, x.z))

def calc_r_cignificants(x_list, range):

  rList = []
  x0List = []
  x1List = []
  xmList = []

  for i in range:
    rList.append(i)
    x0List.append(lamB(x_list[0].x, x_list[0].y, x_list[0].z, i))
    x1List.append(lamB(x_list[1].x, x_list[1].y, x_list[1].z, i))
    xmList.append(lamB(x_list[2].x, x_list[2].y, x_list[2].z, i))

  df = DataFrame({'r': rList, 'x0': x0List, 'x1': x1List, 'xm': xmList})
  df.to_excel('output1.xlsx', sheet_name='sheet1')

def main():
  X0 = vmath.Vector3(0, 0, 0);
  X1 = vmath.Vector3(1, 1, 1);
  a,b,c = 9, 5, 6;
  Xm = vmath.Vector3(a/10, b/10, c/10);
  Xmin = vmath.Vector3(math.sqrt(1/6), math.sqrt(1/6), math.sqrt(1/6));

  compute_min(X0)

  # print(B_with_r)

  #calc_r_cignificants([X0, X1, Xm], range(1, 100))

  #calc(Xm)

if __name__ == "__main__":
    main()
