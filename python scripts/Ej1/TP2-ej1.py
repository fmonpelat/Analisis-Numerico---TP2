#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
TP2 ANALISIS NUMERICO
Rocío Gallo
Facundo Monpelat
Metodo de Runge Kutta orden 4
Ejercicio 1
"""

import numpy as np #Manejo de arrays
import sys
sys.setrecursionlimit(10000000)
import math 
# para graficar
import plotly
import plotly.plotly as py
plotly.tools.set_credentials_file(username='fmonpelat', api_key='YpD7z4O34340q0GZGbC7')
#plotly.tools.set_credentials_file(username='fmonpelat', api_key='xxxxxx')
import plotly.graph_objs as go
# debug
from pprint import pprint

T0 = 293.15 # temperatura inicial 
T1 = 923.15 # temperatura T1 del material
T2 = 923.15 # temperatura T2 del material
rho = 7850 # rho del material
OD = 0.2448 # metros 
WT = 0.01384 # metros
Lt = 12 # metros
Hc= 20 # Constante de transferencia del calor
C= 480 # Constante calorimetrica del material
sigma = 5.6703e-8 #W/(m**2 * K**4)
eps = 0.85 #
L = 50 #m
n_bol = 50 #unidades o pasos
cad = np.round(-10 / 10000 * (97490-90000) + 35, 0)
v_0 = L / (n_bol * cad) #m/s
S = math.pi*OD*Lt # Superficie del material
m =  rho*math.pi*OD*WT*(1-WT/OD)*Lt # Masa del material

print('Datos calculados: ')
print('m [kg] = ' +str(m))
print('S [m2] = ' +str(S))
print('Cadencia [seg] : '+str(cad))
print('v_0 [m/s] = ' +str(v_0))

def T_inf(x):
    if x <= L/2:
        return T1
    else:
        return T2

# definición de la EDO en forma dT/dt = f(t,T)
def f(t,T): 
    x = t * v_0 #aproximacion continua, se mueve de a pasitos en realidad
    return (-((Hc*S)/(m*C))*(T-T_inf(x)))
    
# Respuesta analitica de la EDO (para comparar) T1=T2
def analiticS(t):
    return (T0-T1)*math.exp(-((Hc*S)/(m*C))*t)+T1


def main():
    # Declaramos los arrays de datos
    data_rk=[]
    data_euler=[]
    data_analitic=[]

    x0 = 0. #tiempo inicial
    xf = n_bol * cad #tiempo final
    h = cad
    i=0
    
    x_rk, y_rk = rungeKutta(f,analiticS,h,x0,xf,T0,data_rk,i)
    x_euler, y_euler = euler(f,analiticS,h,x0,xf,T0,data_euler)

    # Para ver por terminal los datos
    # print('Metodo de Runge Kutta\n')
    # print_data(data_rk)
    # print("Metodo de Euler\n")
    # print_data(data_euler)
    
    cargaAnalitic(h,data_rk[-1]['i'],data_analitic)
    # Graficamos los datos (Punto A)
    Graficar(data_rk,data_euler,data_analitic)
    # Graficamos el error para ambos metodos (Punto B)
    GraficarError(data_rk,data_euler)



#----------------------------------------------------------
# FUNCION euler(f,g,h,x0,xf,y0,data):
#
# PARAMETROS
# f:       Datos generados de la funcion rungeKutta()
# g:       Datos generados de la funcion euler()
# h:       Incremento o cadencia
# x0 y y0: Valores iniciales para RK
# xf :     Valor final para RK
# data :   Array de diccionarios que contiene los datos X e Y de cada iteracion
# USO      Calcula los valores con el metodo de Euler
#-----------------------------------------------------------
def euler(f,g,h,x0,xf,y0,data):
    t = x0
    y = y0
    i=0
    while t <= xf:
        delta = abs( y-g(t) )/g(t)
        aux_dict={
        'i':i,
        'X':t,
        'Y':y,
        'E':delta,
        }
        data.append(aux_dict)
        i=i+1
        y += h * f(t,y)
        t += h
    return data[-1]['X'],data[-1]['Y']



#----------------------------------------------------------
# FUNCION rungeKutta(f,g,h,x0,xf,y0,data,i)
#
# PARAMETROS
# f:       Datos generados de la funcion rungeKutta()
# g:       Datos generados de la funcion euler()
# h:       Incremento o cadencia
# x0 y y0: Valores iniciales para RK
# xf :     Valor final para RK
# data :   Array de diccionarios que contiene los datos X e Y de cada iteracion
# i :      Contador de iteraciones 
# USO      Calcula los valores con el metodo de Runge-Kutta
#-----------------------------------------------------------
def rungeKutta(f,g,h,x0,xf,y0,data,i):
    debug = 0
    # K1 = F(Xn,Yn)
    # K2 = F( Xn + h/2 , Yn + h*K1 /2 )
    # K3 = F( Xn + h/2 , Yn + h*K2 /2 )
    # K4 = F( Xn + h/2 , Yn + h*K3 )
    k1 = f(x0, y0)
    k2 = f(x0 + (h/2), y0 + h*(k1/2))
    k3 = f(x0 + (h/2), y0 + h*(k2/2))
    k4 = f(x0 + h, y0 + h*k3)
    # Yn+1 = Yn + h/6( K1 + 2K2 + 2K3 + K4 )
    y = y0 + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    delta = abs( y0-g(x0) )/g(x0)
    if debug: print ('dif('+str(y)+'-'+str(g(x0))+' )= '+str(delta)+'\n')
    aux_dict={
        'i':i,
        'X':x0,
        'Y':y0,
        'E':delta,
    }
    if debug: pprint(aux_dict)
    data.append(aux_dict)
    if( round(x0) == xf ):
        if debug: print('X final: ' + str(data[-1]['X']) + ' -  Y final: '+str(data[-1]['Y']) ) 
        return data[-1]['X'],data[-1]['Y']
    else:
        return rungeKutta( f, g, h, x0+h, xf, y, data, i+1)


#----------------------------------------------------------
# FUNCION cargaAnalitic(h,nCant,data)
#
# PARAMETROS
# h:      Intervalo de incremento de X
# nCant:  Cantidad de iteraciones en incrementos de h
# data:   Array de diccionarios con datos de X, Y y las iteraciones realizadas
# USO     Inicializa el array de datos con la funcion analitica para luego graficar 
#         con incrementos de h nCant de veces
#-----------------------------------------------------------
def cargaAnalitic(h,nCant,data):
   for i in range(nCant+1):
       aux_dict={
           'i':i,
           'X':i*h,
           'Y':analiticS(i*h),
       }
       data.append(aux_dict)
   return data[-1]['X'],data[-1]['Y']

#----------------------------------------------------------
# FUNCION GraficarError(data_rk,data2_euler)
#
# PARAMETROS
# data_rk:      Datos generados de la funcion rungeKutta()
# data2_euler:  Datos generados de la funcion euler()
# USO           Imprime un gráfico con Plotly el error relativo porcentual de los datos
#               del metodo de runge kutta y del metodo de euler. 
#-----------------------------------------------------------
def GraficarError(data_rk,data2_euler):
    datax1 = []
    datay1 = []
    datax2 = []
    datay2 = []
    minutes_conversion=60

    for i in range(len(data_rk)):
      datax1.append(data_rk[i]['X']/minutes_conversion)
      datay1.append(data_rk[i]['E'])

    for i in range(len(data2_euler)):
      datax2.append(data2_euler[i]['X']/minutes_conversion)
      datay2.append(data2_euler[i]['E'])

    trace0 = go.Scatter(
                        x=datax1,
                        y=datay1,
                        name = 'Error relativo porcentual RK',
                        mode = 'lines+markers',
    )
    trace1 = go.Scatter(
                        x=datax2,
                        y=datay2,
                        name = 'Error relativo porcentual Euler',
                        mode = 'lines+markers',
    )
    layout = dict(title = 'Error porcentual por método',
                xaxis = dict(title = 'X (Minutos) '),
                yaxis = dict(title = 'Y (Error %)',
                             autorange=True),
                )
    data = [trace0,trace1]
    fig = dict(data=data, layout=layout)
    plotly.offline.plot(fig, auto_open=True)
    return
   

#----------------------------------------------------------
# FUNCION Graficar(data,data2)
#
# PARAMETROS
# data:  Datos generados de la funcion rungeKutta()
# data2: Datos generados de la funcion euler()
# data3: Datos de 
# USO    Imprime con Plotly un gráfico con 2 datos y uno que se 
#        hace con la solucion analitica.
#-----------------------------------------------------------
def Graficar(data,data2,data3):
    # data1,data2 son arrays de diccionarios
    datax1 = []
    datay1 = []
    datax2 = []
    datay2 = []
    datax3 = []
    datay3 = []

    # transformamos de segundos a minutos para graficar en X y en el eje Y 
    # pasamos de kelvin a grados centigrados
    # en X dividimos por 60 para minutos
    # en Y restamos 273 para pasarlo a centigrados 
    kelvin_conversion=273
    minutes_conversion=60

    for i in range(len(data)):
        datax1.append(data[i]['X']/minutes_conversion)
        datay1.append(data[i]['Y']-kelvin_conversion)
    for i in range(len(data2)):
        datax2.append(data2[i]['X']/minutes_conversion)
        datay2.append(data2[i]['Y']-kelvin_conversion)
    for i in range(len(data3)):
        datax3.append(data3[i]['X']/minutes_conversion)
        datay3.append(data3[i]['Y']-kelvin_conversion)

    trace0 = go.Scatter(
                        x=datax1,
                        y=datay1,
                        name = 'Datos de método Runge-Kutta',
                        mode = 'lines+markers',
    )
    trace1 = go.Scatter(
                        x=datax2,
                        y=datay2,
                        name = 'Datos de método euler',
                        mode = 'lines+markers',
    )
    trace2 = go.Scatter(
                        x=datax3,
                        y=datay3,
                        name = 'Datos de la función analitica',
                        mode = 'lines+markers',
    )
    
    layout = go.Layout(title = 'Gráfica comparativa',
                xaxis = dict(title = 'X (Minutos) '),
                yaxis = dict(title = 'Y (Grados Centigrados)',
                             autorange=True,
                             hoverformat = '.10f'),
                )
    data = [trace0,trace1,trace2]
    #la figura es un array de data(traces) y el layout (que tambien es un array)
    fig = dict(data=data, layout=layout)
    plotly.offline.plot(fig, auto_open=True)
    return


#----------------------------------------------------------
# FUNCION print_data(datos)
#
# PARAMETROS
# printdata: Lista de diccionarios con los valores X, Y i y el error de cada iteracion.
# USO  		Imprime los valores en CSV para facil lectura desde la terminal.
#-----------------------------------------------------------
def print_data(data):
    kelvin_conversion=273
    minutes_conversion=60

    for i in range(len(data)):
        exact = analiticS(data[i]['X'])
        print ('iteracion={0:.2f},X={1:.5f},Y={2:.10f},exact={3:.10f},error(%)={4:.10f}'.format(i,data[i]['X'],data[i]['Y']-kelvin_conversion,exact-kelvin_conversion,data[i]['E']*100))
    return

if(__name__ == "__main__"):
    main()