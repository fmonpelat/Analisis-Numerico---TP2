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

# Runge Kutta orden 4
# Yn+1 = Yn + h/6( K1 + 2K2 + 2K3 + K4 )
# K1 = F(Xn,Yn)
# K2 = F( Xn + h/2 , Yn + h*K1 /2 )
# K3 = F( Xn + h/2 , Yn + h*K2 /2 )
# K4 = F( Xn + h/2 , Yn + h*K3 )
# Xn es la condicion inicial

## NOTA: debemos calcular la superficie Hc y demas ya que debe ser versatil! lo dice el TP

# EDO 
def f(t,y): 

    T1 = 923.15 # temperatura final del material
    rho = 7850 # rho del material
    OD = 0.2448 # metros 
    WT = 0.01384 # metros
    Lt = 12 # metros
    Hc= 20 # Constante de transferencia del calor
    C= 480 # Constante calorimetrica del material

    S = math.pi*OD*Lt # Superficie del material
    m =  rho*math.pi*OD*WT*(1-WT/OD)*Lt # Masa del material
    return -((Hc*S)/(m*C))*(t-T1)

#EDO respuesta analitica (para comparar)
def g(t):
    T0 = 293.15 # temperatura inicial del material
    T1 = 923.15 # temperatura final del material
    rho = 7850 # rho del material
    OD = 0.2448 # metros 
    WT = 0.01384 # metros
    Lt = 12 # metros
    Hc= 20 # Constante de transferencia del calor
    C= 480 # Constante calorimetrica del material

    S = math.pi*OD*Lt # Superficie del material
    m =  rho*math.pi*OD*WT*(1-WT/OD)*Lt # Masa del material
    return (T0-T1)*math.exp(-((Hc*S)/(m*C))*t)+T1

def main():
    #para probar buscamos runge kutta de orden 4 para la funcion f
    printdata=[]
    x0=0
    T0=293.15 # temperatura inicial 

    T1 = 923.15 # temperatura final
    h = 28 # cadencia

    x, y = rungeKutta(f,h,x0,T0,T1,printdata)
    print_data(printdata)
    GraficarDiferencia(printdata)



def rungeKutta(f,h,x0,y0,val,data):
    debug = 1
    aux_dict={
        'X':x0,
        'Y':y0,
    }
    # K1 = F(Xn,Yn)
    # K2 = F( Xn + h/2 , Yn + h*K1 /2 )
    # K3 = F( Xn + h/2 , Yn + h*K2 /2 )
    # K4 = F( Xn + h/2 , Yn + h*K3 )
    x = x0
    k1 = f(x0, y0)
    k2 = f(x0 + (h/2), y0 + h*(k1/2))
    k3 = f(x0 + (h/2), y0 + h*(k2/2))
    k4 = f(x0 + h, y0 + h*k3)
    # Yn+1 = Yn + h/6( K1 + 2K2 + 2K3 + K4 )
    y = y0 + (h/6)*(k1 + 2*k2 + 2*k3 + k4)


    if(x0 > val):
        if (debug==1): print('encontramos el valor!')
        return data[-1]['X'],data[-1]['Y']
    else:
        data.append(aux_dict)
        return rungeKutta(f, h, x0+h, y, val, data)



def GraficarDiferencia(data):
    # data1,data2,data3 son arrays con puntos a ser graficados
    # TODO hay que ajustar estos valores de linspace ya que el gráfico es logaritmico
    datax = []
    datay = []
    xx = np.linspace(0, 1000, 50)
    yy=[]
    for i in range(len(xx)):
        yy.append(g(xx[i]))

    for i in range(len(data)):
      datax.append(data[i]['X'])
      datay.append(data[i]['Y'])

    trace0 = go.Scatter(
                        x=datax,
                        y=datay,
                        name = 'Datos de Runge-Kutta',
                        mode = 'lines+markers',
    )
    trace1 = go.Scatter(
                        x=xx,
                        y=yy,
                        name = 'Datos de grafica analitica',
                        mode = 'lines+markers',
    )
    
    layout = dict(title = 'Gráfica de Metodo Runge-Kutta',
                xaxis = dict(title = 'X '),
                yaxis = dict(title = 'Y ',
                             autorange=True),
                )
    data = [trace0,trace1]
    #la figura es un dicccionario de data(traces) y el layout (que tambien es un dicccionario)
    fig = dict(data=data, layout=layout)
    plotly.offline.plot(fig, auto_open=True)
    return

def print_data(printdata):
    for i in range(len(printdata)):
      print ('iteracion={0:.2f},X={1:.5f},Y={2:.5f}'.format(i,printdata[i]['X'],printdata[i]['Y']))
    return

if(__name__ == "__main__"):
    main()