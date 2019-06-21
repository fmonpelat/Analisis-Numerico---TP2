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

# EDO (no se si es correcta; la variable y seria T)
def f(x,y): 
    return -0.00040653*(y-650)

#EDO respuesta analitica (nose si esta bien)
def g(x):
    return 650*(20-650)*math.exp(-0.00040653*x)

def main():
    #para probar buscamos runge kutta de orden 4 para la funcion x+y
    printdata=[]
    x0=0
    y0=20
    # y(5)? con paso de h=0.5
    val=650
    h=28
    x, y = rungeKutta(f,h,x0,y0,val,printdata)
    print_data(printdata)
    GraficarDiferencia(printdata)



def rungeKutta(f,h,x0,y0,val,data):
    debug = 1
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

    aux_dict={
        'X':x,
        'Y':y,
    }

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
    yy=[]
    xx = np.linspace(0, 100, 100)
    for i in range(len(xx)):
        yy.append(g(xx[i]))

    for i in range(len(data)):
      datax.append(data[i]['X'])
      datay.append(data[i]['Y'])
    pprint(yy)

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
    
    layout = dict(title = 'Gráfica de Runge-Kutta',
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