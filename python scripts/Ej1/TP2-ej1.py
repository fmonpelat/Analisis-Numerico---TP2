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

def edo(t,T):
    x = t * v_0 #aproximacion continua, se mueve de a pasitos en realidad
    return (-1 / (m*C) * ( Hc * S * (T - T_inf(x)) + sigma * eps * S * (T**4 - T_inf(x)**4) ) )

def f(t,T): 
    x = t * v_0 #aproximacion continua, se mueve de a pasitos en realidad
    return (-((Hc*S)/(m*C))*(T-T_inf(x)))
    
# Respuesta analitica de la EDO (para comparar) T1=T2
def analiticS(t):
    return (T0-T1)*math.exp(-((Hc*S)/(m*C))*t)+T1


def main():
    #para probar buscamos runge kutta de orden 4 para la funcion f
    data_rk=[]
    data_euler=[]

    x0 = 0. #tiempo inicial
    xf = n_bol * cad #tiempo final
    h = cad
    i=0

    #x, y = rungeKutta(f,analiticS,h,x0,xf,T0,data_rk,i)
    #print('Metodo de Runge Kutta\n')
    #print_data(data_rk)
    x_euler, y_euler = euler(f,analiticS,h,x0,xf,T0,data_euler)
    print("Metodo de Euler\n")
    print_data(data_euler)
    #Graficar(data_rk,data_euler)
    #GraficarError(data_rk,data_euler)



# Función euler(f,x0,y0,xf,h,data)
# f : funcion EDO => dT/dt = f
# x0 y y0 : valores iniciales para euler
# xf: valor final en x
# h: incremento de euler
# data: array de diccionario que contiene la tupla X e Y de cada iteracion
def euler(f,g,h,x0,xf,y0,data):
    t = x0
    y = y0
    while t <= xf:
        y0=y
        t += h
        y += h * f(t,y)
        delta = abs( y-g(t) )/g(t)
        aux_dict={
        'i':i,
        'X':t,
        'Y':y,
        'E':delta,
        }
        data.append(aux_dict)
    return data[-1]['X'],data[-1]['Y']



# Función rungeKutta(f,h,x0,y0,val,data,i)
# f : funcion EDO => dT/dt = f
# g : funcion analitica para el calculo del error
# x0 y y0: valores iniciales para RK
# xf : valore final para RK
# n_steps: cantidad de pasos al cual dividir la resta de xf y x0
# data : array de diccionarios que contiene la tupla X e Y de cada iteracion
# i : contador de iteraciones 

def rungeKutta(f,g,h,x0,xf,y0,data,i):
    debug = 1
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
    delta = abs( y-g(x0) )/g(x0)

    aux_dict={
        'i':i,
        'X':x0,
        'Y':y0,
        'E':delta,
    }
    if debug: pprint(aux_dict)

    if( round(x0) == xf ):
        if debug: print('X final: ' + str(data[-1]['X']) + '/ Y final: '+str(data[-1]['Y']) ) 
        return data[-1]['X'],data[-1]['Y']
    else:
        data.append(aux_dict)
        return rungeKutta( f, g, h, x0+h, xf, y, data, i)


def GraficarError(data_rk,data2_euler):
    datax1 = []
    datay1 = []
    datax2 = []
    datay2 = []
    kelvin_conversion=273
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
    layout = dict(title = 'Gráfica comparación del Error porcentual',
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
# data: 
# data2: 
# USO    Imprime con Plotly un gráfico con 2 datos y uno que se 
#        hace con la solucion analitica.
#-----------------------------------------------------------
def Graficar(data,data2):
    # data1,data2 son arrays de diccionarios
    datax1 = []
    datay1 = []
    datax2 = []
    datay2 = []

    # transformamos de segundos a minutos para graficar en X y en el eje Y 
    # pasamos de kelvin a grados centigrados
    # en X dividimos por 60 para minutos
    # en Y restamos 273 para pasarlo a centigrados 
    kelvin_conversion=273
    minutes_conversion=60
    xx = np.linspace(0, data[-1]['X'], 50)
    yy=[]
   
    for i in range(len(xx)):
        yy.append(g(xx[i]))
    for i in range(len(xx)):
        xx[i]=xx[i]/minutes_conversion
        yy[i]=yy[i]-kelvin_conversion
    for i in range(len(data)):
        datax1.append(data[i]['X']/minutes_conversion)
        datay1.append(data[i]['Y']-kelvin_conversion)
    for i in range(len(data2)):
        datax2.append(data2[i]['X']/minutes_conversion)
        datay2.append(data2[i]['Y']-kelvin_conversion)

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
                        x=xx,
                        y=yy,
                        name = 'Datos de gráfica analitica',
                        mode = 'lines+markers',
    )
    
    layout = dict(title = 'Gráfica comparación',
                xaxis = dict(title = 'X (Minutos) '),
                yaxis = dict(title = 'Y (Grados Centigrados)',
                             autorange=True),
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

    print ('iteracion={0:.2f},X={1:.5f},Y={2:.10f},exact={3:.10f},error(%)={4:.10f}'.format(i,data[i]['X']/minutes_conversion,data[i]['Y']-kelvin_conversion,exact-kelvin_conversion,data[i]['E']*100))
    return



#----------------------------------------------------------
# FUNCION PrintTabLatex()
#
# PARAMETROS
# aTitulos: Lista con los títulos a imprimir
# aDatos: 	Lista de listas con los datos correspondientes 
# 		  	a cada titulo
# USO  		Imprime una tabla en formato latex para facilidad 
# 			del pasaje de datos al informe
#-----------------------------------------------------------
def PrintTabLatex(aTitulos,aDatos):

	cCPypeC = "|"
	aTab	= [] 
	cLin	= ""
	cTit	= ""
	cIni	= "\hline "
	cFin	= " \\\\"

	for nj in range(len(aTitulos)):
		if nj < len(aTitulos)-1:
			cTit = cTit + aTitulos[nj]+ " & "
		if nj == len(aTitulos)-1:
			cTit = cTit + aTitulos[nj]

	cTit = cIni + cTit + cFin
	for ni in range(len(aDatos)):
		cCPypeC=cCPypeC+"c|"
		cLin = ""
		for nj in range(len(aDatos[ni])):
			if nj < len(aDatos[ni])-1:
				cLin = cLin + aDatos[ni][nj]+ " & "
			if nj == len(aDatos[ni])-1:
				cLin = cLin + aDatos[ni][nj]

		cLin = cIni + cLin + cFin 
		aTab.append(cLin)


	print("\\begin{table}[H]")
	print("		\makegapedcells")
	print("		\centering")
	print("\\resizebox{0.7\\textwidth}{!}{")
	print("			\\begin{tabular}{"+cCPypeC+"}")
	print(cTit)
	for ni in range(len(aTab)):
		print(aTab[ni])
	print("			\hline")
	print("		\end{tabular}}")
	print("\end{table}")
	return

if(__name__ == "__main__"):
    main()