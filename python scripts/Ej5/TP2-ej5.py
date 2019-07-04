"""
TP2 ANALISIS NUMERICO
Rocío Gallo
Facundo Monpelat
Metodo de Runge Kutta orden 4
Ejercicio 3
"""

import numpy as np #Manejo de arrays
import sys
sys.setrecursionlimit(1000)
import math 
# para graficar
import plotly
import plotly.plotly as py
plotly.tools.set_credentials_file(username='fmonpelat', api_key='YpD7z4O34340q0GZGbC7')
#plotly.tools.set_credentials_file(username='fmonpelat', api_key='xxxxxx')
import plotly.graph_objs as go
from pprint import pprint


T0      = 293.15 # temperatura inicial 
T1      = 1002 # temperatura T1 del material
T2      = 922.5 # temperatura T2 del material
rho     = 7850 # rho del material
OD      = 0.2448 # metros 
WT      = 0.01384 # metros
Lt      = 12 # metros
Hc      = 20 # Constante de transferencia del calor
C       = 480 # Constante calorimetrica del material
sigma   = 5.6703e-8 #W/(m**2 * K**4)
eps     = 0.85 #
L       = 50 #m
n_bol   = 50 #unidades o pasos
cad     = (np.round(-10 / 10000 * (97490-90000) + 35, 0))
cad     = 28
v_0     = L / (n_bol * cad) #m/s
S       = math.pi*OD*Lt # Superficie del material
m       = rho*math.pi*OD*WT*(1-WT/OD)*Lt # Masa del material

print('\nDatos calculados: \n')
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
    return (-1 / (m*C) * ( Hc * S * (T - T_inf(x)) + sigma * eps * S * (np.power(T,4) - np.power(T_inf(T),4)) ) )

def analiticS(t):
    return (T0-T1)*math.exp(-((Hc*S)/(m*C))*t)+T1

def skParam(T_1,T_2):
    #T1      = 983.15 # temperatura T1 del material
    #T2      = 983.15 # temperatura T2 del material
    global T1
    T1=T_1
    global T2
    T2 = T_2
    minutes_conversion = 60
    kelvin_conversion = 273.15
    x0 = 0. #tiempo inicial
    xf = n_bol * cad #tiempo final
    h = cad
    i = 0
    data_rk = []
    dictSoak = {"Tsk":0,"Sk":0}

    x, y = rungeKutta(f,analiticS,h,x0,xf,T0,data_rk,i)
    dictSoak = paramSal(data_rk)
    #Impresiones para pruebas
    #print('tiempo Sk: {0:.3f} [min] tiempo obj: {1:.3f} [min] temp mean: {2:.3f} [C] temp obj: {3:.2f} [C]'.format(
    #    dictSoak['Sk']/minutes_conversion,Sk_obj,dictSoak['Tsk']-kelvin_conversion,T_Sk_obj)) 
    return dictSoak




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

    aux_dict={
        'i':i,
        'X':x0,
        'Y':y0,
        'E':delta
    }

    if debug: pprint(aux_dict)

    data.append(aux_dict)

    if( round(x0) == xf ):
        if debug: print('X final: ' + str(data[-1]['X']) + '/ Y final: '+str(data[-1]['Y']) ) 
        return data[-1]['X'],data[-1]['Y']

    else:
        return rungeKutta( f, g, h, x0+h, xf, y, data, i+1)

#----------------------------------------------------------
# FUNCION GraficarSoaking(data)
#
# PARAMETROS
# data:  Datos generados de la funcion rungeKutta()
# USO    Imprime con Plotly un gráfico con 2 datos
#-----------------------------------------------------------
def GraficarSoaking(data,graph_title):
    # data1,data2 son arrays de diccionarios
    datax1 = []
    datay1 = []


    # transformamos de segundos a minutos para graficar en X y en el eje Y 
    # pasamos de kelvin a grados centigrados
    # en X dividimos por 60 para minutos
    # en Y restamos 273 para pasarlo a centigrados 
    kelvin_conversion=273
    minutes_conversion=60
    
    for i in range(len(data)):
        datax1.append(data[i]['X']/minutes_conversion)
        datay1.append(data[i]['Y']-kelvin_conversion)
        
    # buscamos cuando finaliza la T1 para la linea horizontal
    i = 0
    x = data[0]['X'] * v_0
    while T_inf(x) == T1:
        x = data[i]['X'] * v_0
        i = i + 1
    T1end = data[i]['X']/minutes_conversion
    
    trace0 = go.Scatter(
                        x=datax1,
                        y=datay1,
                        name = 'Datos de método Runge-Kutta',
                        mode = 'lines+markers',
    )
    
    
    layout = dict(
                annotations=[
                    dict(
                        x=abs((datax1[0]-T1end))/2,
                        y=(T1-kelvin_conversion)*0.95,
                        xref='x',
                        yref='y',
                        text='T1 = '+str(T1-kelvin_conversion)+'[C]',
                        showarrow=False,
                        arrowhead=7,
                        ax=0,
                        ay=20
                    ),
                    dict(
                        x=T1end+ (datax1[-1]-T1end)/2,
                        y=(T2-kelvin_conversion)*0.95,
                        xref='x',
                        yref='y',
                        text='T2 = '+str(T2-kelvin_conversion)+'[C]',
                        showarrow=False,
                        arrowhead=7,
                        ax=0,
                        ay=20
                    )
                ],
                shapes= [
                           {
                                'type': 'line',
                                'x0': datax1[0],
                                'y0': T1-kelvin_conversion,
                                'x1': T1end,
                                'y1': T1-kelvin_conversion,
                                'line': {
                                    'color': 'rgb(233, 69, 34)',
                                    'width': 1.5,
                                    'dash': 'solid',
                                },
                            },
                            {
                                'type': 'line',
                                'x0': T1end,
                                'y0': T2-kelvin_conversion,
                                'x1': datax1[-1],
                                'y1': T2-kelvin_conversion,
                                'line': {
                                    'color': 'rgb(25, 145, 214)',
                                    'width': 1.5,
                                    'dash': 'solid',
                                },
                            },
                ],
                title = graph_title,
                xaxis = dict(title = 'X (Minutos) '),
                yaxis = dict(title = 'Y (Grados Centigrados)',
                             autorange=True),
                )
    data = [trace0]
    #la figura es un array de data(traces) y el layout (que tambien es un array)
    fig = dict(data=data, layout=layout)
    plotly.offline.plot(fig, auto_open=True)
    return

#----------------------------------------------------------
# FUNCION print_data(datos)
#
# PARAMETROS
# printdata: Lista de diccionarios con los valores X, Y i y el error de cada iteracion.
# USO       Imprime los valores en CSV para facil lectura desde la terminal.
#-----------------------------------------------------------
def print_data(data):
    kelvin_conversion=273
    minutes_conversion=60

    for i in range(len(data)):
        print ('iteracion={0:.2f},X={1:.5f},Y={2:.10f}'.format(i,data[i]['X']/minutes_conversion,data[i]['Y']-kelvin_conversion))
   

    return

#----------------------------------------------------------
# FUNCION paramSal(datos)
#
# PARAMETROS
# datos:  Lista de diccionarios con los datos de la función rungeKutta()
# USO     Imprime los datos de soaking y tiempos de soaking
#-----------------------------------------------------------
def paramSal(datos):

    debug = 0
    kelvin_conversion = 273
    minutes_conversion = 60
    SoakingTemp = 10 # temp de diferencia soaking en grados
    Tfin = datos[-1]['Y']
    tfin = datos[-1]['X']
    sumT = 0
    Tsk  = 0
    Sk   = 0
    nInc = 1

    #recorro T(t) de forma decreciente
    for i in range(len(datos)-1,0,-1):
        if debug: print('i= '+str(i)+' y: '+str(datos[i]['Y'])+' incrementoProm: '+str(nInc)+' tiempoIni: '+str(datos[i]['X']) )
        #Incrementos
        sumT = sumT + datos[i]['Y']
        nInc = nInc + 1

        #Temperatura y tiempo iniciales 
        Tini = datos[i]['Y']
        tini = datos[i]['X']

        if Tfin-Tini < SoakingTemp:
            continue
        else:
            Sk = tfin-tini
            break
    Tsk = sumT/nInc

    #print("Tsk(K):  "+ str(Tsk) +" K")
    #print("Tsk(C):  "+ str(Tsk-kelvin_conversion)+" C")
    #print("Sk(seg):  "+ str(Sk) +" seg")
    #print("Sk(min):  "+ str(Sk/minutes_conversion) +" min")

    return {"Tsk":Tsk-273.5,"Sk":Sk/60}

#---------------------------------------------------------------------------------------------------------------
# FUNCION
# SolPuntoFijo(F,J,nErr,Tsk0,Ssk0)
#
# PARAMETROS
# _F_     Diccionario de 2 elementos: F={"f1":f1,"f2":f2}
# _J_     Matriz de 2x2 Jacobiana inversa 
# _nErr_  Tolerancia entre una iteración y la siguiente
# _Tsk0_  Semilla de temperatura de soaking
# _Ssk0_  Semilla de tiempo de soaking
#
# RESULTADO
#_{"T1":T1 , "T2": T2,"nIter":nIter}_ Diccionario con los dos valores de las raíces y la cantidad de iteraciones
# Observaciones: El método deja prefijados aquellos valores que se conservan
#                en los casos A, B y C
#
#---------------------------------------------------------------------------------------------------------------

def SolPuntoFijo(T1_0,T2_0,Tsk0,Sk0=10,J=[[-0.5,1.5],[1.5,-0.5]],nErr=5):

    nIter=1
    #Semillas
    T1Previo=Tsk0
    T2Previo=Sk0

    aPar = skParam(T1_0,T2_0)
    F = {"f1":aPar["Tsk"]-Tsk0, "f2":aPar["Sk"]-Sk0}

    #Primera iteracion

    T1 = Tsk0-(J[0][0])*F["f1"]-(J[0][1])*F["f2"]
    T2 = Sk0-J[1][0]*F["f1"]-J[1][1]*F["f2"]


    while abs(T1 - T1Previo) > nErr and abs(T2 - T2Previo) > nErr :
    #while abs(T2 - T2Previo) > nErr :

        T1Previo=T1
        T2Previo=T2
        
        # T1Previo y T2Previo estan en celcius pero skParam recibe en kelvin
        aPar = skParam(T1Previo+273.15,T2Previo+273.15)
        F = {"f1":aPar["Tsk"]-Tsk0, "f2":aPar["Sk"]-Sk0}

        #Iteración de punto fijo
        #La multiplicacion de Matriz*Vector fue armada de forma analítica
        T1 = T1Previo-(J[0][0])*F["f1"]-(J[0][1])*F["f2"]
        T2 = T2Previo-J[1][0]*F["f1"]-J[1][1]*F["f2"]

        nIter = nIter+1
        if nIter == 5000:
            print("No hubo convergencia")
            break

    return {"T1":T1 , "T2": T2 ,"nIter":nIter}

    
#----------------------------------------------------------
# FUNCION PrintTabLatex()
#
# PARAMETROS
# aTitulos: Lista con los títulos a imprimir
# aDatos:   Lista de listas con los datos correspondientes 
#           a cada titulo
# USO       Imprime una tabla en formato latex para facilidad 
#           del pasaje de datos al informe
#-----------------------------------------------------------


def PrintTabLatex(aTitulos,aDatos):

    cCPypeC = "|c|"
    aTab    = [] 
    cLin    = ""
    cTit    = ""
    cIni    = "\hline "
    cFin    = " \\\\"

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
    print("     \makegapedcells")
    print("     \centering")
    print("\\resizebox{0.7\\textwidth}{!}{")
    print("         \\begin{tabular}{"+cCPypeC+"}")
    print(cTit)
    for ni in range(len(aTab)):
        print(aTab[ni])
    print("         \hline")
    print("     \end{tabular}}")
    print("\end{table}")
    return 

def Graficar_horno(T_1,T_2,graph_title):
    global T1
    T1=T_1
    global T2
    T2 = T_2
    x0 = 0. #tiempo inicial
    xf = n_bol * cad #tiempo final
    h = cad
    i = 0
    data_rk = []

    x, y = rungeKutta(f,analiticS,h,x0,xf,T0,data_rk,i)
    GraficarSoaking(data_rk,graph_title)
    return

def main():
    #para probar buscamos runge kutta de orden 4 para la funcion f
    data_rk=[]
    data_an=[]
    aInvJ = [[-0.5,1.5],[1.5,-0.5]]
    
    x0 = 0. #tiempo inicial
    xf = n_bol * cad #tiempo final
    h = cad
    i=0
    T1_0 = 1003 # temperatura T1 del material
    T2_0 = 922.5#922.5 # temperatura T2 del material

    aTit = ["Caso","Sk","Tsk","T1","T2","Número"]
    aDat = []

    #print('Sk Objetivo [seg] : '+str(Sk_obj)+' K°')
    #print('Tsk Objetivo [K] = ' +str(T_Sk_obj)+' K°\n')

    Tsk_a=602
    print('\n\n\tCASO A ---> Sk: 10 min, Tsk: {0:} C°\n'.format(Tsk_a))
    aSol = SolPuntoFijo(T1_0,T2_0,Tsk_a)
    aDat.append(
        ["A","10",str(Tsk_a),str(round(aSol["T1"],3)),str(round(aSol["T2"],3)),str(aSol["nIter"])]
        )
    Graficar_horno(aSol["T1"],aSol["T2"],"Caso A")
    print('T1 '+str(aSol["T1"])+'C°')
    print('T2 '+str(aSol["T2"])+'C°\n')
    print(str(aSol["nIter"]))

    Tsk_b=624.90
    print('\n\tCASO B ---> Sk: 10 min, Tsk: {0:}°\n'.format(Tsk_b))
    aSol = SolPuntoFijo(T1_0,T2_0,Tsk_b)
    aDat.append(
        ["B","10",str(Tsk_b),str(round(aSol["T1"],3)),str(round(aSol["T2"],3)),str(aSol["nIter"])]
        )
    Graficar_horno(aSol["T1"],aSol["T2"],"Caso B")
    print('T1 '+str(aSol["T1"])+'C°')
    print('T2 '+str(aSol["T2"])+'C°\n')
    print(str(aSol["nIter"]))

    Tsk_c=674.90
    print('\n\tCASO C ---> Sk: 10 min, Tsk: {0:} C°\n'.format(Tsk_c))
    aSol = SolPuntoFijo(T1_0,T2_0,Tsk_c)
    aDat.append(
        ["C","10",str(Tsk_c),str(round(aSol["T1"],3)),str(round(aSol["T2"],3)),str(aSol["nIter"])]
        )
    print('T1 '+str(aSol["T1"])+'C°')
    print('T2 '+str(aSol["T2"])+'C°\n')
    Graficar_horno(aSol["T1"],aSol["T2"],"Caso C")
    print(str(aSol["nIter"]))

    #aDat = np.transpose(aDat)

    PrintTabLatex(aTit,aDat)

    return


if(__name__ == "__main__"):
    main()