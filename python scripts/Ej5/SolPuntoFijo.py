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
#---------------------------------------------------------------------------------------------------------------

def SolPuntoFijo(F,J,nErr,Tsk0,Ssk0):

	nIter=1
	#Semillas
	T1Previo=Tsk0
	T2Previo=Ssk0

	#Primera iteracion
	T1 = Tsk0-J[1][1]*F["f1"](Tsk0,Ssk0)-J[1][2]*F["f2"](Tsk0,Ssk0)
	T2 = Ssk0-J[2][1]*F["f1"](Tsk0,Ssk0)-J[2][2]*F["f2"](Tsk0,Ssk0)

	while (T1 - T1Previo) > nErr and (T2 - T2Previo) > nErr :

		T1Previo=T1
		T2Previo=T2

		#Iteración de punto fijo
		#La multiplicacion de Matriz*Vector fue armada de forma analítica
		T1 = T1Previo-J[1][1]*F["f1"](T1Previo,T2Previo)-J[1][2]*F["f2"](T1Previo,T2Previo)
		T2 = T2Previo-J[2][1]*F["f1"](T1Previo,T2Previo)-J[2][2]*F["f2"](T1Previo,T2Previo)
		nIter = nIter+1

	return {"T1":T1 , "T2": T2 ,"nIter":nIter}
