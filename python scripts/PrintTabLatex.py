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
