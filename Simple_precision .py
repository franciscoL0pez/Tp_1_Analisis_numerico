import math
import numpy as np 
from sigfig import round
import warnings 
𝜋 = math.pi
warnings.filterwarnings("ignore", category=UserWarning) #Ignoramos las advertencias al usar sigfigs

# Armamos el codigo con simple precision 
def simple_precision(numero:float):
    
    simple_pres = round(numero,sigfigs = 4) #Trabajamos con simple presicion, redondeandoa  8 digitos significativas
        
    return simple_pres


def crear_matriz ( tam_matriz:int )->None :
    k = 1
    m = 1

    matriz = []
    for i in range(tam_matriz):
        k = k + i
        matriz.append([])

        for j in range(tam_matriz):
            m = m + j  
            a = simple_precision(1 / (k + m - 1))

            matriz[i].append(a) 
            m = 1
            
        k = 1 
    
    return matriz


def mostrar_matriz( matriz:list ) -> None :
    for fila in matriz:
        for valor in fila:

            cantidad_cifras = len(str(valor))
            valor_espaciado = str(valor)

            for i in range(22-cantidad_cifras):
                valor_espaciado += " "
            print( valor_espaciado, end="|")

        print()


def calcular_x_simple_pres( n:int  ) -> list:
    lista_de_incognitas = []
    for i in range(n):
        i = i + 1  
        a =  simple_precision((2*(𝜋) / (n - 1 ) ))

        x = simple_precision(math.cos(a) * (i - 1))

        lista_de_incognitas.append(x)
    
   
    return lista_de_incognitas 


def hallar_b_con_simple_press(matriz:list , vector_x:list) -> list :
    vector_b = []

    for i in range(len(vector_x)):
        lista = []
        for j in range(len(vector_x)):

            a = (matriz[i][j]) * (vector_x[j])
            
            lista.append( simple_precision(a) )
        
        c = sum(lista)
    
        vector_b.append([])
        vector_b[i].append(c)
    
    return np.array(vector_b)



def crear_matriz_aumentada(matriz:list, vector_b: list ) -> list :
    matriz_aumentada = np.concatenate((matriz,vector_b),axis=1)
    

    for i in range(len(vector_b)):
        matriz_aumentada[i][-1] = simple_precision(matriz_aumentada[i][-1])    #Ponemos en simple press los ultimos vecs de la mat ab
    
    return matriz_aumentada


def aplicar_s_p_a_filas(fila:list,cant_de_columnas:int)-> list: 
    for i in range(cant_de_columnas):
        fila[i] = simple_precision(fila[i])

    return fila 


def aplicar_s_p_a_columnas(columna:list,cant_de_filas:int)-> list:
    for i in range(cant_de_filas):
        columna[i][0] = simple_precision(columna[i][0])
    return columna


def aplicar_gauss_simple_press(mat_ab:list) -> list :
    tam_de_mat_ab = np.shape(mat_ab)
    n = tam_de_mat_ab[0]  # filas
    m = tam_de_mat_ab[1]  # columnas 
    
    lista_multiplicadores = []
    
    for i in range(0,n-1,1):
        pivote   = mat_ab[i,i]
        adelante = i + 1

        for k in range(adelante,n,1):
            multiplicador  = simple_precision(mat_ab[k,i]/pivote) 

            x = aplicar_s_p_a_filas(mat_ab[k,:], m) 
            
            y = aplicar_s_p_a_filas(mat_ab[i,:]*multiplicador, m)
           


            mat_ab[k,:] =  x - y
            mat_ab[k,:] = aplicar_s_p_a_filas(mat_ab[k,:], m  ) 


            lista_multiplicadores.append(multiplicador)
        
    
    return lista_multiplicadores


def ceros_debajo_diagonal(u):
    for i in range(len(u)):
        for j in range(len(u)):
            if i>j:
                u[i][j]=0


def armar_L(lista_multiplicadores:list, dimension:int):
 
    matriz_i =  np.identity(dimension)
    for j in range (dimension): 

        for i in range (dimension):

            if i >= 1 and matriz_i[i-1][j] != 0  :
                matriz_i[i][j] = lista_multiplicadores[0]

                lista_multiplicadores.pop(0)
   
    return matriz_i


def sustitucion_hacia_atras_simple_press(matriz:list, vector_b_escalonado:list):
    n = len(matriz)
    
    x = [0] * n
    
    for i in range(n-1, -1, -1):
        
        suma = 0
        for j in range(i+1, n):
            suma += matriz[i][j] * x[j]
            suma[0] = simple_precision(suma[0])
            

        x[i] = (vector_b_escalonado[i] - suma) / matriz[i][i]
       
        x[i][0] = simple_precision(x[i][0])
        
   
    return x


def sustitucion_hacia_adelante(matriz_L:list, vector_b):
    n = len(vector_b)
    x = [0] * n

    for i in range(n):
        suma = 0

        for j in range(i):

            suma += matriz_L[i][j] * x[j]
            suma[0] = simple_precision(suma[0])


        x[i] = (vector_b[i] - suma) /  matriz_L[i][i]
        x[i][0] = simple_precision(x[i][0])

    
    return x



def calcular_residuo_relativo(vector_b, matriz, soluciones):
    resultado = np.dot(matriz, soluciones)
    r = np.linalg.norm(vector_b-resultado)/np.linalg.norm(vector_b)
    return r


def numero_de_condicion(𝛿x,x):
    n =np.linalg.norm(𝛿x)/np.linalg.norm(x)

    return n*(10**8)


def main() -> None:
    n = 12

    matriz = crear_matriz(n)

    soluciones = calcular_x_simple_pres(n)
   
    vector_b = hallar_b_con_simple_press(matriz,soluciones)
  
    vector_b = aplicar_s_p_a_columnas(vector_b,n)
    
    mat_ab = crear_matriz_aumentada(matriz, vector_b)
  
    lista_multiplicadores = aplicar_gauss_simple_press(mat_ab)
   
    U, y = np.hsplit(mat_ab, [-1])
    
    ceros_debajo_diagonal(U)
   
    mat_L = armar_L(lista_multiplicadores, n )
    
    vector_b_escalonado = sustitucion_hacia_adelante(mat_L,vector_b)

    vector_x = sustitucion_hacia_atras_simple_press(U, y)
  
    calcular_residuo_relativo(vector_b, matriz, vector_x)
 
    r_absoluto =(vector_b - np.dot(matriz, vector_x))
 
    
    𝛿y = sustitucion_hacia_adelante(mat_L, r_absoluto)
 
    𝛿x = sustitucion_hacia_atras_simple_press(U, 𝛿y)
  
    

    


main()