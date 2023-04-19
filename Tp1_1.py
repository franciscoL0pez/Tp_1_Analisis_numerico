import math
import numpy as np 
from sigfig import round

𝜋 = math.pi
def crear_matriz ( tam_matriz:int )->None :
    k = 1
    m = 1

    matriz = []
    for i in range(tam_matriz):
        k = k + i
        matriz.append([])

        for j in range(tam_matriz):
            m = m + j  
            a = 1 / (k + m - 1)

            matriz[i].append(a) 
            m = 1
            
        k = 1 

    return matriz


def calcular_x( n:int  ) -> list:
    lista_de_incognitas = []
    for i in range(n):
        i = i + 1  
        a =  (2*(𝜋) / (n - 1 ) )

        x = math.cos(a) * (i - 1)

        lista_de_incognitas.append(x)
    
    
    return lista_de_incognitas 

def mostrar_matriz( matriz:list ) -> None :
    for fila in matriz:
        for valor in fila:

            cantidad_cifras = len(str(valor))
            valor_espaciado = str(valor)

            for i in range(22-cantidad_cifras):
                valor_espaciado += " "
            print( valor_espaciado, end="|")

        print()

def hallar_b(matriz:list , vector_x:list) -> list :
    vector_b = []

    for i in range(len(vector_x)):
        lista = []
        for j in range(len(vector_x)):
            a = matriz[i][j] *  vector_x[j]
            lista.append(a)
        
        c = sum(lista)
    
        vector_b.append([])
        vector_b[i].append(c)
    
    return np.array(vector_b)


def crear_matriz_aumentada(matriz:list, vector_b: list ) -> list :
    matriz_aumentada = np.concatenate((matriz,vector_b),axis=1)
    #matriz_aumentada_copia = np.copy(matriz_aumentada)

    return matriz_aumentada


def aplicar_gauss(mat_ab:list) -> list :
    tam_de_mat_ab = np.shape(mat_ab)
    n = tam_de_mat_ab[0]  # filas
    m = tam_de_mat_ab[1]  # columnas 
    lista_multiplicadores = []
    
    for i in range(0,n-1,1):
        pivote   = mat_ab[i,i]
        adelante = i + 1

        for k in range(adelante,n,1):
            multiplicador  = mat_ab[k,i]/pivote
            mat_ab[k,:] = mat_ab[k,:] - mat_ab[i,:]*multiplicador

            lista_multiplicadores.append(multiplicador)
        

    return lista_multiplicadores


def sustitucion_hacia_atras(matriz:list, vector_b_escalonado:list):
    n = len(matriz)
   
    x = [0] * n

    for i in range(n-1, -1, -1):

        suma = 0
        for j in range(i+1, n):
            suma += matriz[i][j] * x[j]

        x[i] = (vector_b_escalonado[i] - suma) / matriz[i][i]
    
    mostrar_matriz(x)

    return x 


def sustitucion_hacia_adelante(matriz_L:list, vector_b):
    n = len(vector_b)
    x = [0] * n

    for i in range(n):
        suma = 0

        for j in range(i):

            suma += matriz_L[i][j] * x[j]

        x[i] = (vector_b[i] - suma) / matriz_L[i][i]

    mostrar_matriz(x)
    return x


def armar_L(lista_multiplicadores:list, dimension:int):
 
    matriz_i =  np.identity(dimension)
    for j in range (dimension): 

        for i in range (dimension):

            if i >= 1 and matriz_i[i-1][j] != 0  :
                matriz_i[i][j] = lista_multiplicadores[0]

                lista_multiplicadores.pop(0)
    return matriz_i

  
def main() -> None:
    "Tenemos una precision de 15 digitos"
    matriz = crear_matriz(4)
    soluciones_del_sistema = calcular_x(4)
   
    vector_b = hallar_b(matriz, soluciones_del_sistema)
   

    mat_ab = crear_matriz_aumentada(matriz,vector_b) #Concatenamos A y b para escalonar toda la mat junta
    lista_multiplicadores = aplicar_gauss(mat_ab) # Aplicamos gauss y guardamos los multiplicadores


    #mostrar_matriz(lista_multiplicadores)

    
    #U, y = np.hsplit(mat_ab, [-1]) # separamos a U y Y siendo U nuestra mas triangulada y Y nuestro sol indep
    mat_L = armar_L(lista_multiplicadores, 4) #Creamos la matriz L
    #mostrar_matriz(mat_L)
   

    #x = sustitucion_hacia_adelante(mat_L, vector_b) #Sustitucion hacia adelante para la mat L
    #sustitucion_hacia_atras(U,y)   Sustitucion hacia atras para la mat U
   

    #mostrar_matriz(mat_ab)
   
main()