import re
import numpy # type: ignore
from numpy.polynomial import polynomial as P # type: ignore

archivo = open("polinomio.txt")     #Abrimos el archivo de texto
contenidoArchivo = archivo.read()   #Leemos y guardamos lo que tenga
contenidoArchivo += " "     #Anadimos un espacio de mas para que almacenene el ultimo numero
print("Polinomio inicial, extremo izquierdo y derecho : " + contenidoArchivo)  #Imprimimos lo que tenia el archivo

#Subrutinas
def SeparaContenido():  #Para identificar el polinomio a analizar y el intervalo que se requiere
    SeparadorPolinomio = contenidoArchivo.split()   #Separa en un tipo lista el polinomio y el intervalo, en [0] -> polinomio, en [1] -> intervalo
    polinomio = SeparadorPolinomio[0]   #Se almacena en un string el polinomio
    SeparadorIntervalo = SeparadorPolinomio[1].split(",")   #En [1] esta el intervalo separado por "," asi que se separa de igual forma, indicando ahora que el separador sera la misma ","
    a = float(SeparadorIntervalo[0])    #Se almacena en tipo float el extremo izquierdo
    b = float(SeparadorIntervalo[1])    #Se almacena en tipo float el extremo derecho
    DefinePolinomio(polinomio,a,b)

def DefinePolinomio(polinomio,a,b):     #Procesara el polinomio, almacenando sus coeficientes y sus exponentes
    ExpresionRegular = r"(-?\d*\.?\d*)(x?)(?:(?:\^)(\d))?"  #ER que identifica de a "un termino" del polinomio
    coeficientes = []   #Almacenara los coeficientes del polinomio
    exponentes = []     #Almacenara los exponentes del polinomio
    for coef, x, exp in re.findall(ExpresionRegular, polinomio):    #Se tratara cada caso por separado
        if not coef and not x:
            continue
        if x and not coef:  #para el caso explicito de -> +x agregar el coeficiente como ->1
            coef = '1'
        if x and coef == "-":   #para el caso explicito de -> -x agregar el coeficiente como ->-1
            coef = "-1"
        if x and not exp:   #para el caso explicito de -> x^1 agregar el exponente como ->1
            exp = '1'
        if coef and not x:  #para identificar el ultimo coeficiente con grado 0
            exp = '0'
        coeficientes.append(float(coef))    #Almacenar los coeficientes
        exponentes.append(float(exp))   #Almacenar exponentes
    Derivada(coeficientes,exponentes,a,b)

def Derivada(coeficientes,exponentes,a,b):
    coeficientes.reverse()
    exponentes.reverse()
    Diferenciacion = P.polyder(coeficientes)
    CuadradoPolinomio = P.polypow(Diferenciacion,2)
    CuadradoPolinomio[0] = CuadradoPolinomio[0]+1
    ReglaTrapecio(CuadradoPolinomio,a,b)

def ReglaTrapecio(CuadradoPolinomio,a,b):
    Flip = numpy.flip(CuadradoPolinomio,0)
    X0 = numpy.polyval(Flip,a)
    RaizX0 = numpy.sqrt(X0)
    Xn = numpy.polyval(Flip,b)
    RaizXn = numpy.sqrt(Xn)
    DeltaX = (b-a)/500000   #500
    Inc = DeltaX/2
    Sumatoria = 0
    for contador in range(1,500000):    #500
        xi = a + (contador*DeltaX)
        Evaluacion = numpy.polyval(Flip,xi)
        RaizEv = numpy.sqrt(Evaluacion)
        Sumatoria = Sumatoria + RaizEv
    Sigma = 2*Sumatoria
    Corchete = RaizX0+RaizXn+Sigma
    Trapecio = Inc*Corchete
    print("La longitud de Arco: ",Trapecio)

SeparaContenido()
