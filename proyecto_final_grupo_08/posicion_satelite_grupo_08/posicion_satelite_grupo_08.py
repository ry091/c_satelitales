import math #libreria usada para realizar operaciones matematicas (sen, cos,tan)
import datetime #libreria usada para manejar fechas y horas


altitud= 520/1000.0 # altura fija para santiago en km

def tiempo_siderio(hora_local): #función que calcula el tiempo sideral en radianes
    a = int((14 - hora_local.month) / 12) # Ajuste para el cálculo del año juliano
    y = hora_local.year + 4800 - a #año juliano
    m = hora_local.month + 12 * a - 3 #mes juliano
    
    jdn = (hora_local.day + int((153 * m + 2) / 5) + 365 * y + int(y / 4) - int(y / 100) + int(y / 400) - 32045) #calculo del número de día juliano
    fraccion = (hora_local.hour + hora_local.minute/60.0 +(hora_local.second + hora_local.microsecond/1000000.0)/3600.0) / 24.0 #fracción del día en decimal, para mayor precisión
    jd = jdn + fraccion - 0.5 #calculo del dia juliano
    d = jd - 2451545.0 #ajuste para el calculo, desde el 1-1-2000
    T = d / 36525.0 #tiempo en siglos julianos desde el 1-1-2000
    #calculo del tiempo siderio en segundos
    gmst = (24110.54841 + 8640184.812866*T + 0.093104*T*T - 6.2e-6*T*T*T +(hora_local.hour * 3600 + hora_local.minute * 60 + hora_local.second + hora_local.microsecond/1000000.0) * 1.00273790935)
    gmst = gmst % 86400 #asegurarse de que el tiempo sideral esté en el rango de 0 a 86400 segundos
    if gmst < 0: #cso de ajuste, si el tiempo es negativo
        gmst += 86400 #ajuste para que sea positivo
        
    gmst_rad = gmst * 2 * math.pi / 86400 #convertir a radianes
    return gmst_rad #retorna el tiempo siderio en radianes

def anomalia_verdadera(n,t,M_0,epoca,e): #funcion que calcula la anomalía verdadera de un satélite
    tiempo_dias = (t - epoca).total_seconds() / 86400 # convertir el tiempo a días
    M = (M_0 + n * tiempo_dias) % 360  # Anomalía media en grados
    M_rad = math.radians(M) #convertir a radianes
    
    E = M_rad #anomalía excéntrica inicial
    aux = 1 #variable aux para usar el metodo de Newton-Raphson
    i = 0 #n° de iteraciones
    while not (aux < 10e-8): #formula de convergencia de Newton-Raphson
        E_x = E - ((E - e * math.sin(E) - M_rad) / (1 - e * math.cos(E))) #calculo de la nueva anomalia
        aux = abs(E_x - E) # ajuste de diferencia entre la nueva y la anterior anomalía
        E = E_x # actualizar E
        i += 1 #contar iteración
    
    V_1 = math.sqrt((1 + e) / (1 - e)) * math.tan(E / 2) # calcular la anomalía verdadera
    verdadera = 2 * math.atan(V_1) % (2 * math.pi)  # asegurarse de que esté en el rango de 0, 2pi
    return verdadera #retorna la anomalía verdadera en radianes

def d_radial(verdadera, e, n): #funcion que calcula la distancia radial del satélite
    mu = 398600.4418 #constante de kepler, en km^3/s^2
    n_rad = n * 2 * math.pi / 86400  # Movimiento medio en radianes por segundo
    a = (mu / (n_rad ** 2)) ** (1 / 3) #obtención del semieje mayor en km
    r = (a * (1 - e**2)) / (1 + e * math.cos(verdadera)) # distancia radial en km
    return r #retorna la distancia radial en km


def coordenadas_satelite(r, verdadera, nodo_ascendente, inclinacion, argumento_perigeo): #funcion que calcula el vector de posición del satelite en sistemas eci
    x_pqw =r*math.cos(verdadera) # coordenadas en el plano PQW- inicialmente para x
    y_pqw =r *math.sin(verdadera) # coordenadas en el plano PQW- inicialmente para y
    z_pqw =0 # coordenadas en el plano PQW- inicialmente para z
    
    # argumentos trigonometricos usados para facilitar el calculo
    cos_omega=math.cos(argumento_perigeo) #cos del argumento del perigeo
    sin_omega=math.sin(argumento_perigeo) #sin del argumento del perigeo
    cos_i=math.cos(inclinacion) #cos de la inclinación
    sin_i=math.sin(inclinacion) #sin de la inclinación
    cos_nodo_asc=math.cos(nodo_ascendente) #cos del nodo ascendente
    sin_nodo_asc=math.sin(nodo_ascendente) #sin del nodo ascendente

    # realizar la matriz de transformaciones, con proceso de rotación- derivando 3 veces de acuerdo  a PQW
    r1_1=cos_nodo_asc* cos_omega-sin_nodo_asc *sin_omega*cos_i # 1° rotacion (1 derivada de la matriz de rotación)
    r1_2=-cos_nodo_asc* sin_omega- sin_nodo_asc*cos_omega*cos_i # 1° rotacion (1 derivada de la matriz de rotación)
    r1_3=sin_nodo_asc* sin_i # 1° rotacion (1 derivada de la matriz de rotación)

    r2_1= sin_nodo_asc* cos_omega +cos_nodo_asc*sin_omega *cos_i # 2° rotacion (2 derivada de la matriz de rotación)
    r2_2=-sin_nodo_asc*sin_omega + cos_nodo_asc *cos_omega*cos_i # 2° rotacion (2 derivada de la matriz de rotación)
    r2_3=-cos_nodo_asc* sin_i # 2° rotacion (2 derivada de la matriz de rotación)

    r3_1 = sin_omega *sin_i # 3° rotacion (3 derivada de la matriz de rotación)
    r3_2 = cos_omega *sin_i # 3° rotacion (3 derivada de la matriz de rotación)
    r3_3 = cos_i # 3° rotacion (3 derivada de la matriz de rotación)

    # Aplicar transformación de coordenadas de PQW a ECI
    x_eci=r1_1*x_pqw+r1_2 * y_pqw +r1_3* z_pqw # coordenadas ECI para x
    y_eci= r2_1* x_pqw+ r2_2 * y_pqw+ r2_3* z_pqw #cooredenadas ECI para y
    z_eci=r3_1* x_pqw +r3_2 * y_pqw+r3_3 *z_pqw #coordenadas ECI para z
    
    return x_eci, y_eci, z_eci #retorna las coordenadas ECI del satélite en km

def coordenadas_observador(latitud, longitud, GMST_rad):  #función que calcula las coordenadas ECI del observador
    #uso de la orbita WGS84 para el calculo de las coordenadas ECEF y tener mayor precisión
    a = 6378.137  # parámetro emieje mayor en km
    f = 1/298.257223563  # parámetro de achatamiento
    e2 = f * (2 - f)  # excentricidad al cuadrado
    lat_rad = math.radians(latitud) #converción de latitud a radianes
    lon_rad = math.radians(longitud)  #converción de longitud a radianes
    N =a/ math.sqrt(1 - e2 * math.sin(lat_rad)**2) #radio de curvatura en el primer vertical

    x_ecef =(N+altitud)*math.cos(lat_rad)*math.cos(lon_rad) #coordenadas ECEF para x
    y_ecef =(N+altitud)*math.cos(lat_rad)*math.sin(lon_rad) #coordenadas ECEF para y
    z_ecef =(N*(1 - e2)+altitud)*math.sin(lat_rad) #coordenadas ECEF para z

    # Rotación ECEF → ECI
    x_eci=x_ecef*math.cos(GMST_rad)-y_ecef*math.sin(GMST_rad) # coordenadas ECI para x
    y_eci=x_ecef*math.sin(GMST_rad)+y_ecef*math.cos(GMST_rad) # coordenadas ECI para y
    z_eci=z_ecef  #coordenadas ECI para z
    
    return x_eci, y_eci, z_eci #retorna las coordenadas ECI del observador en km

def sistema_local(vector_relativo, latitud, longitud, GMST_rad): #funcion que calcula el sistema local del observador
    x, y,z = vector_relativo #coordernadas del vector relativo entre el satélite y el observador (v_satelite - v_observador)
    lat_rad = math.radians(latitud) #converción de latitud a radianes
    lon_rad = math.radians(longitud)#converción de longitud a radianes
    H = (GMST_rad +lon_rad)%(2 * math.pi) #hora sideral a radianes
      #transformación de coordenadas ECI a sistema local SEZ, uso de matriz respecto al vector relativo
    s_sez= x*(math.sin(lat_rad)*math.cos(H)) + y*(math.sin(lat_rad)*math.sin(H))- z*(math.cos(lat_rad) ) #coordenada sur en el sistema SEZ
    e_sez= x*(-math.sin(H)) + y*(math.cos(H))#coordenada este en el sistema SEZ
    z_sez=x*(math.cos(lat_rad)* math.cos(H)) + y*(math.cos(lat_rad)*math.sin(H)) + z*(math.sin(lat_rad))#coordenada zenith en el sistema SEZ

    elevaccion= math.asin(z_sez/math.sqrt(s_sez**2 + e_sez**2 + z_sez**2)) #obtencion de la elevación en radianes
    azimut= math.atan2(e_sez,-s_sez) #obtencion del azimut en radianes
    azimut_grados= math.degrees(azimut) % 360 #conversion azimut a grados con ajuste
    elevaccion_grados= math.degrees(elevaccion) #conversion elevación a grados

    print(f"Azimut: {azimut_grados:.4f}°, Elevación: {elevaccion_grados:.4f}°"    ) #imprime los resultados obtenidos
    return 0 #se acaba la función sin retorno

def main():
    with open("tle.txt", "r") as f: #lee un archivo tle.txt dentro del mismo directorio
        lineas = [linea.strip() for linea in f.readlines()] #separa los resultados por linea

    linea1= lineas[0] #obtiene la primera linea que esta en el archivo
    parte1 = linea1.split() #separa la primera linea por espacios
    epoca = (parte1[3]) #obtiene la epoca del satelite, del tle
    #separación de la epoca en sus partes
    epoca_a=int(epoca[:2]) + 2000 # obtención del año, se asume que es desde el siglo 21 y se le suma 2000
    epoca_d= int(epoca[2:5]) #obtencion del dia juliano, desde el 1 de enero
    epoca_fraccion = float(epoca[5:]) #obtencion de la fracción del dia, en formato decimal
    epoca_dt= datetime.datetime(epoca_a,1,1) + datetime.timedelta(days=epoca_d - 1 + epoca_fraccion) #obtención de la fecha y hora de la epoca, en formato datetime
    linea2 = lineas[1] #obtiene la segunda linea que esta en el archivo del tle
    parte2 = linea2.split() #separa la segunda linea por espacios
    inclinacion = math.radians(float(parte2[2])) #obtiene la inclinacion del satelite y lo convierte a radianes
    nodo_ascendente = math.radians(float(parte2[3])) #obtiene el nodo ascendente del satelite y lo convierte a radianes
    excentricidad = float("0."+parte2[4]) #obtiene la excentricidad del satelite, y se realiza un casteo
    argumento_perigeo = math.radians(float(parte2[5])) #obtiene el argumento del perigeo del satelite y lo convierte a radianes
    anomalia_media_0 = float(parte2[6]) #obtiene la anomalía media inicial del satélite
    movimiento_medio = float(parte2[7])#obtiene el movimiento medio del satélite en revoluciones por día
    hora_utc = datetime.datetime.utcnow()# obtiene la hora UTC actual
    fecha_local = hora_utc + datetime.timedelta(hours=-4) #ajusta la hora UTC a la hora local en chile
    hora_local = fecha_local.strftime("%H:%M:%S") #muestra solo la hora
    latitud= float(lineas[2]) #obtiene la latitud de la tercera linea del archivo tle
    longitud= float(lineas[3]) #obtiene la longitud de la cuarta linea del archivo tle
    GMST_rad = tiempo_siderio(hora_utc)#obtiene la hora sideral en radianes
    print(f"hora local: {hora_local}") #imprime la hora local


    a_verdadera = anomalia_verdadera(movimiento_medio * 360, hora_utc, anomalia_media_0, epoca_dt, excentricidad) #calcula la anomalía verdadera del satélite
    distancia_radial= d_radial(a_verdadera, excentricidad, movimiento_medio) #calcula la distancia radial del satélite
    x_sat,y_sat,z_sat= coordenadas_satelite(distancia_radial, a_verdadera, nodo_ascendente, inclinacion, argumento_perigeo) #calcula las coordenadas del satelite
    x_obs, y_obs, z_obs = coordenadas_observador(latitud, longitud, GMST_rad) #calcula las coordenadas del observador

    vector_relativo = [x_sat - x_obs, y_sat - y_obs, z_sat - z_obs] #resta ambos vectores (satelite - observador) para obtener el vector relativo
    sistema_local(vector_relativo, latitud, longitud, GMST_rad) #obtiene el azimut y elevación del satélite respecto al observador

if __name__ == "__main__":
    main()
