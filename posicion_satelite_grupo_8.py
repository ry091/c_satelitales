import math
import datetime


def tiempo_siderio(hora_local):
    ano= hora_local.year
    mes= hora_local.month
    dia= hora_local.day
    hora= hora_local.hour
    minuto= hora_local.minute
    segundo= hora_local.second
    if mes <= 2:
        ano -= 1
        mes += 12
    
    jd= 367 * ano - int((7 * (ano + int((mes + 9) / 12))) // 4) + int((275 * mes) / 9) + dia + 1721013.5 + (hora + minuto / 60 + segundo / 3600) / 24
    D = jd - 2451545.0  # desde 1-2000
    
    GMST_grados = (280.46061837 + 360.98564736629 * D) % 360

    return math.radians(GMST_grados)

def anomalia_verdadera(n,t,M_0,epoca,e):
    M= (float(M_0 +n*((t - epoca).total_seconds()/86400))) %360 #valor de la anomalia media en grados
    M_rad = math.radians(M) 
    print(f"anomalia media: {M_rad}")

    #anomalia excentrica
    E= M_rad
    aux= 1
    i=0
    while not (aux< 10e-8):
        E_x= E - ((E - e * math.sin(E) - M_rad) / (1 - e * math.cos(E)))
        aux= abs(E_x -E)
        E= E_x
        i=i+1
        #print(f"valor anomalia: {E}  ; valor tramo :{aux}")

    print(f"valor anomala excentrica: {E}")

    #anomalia verdadera
    V_1 = math.sqrt((1+e)/(1-e)) * math.tan(E / 2)
    verdadera = 2 * math.atan(V_1) % (2 * math.pi)  # normalizar a [0, 2π)
    print(f"anomalia verdadera: {verdadera}")
    return verdadera

def d_radial(verdadera,e,n):
    periodo= 86400/ n # periodo en segundos
    a= ((math.pow(periodo,2)* 3.9857128e14)/(4*math.pow(math.pi,2)))**(1/3) # semieje mayor (a)
    r= (a*(1-math.pow(e,2)))/(1 + e * math.cos(verdadera)) # distancia radial
    return r

def cordenadas_satelite(r, verdadera, nodo_ascendente, inclinacion, argumento_perigeo):
    x = (r * (math.cos(nodo_ascendente) * math.cos(verdadera + argumento_perigeo) - math.sin(nodo_ascendente) * math.sin(verdadera + argumento_perigeo) * math.cos(inclinacion))) /1000
    y = (r * (math.sin(nodo_ascendente) * math.cos(verdadera + argumento_perigeo) + math.cos(nodo_ascendente) * math.sin(verdadera + argumento_perigeo) * math.cos(inclinacion)))/1000
    z = (r * (math.sin(inclinacion) * math.sin(verdadera + argumento_perigeo)))/1000

    print(f"coordenadas del satelite: x={x}, y={y}, z={z}")
    return x, y, z

def coordenadas_observador(latitud,longitud,GMST_rad):
    #uso del modelo WGS84 para coordenadas geodésicas
    altitud = 520 /1000 #santiago, Chile
    a= 6378.137  # semieje mayor en km
    f= 1/298.257223563  # achatamiento
    e2= f * (2 - f)  # excentricidad al cuadrado
    lat_rad = math.radians(latitud)
    lon_rad = math.radians(longitud)
    N = a / math.sqrt(1 - e2 * math.sin(lat_rad)**2)  # radio de curvatura en el meridiano
    x = (N + altitud) * math.cos(lat_rad) * math.cos(lon_rad)
    y = (N + altitud) * math.cos(lat_rad) * math.sin(lon_rad)
    z = ((1 - e2) * N + altitud) * math.sin(lat_rad)

    x_eci= x * math.cos(GMST_rad) - y * math.sin(GMST_rad)
    y_eci= x * math.sin(GMST_rad) + y * math.cos(GMST_rad)
    z_eci = z  

    print(f"coordenadas del observador: x={x_eci}, y={y_eci}, z={z_eci}")

    return x_eci, y_eci, z_eci


def sistema_local(vector_relativo, latitud, longitud, GMST_rad):
    x= vector_relativo[0]
    y= vector_relativo[1]
    z= vector_relativo[2]
    lat_rad = math.radians(latitud)
    lon_rad = math.radians(longitud)
    H= (GMST_rad + lon_rad) % (2 * math.pi)  # hora sideral en radianes
    # rotacion sistema sez
    sez_s = (-math.sin(lat_rad) * math.cos(H)) * x + (-math.sin(lat_rad) * math.sin(H) )* y + (math.cos(lat_rad) )* z
    sez_e= (-math.sin(H))*x + (math.cos(H))*y
    sez_z= (math.cos(lat_rad) * math.cos(H))*x + (math.cos(lat_rad) * math.sin(H))*y + (math.sin(lat_rad))*z

    #calculo de azimut y elevacion
    elevacion = math.asin(sez_z / math.sqrt(sez_s**2 + sez_e**2 + sez_z**2))
    azimut = math.atan2(sez_e, -sez_s)

    elevacion_grados = math.degrees(elevacion)
    azimut_grados = math.degrees(azimut)% 360
    print(f"Azimut: {azimut_grados:.2f} grados, Elevacion: {elevacion_grados:.2f} grados")



def main():
    with open("tle.txt", "r") as f:
        lineas = [linea.strip() for linea in f.readlines()]

    linea1= lineas[0]
    parte1 = linea1.split()
    epoca = (parte1[3])
    epoca_a=int(epoca[:2]) + 2000 # se asume desde el año 2000
    epoca_d= int(epoca[2:5]) 
    epoca_fraccion = float(epoca[5:]) 
    epoca_dt= datetime.datetime(epoca_a,1,1) + datetime.timedelta(days=epoca_d - 1 + epoca_fraccion)
    print(f"epoca: {epoca_dt}")

    linea2 = lineas[1]
    latitud= float(lineas[2])
    longitud= float(lineas[3])
    parte2 = linea2.split()

    inclinacion = math.radians(float(parte2[2]) ) 
    nodo_ascendente = math.radians(float(parte2[3]) )
    excentricidad = float("0."  + parte2[4])
    argumento_perigeo = math.radians(float(parte2[5]))
    anomalia_media_0 = float(parte2[6])
    movimiento_medio = float(parte2[7] )
    hora_local = datetime.datetime.utcnow()
    GMST_rad = tiempo_siderio(hora_local)

   

    print(f"hora local: {hora_local}")

    a_verdadera=anomalia_verdadera(movimiento_medio *360, hora_local, anomalia_media_0, epoca_dt,excentricidad)
    distancia_radial= d_radial(a_verdadera, excentricidad, movimiento_medio)

    x_sat,y_sat,z_sat= cordenadas_satelite(distancia_radial, a_verdadera, nodo_ascendente, inclinacion, argumento_perigeo)
    x_obs, y_obs, z_obs = coordenadas_observador(latitud, longitud, GMST_rad)

    vector_relativo = [x_sat - x_obs, y_sat - y_obs, z_sat - z_obs] 
    sistema_local(vector_relativo, latitud, longitud, GMST_rad,)






if __name__ == "__main__":
    main()
