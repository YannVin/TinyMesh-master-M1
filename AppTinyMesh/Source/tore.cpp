// Tore

//include
#include "tore.h"

//constructeur par défaut
Tore::Tore() {
    vec = Vector(); //Vecteur position du tore
    raytube = 1;  //rayon du tube
    raycercle = 4; //rayon du cercle
}

Tore::Tore(Vector v ,double rt ,double rc) {
    vec = v; //Vecteur position du tore
    raytube = rt;  //rayon du tube
    raycercle = rc; //rayon du cercle
}

