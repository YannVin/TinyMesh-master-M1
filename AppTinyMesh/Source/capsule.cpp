// Capsule

//include
#include "capsule.h"

//constructeur par défaut
Capsule::Capsule() {
    vec = Vector(); //Vecteur position de la Capsule
    lon = 5; //Longueur de la capsule
    lar = 3; //Largeur de la capsule
}

Capsule::Capsule(Vector v, double lo, double la) {
    vec = v; //Vecteur position de la Capsule
    lon = lo; //Longueur de la capsule
    lar = la; //Largeur de la capsule
}

