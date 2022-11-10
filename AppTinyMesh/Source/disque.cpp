//Disque

//include
#include "disque.h"

//constructeur par défaut
Disque::Disque() {
    ray = 10; //taille du rayon du disque
    vec = Vector(); //Vecteur position du disque
}

Disque::Disque(Vector v, double r){
    ray = r; //taille du rayon du disque
    vec = v; //Vecteur position du disque
}
