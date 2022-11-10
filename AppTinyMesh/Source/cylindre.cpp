//Cylindre

//include
#include "cylindre.h"

//constructeur par défaut
Cylindre::Cylindre() {
    ray = 1; //taille du rayon du cylindre
    vec = Vector(); //Vecteur position du cylindre
    hau = Vector(0,0,1); //hauteur du cylindre
}

Cylindre::Cylindre(Vector v, Vector h, double r)
{
    ray = r; //taille du rayon du cylindre
    vec = v; //Vecteur position du cylindre
    hau = h; //hauteur du cylindre
}

