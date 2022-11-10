// Sphere

//include
#include "sphere.h"

//constructeur par défaut
Sphere::Sphere() {
    vec = Vector(); //Vecteur position de la sphere
    rayl = 5;  //taille du rayon de largeur
    rayh = 5;  //taille du rayon en hauteur
}

Sphere::Sphere(Vector v, double rl, double rh)
{
    vec = v; //Vecteur position de la sphere
    rayl = rl;  //taille du rayon de largeur
    rayh = rh;  //taille du rayon en hauteur
}

