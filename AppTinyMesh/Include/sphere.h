// Sphere

#ifndef SPHERE_H
#define SPHERE_H

#pragma once

//includes

#include <vector>
#include <iostream>

#include "mathematics.h"

//class

class Sphere {

public :


    //constructeur par défaut
    Sphere();
    Sphere(Vector, double, double);

    Vector vec; //Vecteur position de la sphere
    double rayl; //taille du rayon de largeur
    double rayh; //taille du rayon de hauteur
};

#endif // SPHERE_H
