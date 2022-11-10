// Capsule

#ifndef CAPSULE_H
#define CAPSULE_H

#pragma once

//includes

#include <vector>
#include <iostream>

#include "mathematics.h"

//class

class Capsule {

public :


    //constructeur par défaut
    Capsule();
    Capsule(Vector, double, double);

    Vector vec; //Vecteur position de la Capsule
    double lon; //Longueur de la capsule
    double lar; //Largeur de la capsule
};



#endif // CAPSULE_H
