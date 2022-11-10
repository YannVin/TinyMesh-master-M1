// Tore

#ifndef TORE_H
#define TORE_H
#pragma once

//includes

#include <vector>
#include <iostream>

#include "mathematics.h"

//class

class Tore {

public :


    //constructeur par défaut
    Tore();
    Tore(Vector,double,double);

    Vector vec; //Vecteur position du tore
    double raytube; //rayon du tube
    double raycercle; //rayon du cercle

};

#endif // TORE_H

