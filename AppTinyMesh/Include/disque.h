// Disque

#ifndef DISQUE_H
#define DISQUE_H


#pragma once

//includes

#include <vector>
#include <iostream>

#include "mathematics.h"

//class

class Disque {

public :


    //constructeur par défaut
    Disque();
    Disque(Vector, double);

    double ray; //taille du rayon
    Vector vec; //Vecteur position du disque
};


#endif // DISQUE_H

