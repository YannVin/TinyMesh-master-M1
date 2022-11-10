// Cylindre

#ifndef CYLINDRE_H
#define CYLINDRE_H

#pragma once

//includes

#include <vector>
#include <iostream>

#include "mathematics.h"

//class

class Cylindre {

public :


    //constructeur par défaut
    Cylindre();

    Cylindre(Vector,Vector,double);

    double ray; //taille du rayon
    Vector vec; //Vecteur position du cylindre
    Vector hau; //hauteur du cylindre
};



#endif // CYLINDRE_H
