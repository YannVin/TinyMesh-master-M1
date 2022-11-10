// HeightField

#ifndef HEIGHTFIELD_H
#define HEIGHTFIELD_H


#include <vector>
#include <iostream>
#include <QImage>
#include "mathematics.h"

//class

class HeightField {

public :


    //constructeur par défaut
    HeightField();
    HeightField(QImage,int);

    Vector vec; //Vecteur position du HeightField
    QImage image; //image en niveaux de gris
    int echelle; //echelle pour limiter la taille globale

};


#endif // HEIGHTFIELD_H
