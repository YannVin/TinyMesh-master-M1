// HeightField

//include
#include "heightField.h"

//constructeur par défaut
HeightField::HeightField() {
    vec = Vector(-image.width()/2,-image.height()/2,0); //Vecteur position du HeightField
    image = QImage("C:/Users/ilies/OneDrive/Documents/Texture_terrain.png","png"); //image en niveaux de gris
    echelle = 10; //echelle pour limiter la taille globale
}

HeightField::HeightField(QImage a, int e) {
    image = a; //image en niveaux de gris
    vec = Vector(-image.height()/2,-image.width()/2,0); //Vecteur position du HeightField
    echelle = e; //echelle pour limiter la taille globale
}

