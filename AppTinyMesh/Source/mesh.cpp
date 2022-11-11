#include "mesh.h"

/*!
\class Mesh mesh.h

\brief Core triangle mesh class.
*/



/*!
\brief Initialize the mesh to empty.
*/
Mesh::Mesh()
{
}

/*!
\brief Initialize the mesh from a list of vertices and a list of triangles.

Indices must have a size multiple of three (three for triangle vertices and three for triangle normals).

\param vertices List of geometry vertices.
\param indices List of indices wich represent the geometry triangles.
*/
Mesh::Mesh(const std::vector<Vector>& vertices, const std::vector<int>& indices) :vertices(vertices), varray(indices)
{
  normals.resize(vertices.size(), Vector::Z);
}

/*!
\brief Create the mesh.

\param vertices Array of vertices.
\param normals Array of normals.
\param va, na Array of vertex and normal indexes.
*/
Mesh::Mesh(const std::vector<Vector>& vertices, const std::vector<Vector>& normals, const std::vector<int>& va, const std::vector<int>& na) :vertices(vertices), normals(normals), varray(va), narray(na)
{
}

/*!
\brief Reserve memory for arrays.
\param nv,nn,nvi,nvn Number of vertices, normals, vertex indexes and vertex normals.
*/
void Mesh::Reserve(int nv, int nn, int nvi, int nvn)
{
  vertices.reserve(nv);
  normals.reserve(nn);
  varray.reserve(nvi);
  narray.reserve(nvn);
}

/*!
\brief Empty
*/
Mesh::~Mesh()
{
}

/*!
\brief Smooth the normals of the mesh.

This function weights the normals of the faces by their corresponding area.
\sa Triangle::AreaNormal()
*/
void Mesh::SmoothNormals()
{
  // Initialize 
  normals.resize(vertices.size(), Vector::Null);

  narray = varray;

  // Accumulate normals
  for (int i = 0; i < varray.size(); i += 3)
  {
    Vector tn = Triangle(vertices[varray.at(i)], vertices[varray.at(i + 1)], vertices[varray.at(i + 2)]).AreaNormal();
    normals[narray[i + 0]] += tn;
    normals[narray[i + 1]] += tn;
    normals[narray[i + 2]] += tn;
  }

  // Normalize 
  for (int i = 0; i < normals.size(); i++)
  {
    Normalize(normals[i]);
  }
}

/*!
\brief Add a smooth triangle to the geometry.
\param a, b, c Index of the vertices.
\param na, nb, nc Index of the normals.
*/
void Mesh::AddSmoothTriangle(int a, int na, int b, int nb, int c, int nc)
{
  varray.push_back(a);
  narray.push_back(na);
  varray.push_back(b);
  narray.push_back(nb);
  varray.push_back(c);
  narray.push_back(nc);
}

/*!
\brief Add a triangle to the geometry.
\param a, b, c Index of the vertices.
\param n Index of the normal.
*/
void Mesh::AddTriangle(int a, int b, int c, int n)
{
  varray.push_back(a);
  narray.push_back(n);
  varray.push_back(b);
  narray.push_back(n);
  varray.push_back(c);
  narray.push_back(n);
}

/*!
\brief Add a smmoth quadrangle to the geometry.

Creates two smooth triangles abc and acd.

\param a, b, c, d  Index of the vertices.
\param na, nb, nc, nd Index of the normal for all vertices.
*/
void Mesh::AddSmoothQuadrangle(int a, int na, int b, int nb, int c, int nc, int d, int nd)
{
  // First triangle
  AddSmoothTriangle(a, na, b, nb, c, nc);

  // Second triangle
  AddSmoothTriangle(a, na, c, nc, d, nd);
}

/*!
\brief Add a quadrangle to the geometry.

\param a, b, c, d  Index of the vertices and normals.
*/
void Mesh::AddQuadrangle(int a, int b, int c, int d)
{
  AddSmoothQuadrangle(a, a, b, b, c, c, d, d);
}

/*!
\brief Compute the bounding box of the object.
*/
Box Mesh::GetBox() const
{
  if (vertices.size() == 0)
  {
    return Box::Null;
  }
  return Box(vertices);
}

/*!
\brief Fonction qui concatene 2 mesh.
\param m Mesh, Le mesh que l'on veux concatener.
*/
void Mesh::Merge(const Mesh& m)
{
    for(int i= 0; i< m.NormalIndexes().size(); i++) //on ajoute au tableau narray les futures indices des points du mesh qu'on ajoute.
    {
        narray.push_back(Normales() + m.NormalIndexes()[i]); //on ajoute a tous les indices la taille du tableau pour pointer sur la partie que l'on ajoute au tableau normals (sans l'ajout de de Normales() les triangle du 2ieme mesh ne seront jamais fait car les indices de ses vecteur ne serait jamais pointé)
    }

    for(int i= 0; i< m.VertexIndexes().size(); i++) //meme chose pour varray
    {
        varray.push_back(Vertexes() + m.VertexIndexes()[i]);
    }

    // on concatene les tableau de vecteurs et de normal
    for(int i= 0; i< m.Normales(); i++)
    {
        normals.push_back(m.Normal(i));
    }

    for(int i= 0; i< m.Vertexes(); i++)
    {
        vertices.push_back(m.Vertex(i));
    }
}

/*!
\brief Fonction de deformation d'une sphere.
\param s Sphere de deformation, x Vecteur de deforamtion.
*/
void Mesh::Warp(const Sphere& s,Vector x)
{
    Vector v;
    for(int i=0;i<vertices.size();i++)
    {
        double d= SquaredNorm(vertices[i]-s.vec);

        if(d<(s.rayl*s.rayl))
        {
            v = x * (1-(d/(s.rayl*s.rayl))) ;
            vertices[i]= vertices[i] + v;
        }
    }
}


/*!
\brief Creates an axis aligned box.

The object has 8 vertices, 6 normals and 12 triangles.
\param box The box.
*/
Mesh::Mesh(const Box& box)
{
  // Vertices
  vertices.resize(8);

  for (int i = 0; i < 8; i++)
  {
    vertices[i] = box.Vertex(i);
  }

  // Normals
  normals.push_back(Vector(-1, 0, 0));
  normals.push_back(Vector(1, 0, 0));
  normals.push_back(Vector(0, -1, 0));
  normals.push_back(Vector(0, 1, 0));
  normals.push_back(Vector(0, 0, -1));
  normals.push_back(Vector(0, 0, 1));

  // Reserve space for the triangle array
  varray.reserve(12 * 3);
  narray.reserve(12 * 3);

  AddTriangle(0, 2, 1, 4);
  AddTriangle(1, 2, 3, 4);

  AddTriangle(4, 5, 6, 5);
  AddTriangle(5, 7, 6, 5);

  AddTriangle(0, 4, 2, 0);
  AddTriangle(4, 6, 2, 0);

  AddTriangle(1, 3, 5, 1);
  AddTriangle(3, 7, 5, 1);

  AddTriangle(0, 1, 5, 2);
  AddTriangle(0, 5, 4, 2);

  AddTriangle(3, 2, 7, 3);
  AddTriangle(6, 7, 2, 3);
}


/*!
\brief Consctructeur Mesh du Disque.

\param d disque.
\param nbrpt nombre de point du disque.
*/
Mesh::Mesh(const Disque& d, const int& nbrpt)
{
    // Vertices
    vertices.resize(nbrpt + 2);
    float alpha;
    float step= 2.0 * M_PI / (nbrpt);
    vertices[0] = Vector()+d.vec;

    // Reserve space for the triangle array
    varray.reserve(nbrpt * 3);
    narray.reserve(nbrpt * 3);

    // Variation de l’angle de 0 à 2pi
    for (int i = 0; i <= nbrpt; i++)
    {
        alpha = i*step; //Evolution de l'angle
        vertices[i+1] = Vector(cos(alpha),0,sin(alpha))*d.ray + d.vec; // Sommet des coté du disque
        normals.push_back(Vector(cos(alpha),0,sin(alpha))*d.ray + d.vec); // Ajout des normals
        AddTriangle(0, i, i+1, i); //ajout de chaque triangle entre les points
    }

}


/*!
\brief Consctructeur Mesh du Cylindre.

\param c Cylindre.
\param nbrpt nombre de point du Cylindre.
*/
Mesh::Mesh(const Cylindre& c, const int& nbrpt) 
{
    // Vertices
    vertices.resize((nbrpt + 2)*2-1);
    float alpha;
    float step= 2.0 * M_PI / (nbrpt);
    vertices[nbrpt*2] = Vector() + c.vec;
    vertices[(nbrpt + 1)*2] = c.hau + c.vec;

    // Reserve space for the triangle array
    varray.reserve(nbrpt * 6);
    narray.reserve(nbrpt * 6);

    // Variation de l’angle de 0 à 2pi
    for (int i = 0; i < nbrpt; i++)
    {
        alpha = i*step; //Evolution de l'angle
        vertices[i] = Vector(cos(alpha)*c.ray,sin(alpha)*c.ray,0) + c.vec; // Sommet des coté du disque
        normals.push_back(Vector(cos(alpha),sin(alpha),0)); //normals
        AddTriangle(nbrpt*2, i, (i+1)%nbrpt, i*2); //ajout de chaque triangle entre les points

        vertices[i+nbrpt] = Vector(cos(alpha)*c.ray,sin(alpha)*c.ray,0)+ c.hau + c.vec; // Sommet des coté du disque 2
        normals.push_back(Vector(cos(alpha),sin(alpha),0)); //normals du disque 2
        AddTriangle((nbrpt + 1)*2, i+nbrpt, ((i+1)%nbrpt+nbrpt), i*2); //ajout de chaque triangle entre les points

        AddTriangle(i,       i+nbrpt,     (i+1)%nbrpt,           i*2); //ajout de chaque triangle entre les points des 2 disques
        AddTriangle(i+nbrpt, (i+1)%nbrpt, (((i+1)%nbrpt)+nbrpt), i*2); //ajout de chaque triangle entre les points des 2 disques
    }

}

/*!
\brief Consctructeur Mesh de la Sphere.

\param s Sphere.
\param nbrpt nombre de point du Sphere.
*/
Mesh::Mesh(const Sphere& s, const int& nbrpt)
{
    // Vertices
    const int divBeta= nbrpt; //Déclaration du nombre de division Beta
    const int divAlpha= divBeta/2; //Déclaration du nombre de division alpha
    float beta, alpha, alpha2; //Declaration de 3 angles

    vertices.resize((divAlpha+1) *(divBeta));

    // Reserve space for the triangle array
    varray.reserve((divAlpha+1) *divBeta * 6);
    narray.reserve((divAlpha+1) *divBeta * 6);

    // Variation des angles alpha et beta
    for(int i=0; i< divAlpha; i++)
    {
        alpha= -0.5f * M_PI + float(i) * M_PI / divAlpha; //evolution de l'angle alpha
        alpha2= -0.5f * M_PI + float(i+1) * M_PI / divAlpha; //evolution de l'angle alpha2

        for(int j=0; j< divBeta; j++)
        {
            beta= float(j) * 2.f * M_PI / (divBeta); //evolution de l'angle beta

            normals.push_back( Vector(cos(alpha)*cos(beta)*s.rayl,sin(alpha)*s.rayl, cos(alpha)*sin(beta)*s.rayh)); //Normale pour le bas
            vertices[i*divBeta+ j] = Vector(cos(alpha)*cos(beta)*s.rayl,sin(alpha)*s.rayl, cos(alpha)*sin(beta)*s.rayh) + s.vec; //Sommets des points des cercles du bas

            normals.push_back( Vector(cos(alpha2)*cos(beta)*s.rayl,sin(alpha2)*s.rayl, cos(alpha2)*sin(beta)*s.rayh)); //Normale pour le haut
            vertices[i*divBeta+ j+divBeta] = Vector(cos(alpha2)*cos(beta)*s.rayl,sin(alpha2)*s.rayl, cos(alpha2)*sin(beta)*s.rayh) + s.vec; //Sommets des points des cercles du haut

            //ajout de chaque triangle entre les 2 cercles qui se deplaces
            AddTriangle(i*divBeta+ j,             i*divBeta+ j+divBeta,i*divBeta+ (j+1)%divBeta,        (i*divBeta+ j)*2 );
            AddTriangle(i*divBeta+ (j+1)%divBeta, i*divBeta+ j+divBeta,i*divBeta+ (j+1)%divBeta+divBeta,(i*divBeta+ j+divBeta)*2 );

        }
    }
}

/*!
\brief Consctructeur Mesh de la Capsule.

\param c Capsule.
\param nbrpt nombre de point de la Capsule.
*/
Mesh::Mesh(const Capsule& c, const int& nbrpt) //Consctructeur Mesh de la sphere
{
    // Vertices
    const int divBeta= nbrpt; //Déclaration du nombre de division Beta
    const int divAlpha= divBeta/2; //Déclaration du nombre de division alpha
    float beta, alpha, alpha2; //Declaration de 3 angles

    vertices.resize((divAlpha+2)*(divBeta+2)-4 );

    // Reserve space for the triangle array
    varray.reserve((divAlpha+1) *divBeta * 6 );
    narray.reserve((divAlpha+1) *divBeta * 6 );


    //Creation d'un demi sphere (chapeau de la capsule)

    // Variation des angles alpha et beta
    for(int i=0; i< divAlpha/2; i++)
    {
        alpha= -0.5f * M_PI + float(i) * M_PI / divAlpha; //evolution de l'angle alpha
        alpha2= -0.5f * M_PI + float(i+1) * M_PI / divAlpha; //evolution de l'angle alpha2

        for(int j=0; j< divBeta; j++)
        {
            beta= float(j) * 2.f * M_PI / (divBeta); //evolution de l'angle beta

            normals.push_back( Vector(cos(alpha)*cos(beta),sin(alpha), cos(alpha)*sin(beta))); //Normale pour le bas
            vertices[i*divBeta+ j] = Vector(cos(alpha)*cos(beta)*c.lar,sin(alpha)*c.lar, cos(alpha)*sin(beta)*c.lar) + c.vec; //Sommets des points des cercles du bas

            normals.push_back( Vector(cos(alpha2)*cos(beta), sin(alpha2), cos(alpha2)*sin(beta))); //Normale pour le haut
            vertices[i*divBeta+ j+divBeta] = Vector(cos(alpha2)*cos(beta)*c.lar,sin(alpha2)*c.lar, cos(alpha2)*sin(beta)*c.lar) + c.vec; //Sommets des points des cercles du haut

            //ajout de chaque triangle entre les 2 cercles qui se deplaces
            AddTriangle(i*divBeta+ j,             i*divBeta+ j+divBeta,i*divBeta+ (j+1)%divBeta,        (i*divBeta+ j)*2 );
            AddTriangle(i*divBeta+ (j+1)%divBeta, i*divBeta+ j+divBeta,i*divBeta+ (j+1)%divBeta+divBeta,(i*divBeta+ j+divBeta)*2 );

        }
    }

    //Creation du deuxieme demi sphere (chapeau de la capsule). Les triangles pour relier les 2 demis sphere (cylindre) se feront automatiquement avec les AddTriangle. Les 2 mesh sont concaténés.

    // Variation des angles alpha et beta
    for(int i=divAlpha/2; i<=divAlpha; i++)
    {
        alpha= -0.5f * M_PI + float(i) * M_PI / divAlpha; //evolution de l'angle alpha
        alpha2= -0.5f * M_PI + float(i+1) * M_PI / divAlpha; //evolution de l'angle alpha2

        for(int j=0; j< divBeta; j++)
        {
            beta= float(j) * 2.f * M_PI / (divBeta); //evolution de l'angle beta

            normals.push_back( Vector(cos(alpha)*cos(beta),sin(alpha), cos(alpha)*sin(beta))); //Normale pour le bas
            vertices[i*divBeta+ j+divBeta] = Vector(cos(alpha)*cos(beta)*c.lar,sin(alpha)*c.lar+c.lon, cos(alpha)*sin(beta)*c.lar) + c.vec; //Sommets des points des cercles du bas

            normals.push_back( Vector(cos(alpha2)*cos(beta), sin(alpha2), cos(alpha2)*sin(beta))); //Normale pour le haut
            vertices[i*divBeta+ j+divBeta*2] = Vector(cos(alpha2)*cos(beta)*c.lar,sin(alpha2)*c.lar+c.lon, cos(alpha2)*sin(beta)*c.lar) + c.vec; //Sommets des points des cercles du haut (avec ajout de la taille de la capsule "c.lon")

            //ajout de chaque triangle entre les 2 cercles qui se deplaces
            AddTriangle(i*divBeta+ j,             i*divBeta+ j+divBeta,i*divBeta+ (j+1)%divBeta,        (i*divBeta+ j)*2 );
            AddTriangle(i*divBeta+ (j+1)%divBeta, i*divBeta+ j+divBeta,i*divBeta+ (j+1)%divBeta+divBeta,(i*divBeta+ j+divBeta)*2 );



        }
    }

}

/*!
\brief Consctructeur Mesh du Tore.

\param tore Tore.
\param nbrpt nombre de point du Tore.
*/
Mesh::Mesh(const Tore& tore,const int& nbpoint)
{
    // Vertices
    vertices.resize((nbpoint+1)*(nbpoint+1) );
    float alpha;
    float beta;
    float step= 2.0 * M_PI / (nbpoint);

    // Reserve space for the triangle array
    varray.reserve((nbpoint * 3)*nbpoint);
    narray.reserve((nbpoint * 3)*nbpoint);

    // Variation de l’angle de 0 à 2pi
    for (int i = 0; i <= nbpoint; i++)
    {
        alpha = i*step; //Evolution de l'angle du grand cercle par rapport au centre
        for(int j=0; j<= nbpoint ; j++)
        {
            beta = j*step;
            double x = (tore.raycercle + tore.raytube * cos(beta))*cos(alpha);
            double y =  (tore.raytube * sin(beta));
            double z = (tore.raycercle + tore.raytube * cos(beta))*sin(alpha);
            normals.push_back(Vector(x, y, z));
            vertices[i*nbpoint+j] = Vector(x,y,z) + tore.vec; // Sommet des coté du disque

            AddTriangle(i*nbpoint+j , i%nbpoint*nbpoint+nbpoint+(j)%nbpoint, i*nbpoint+(j+1)%nbpoint, i*nbpoint+j); //ajout de chaque triangle entre les points
            AddTriangle(i%nbpoint*nbpoint+nbpoint+(j)%nbpoint , i*nbpoint+(j+1)%nbpoint, i%nbpoint*nbpoint+nbpoint+(j+1)%nbpoint,i*nbpoint+j ); //ajout de chaque triangle entre les points

        }
    }
}

/*!
\brief Consctructeur Mesh du HeightField.

\param terrain HeightField.
*/
Mesh::Mesh(const HeightField& terrain)
{
    // Vertices
    vertices.resize(terrain.image.width()*terrain.image.width() + terrain.image.height());

    // Reserve space for the triangle array
    varray.reserve(terrain.image.width()*terrain.image.width() + terrain.image.height());
    narray.reserve(terrain.image.width()*terrain.image.width() + terrain.image.height());

    int dif = terrain.image.width() - terrain.image.height(); //difference entre la longeur et la largeur de l'image

    for(int i=0;i<terrain.image.width()-1;i++) // on parcourt la largeur de l'image
    {
        for(int j=0;j<terrain.image.height()-1;j++) // on parcourt la hauteur de l'image
        {
            int a = terrain.image.pixelColor(i,j).black(); // on recupere l'entier correspondant à la couleur noir du pixel de l'image choisi

            vertices[i*terrain.image.width()+j] = Vector(j,i,-a/terrain.echelle)+terrain.vec; // point par pixel avec une hauteur en fonction de l'entier de sa couleur noir
            normals.push_back(Vector(j, i, 1)+terrain.vec); 

            //ajout de chaque triangle entre les points
            if((i != terrain.image.width()-2) )
            {
                AddTriangle(i*terrain.image.width()+j, (i+1)*terrain.image.width()+j+1, i*terrain.image.width()+j+1,    i*(terrain.image.width()-1-dif)+j);
            }
            if( j != terrain.image.height()-2)
            {
                AddTriangle(i*terrain.image.width()+j, (i+1)*terrain.image.width()+j+1, (i+1)*terrain.image.width()+j , i*(terrain.image.width()-1-dif)+j);
            }

        }
    }
}


/*!
\brief Rotation du Mesh dans l'axe X.

\param angle double, angle de rotation.
*/
void Mesh::RotationX(double angle)
{
    double angle2 = Math::DegreeToRadian(angle);
    Matrix rotX;
    rotX[0]= 1;
    rotX[1]= 0;
    rotX[2]= 0;

    rotX[3]= 0;
    rotX[4]= cos(angle2);
    rotX[5]= -sin(angle2);

    rotX[6]= 0;
    rotX[7]= sin(angle2);
    rotX[8]= cos(angle2);

    for(int i=0; i < vertices.size();i++)
    {
        double a = (rotX[0] * vertices[i][0]) + (rotX[3] * vertices[i][1]) + (rotX[6] * vertices[i][2]);
        double b= (rotX[1] * vertices[i][0]) + (rotX[4] * vertices[i][1]) + (rotX[7] * vertices[i][2]);
        double c= (rotX[2] * vertices[i][0]) + (rotX[5] * vertices[i][1]) + (rotX[8] * vertices[i][2]);
        vertices[i] = Vector(a,b,c);
    }
}

/*!
\brief Rotation du Mesh dans l'axe Y.

\param angle double, angle de rotation.
*/
void Mesh::RotationY(double angle)
{
    double angle2 = Math::DegreeToRadian(angle);
    Matrix rotY;
    rotY[0]= cos(angle2);
    rotY[1]= 0;
    rotY[2]= sin(angle2);
    rotY[3]= 0;
    rotY[4]= 1;
    rotY[5]= 0;
    rotY[6]= -sin(angle2);
    rotY[7]= 0;
    rotY[8]= cos(angle2);
    for(int i=0; i < vertices.size();i++)
    {
        double a= (rotY[0] * vertices[i][0]) + (rotY[3] * vertices[i][1]) + (rotY[6] * vertices[i][2]);
        double b= (rotY[1] * vertices[i][0]) + (rotY[4] * vertices[i][1]) + (rotY[7] * vertices[i][2]);
        double c= (rotY[2] * vertices[i][0]) + (rotY[5] * vertices[i][1]) + (rotY[8] * vertices[i][2]);
        vertices[i] = Vector(a,b,c);
    }
}

/*!
\brief Rotation du Mesh dans l'axe Z.

\param angle double, angle de rotation.
*/
void Mesh::RotationZ(double angle)
{
    Matrix rotZ;
    double angle2 = Math::DegreeToRadian(angle);
    rotZ[0]= cos(angle2);
    rotZ[1]= -sin(angle2);
    rotZ[2]= 0;
    rotZ[3]= sin(angle2);
    rotZ[4]= cos(angle2);
    rotZ[5]= 0;
    rotZ[6]= 0;
    rotZ[7]= 0;
    rotZ[8]= 1;


    for(int i=0; i < vertices.size();i++)
    {
        double a = (rotZ[0] * vertices[i][0]) + (rotZ[3] * vertices[i][1]) + (rotZ[6] * vertices[i][2]);
        double b = (rotZ[1] * vertices[i][0]) + (rotZ[4] * vertices[i][1]) + (rotZ[7] * vertices[i][2]);
        double c = (rotZ[2] * vertices[i][0]) + (rotZ[5] * vertices[i][1]) + (rotZ[8] * vertices[i][2]);
        vertices[i] = Vector(a,b,c);
    }

}

/*!
\brief Homothetie du Mesh.

\param v Vector vecteur d'agrandissement.
*/
void Mesh::Homothetie(Vector v)
{
    Matrix h;
    h[0]= v[0];
    h[1]= 0;
    h[2]= 0;
    h[3]= 0;
    h[4]= v[1];
    h[5]= 0;
    h[6]= 0;
    h[7]= 0;
    h[8]= v[2];

    for(int i=0; i < vertices.size();i++)
    {
        double a = vertices[i][0]= (h[0] * vertices[i][0]) + (h[3] * vertices[i][1]) + (h[6] * vertices[i][2]);
        double b = vertices[i][1]= (h[1] * vertices[i][0]) + (h[4] * vertices[i][1]) + (h[7] * vertices[i][2]);
        double c = vertices[i][2]= (h[2] * vertices[i][0]) + (h[5] * vertices[i][1]) + (h[8] * vertices[i][2]);
        vertices[i] = Vector(a,b,c);
    }
}

/*!
\brief Applatissement d'un mesh general.
\param t double factor.
*/
void Mesh::Applatissement(double t)
{
    for(int i=0; i<vertices.size();i++)
    {
        if(vertices[i][2]> t)
        {
            vertices[i][2]= t;
        }
    }
}

/*!
\brief Translation du Mesh.

\param v Vector vecteur de translation.
*/
void Mesh::Translation(Vector v)
{

    for(int i=0; i < vertices.size();i++)
    {
        vertices[i] = vertices[i] + v;
    }
}


/*!
\brief Scale the mesh.
\param s Scaling factor.
*/
void Mesh::Scale(double s)
{
    // Vertexes
    for (int i = 0; i < vertices.size(); i++)
    {
        vertices[i] *= s;
    }

    if (s < 0.0)
    {
        // Normals
        for (int i = 0; i < normals.size(); i++)
        {
            normals[i] = -normals[i];
        }
    }
}



#include <QtCore/QFile>
#include <QtCore/QTextStream>
#include <QtCore/QRegularExpression>
#include <QtCore/qstring.h>

/*!
\brief Import a mesh from an .obj file.
\param filename File name.
*/
void Mesh::Load(const QString& filename)
{
  vertices.clear();
  normals.clear();
  varray.clear();
  narray.clear();

  QFile data(filename);

  if (!data.open(QFile::ReadOnly))
    return;
  QTextStream in(&data);

  // Set of regular expressions : Vertex, Normal, Triangle
  QRegularExpression rexv("v\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)");
  QRegularExpression rexn("vn\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)");
  QRegularExpression rext("f\\s*(\\d*)/\\d*/(\\d*)\\s*(\\d*)/\\d*/(\\d*)\\s*(\\d*)/\\d*/(\\d*)");
  while (!in.atEnd())
  {
    QString line = in.readLine();
    QRegularExpressionMatch match = rexv.match(line);
    QRegularExpressionMatch matchN = rexn.match(line);
    QRegularExpressionMatch matchT = rext.match(line);
    if (match.hasMatch())//rexv.indexIn(line, 0) > -1)
    {
      Vector q = Vector(match.captured(1).toDouble(), match.captured(2).toDouble(), match.captured(3).toDouble()); vertices.push_back(q);
    }
    else if (matchN.hasMatch())//rexn.indexIn(line, 0) > -1)
    {
      Vector q = Vector(matchN.captured(1).toDouble(), matchN.captured(2).toDouble(), matchN.captured(3).toDouble());  normals.push_back(q);
    }
    else if (matchT.hasMatch())//rext.indexIn(line, 0) > -1)
    {
      varray.push_back(matchT.captured(1).toInt() - 1);
      varray.push_back(matchT.captured(3).toInt() - 1);
      varray.push_back(matchT.captured(5).toInt() - 1);
      narray.push_back(matchT.captured(2).toInt() - 1);
      narray.push_back(matchT.captured(4).toInt() - 1);
      narray.push_back(matchT.captured(6).toInt() - 1);
    }
  }
  data.close();
}

/*!
\brief Save the mesh in .obj format, with vertices and normals.
\param url Filename.
\param meshName %Mesh name in .obj file.
*/
void Mesh::SaveObj(const QString& url, const QString& meshName) const
{
  QFile data(url);
  if (!data.open(QFile::WriteOnly))
    return;
  QTextStream out(&data);
  out << "g " << meshName << Qt::endl;
  for (int i = 0; i < vertices.size(); i++)
    out << "v " << vertices.at(i)[0] << " " << vertices.at(i)[1] << " " << vertices.at(i)[2] << QString('\n');
  for (int i = 0; i < normals.size(); i++)
    out << "vn " << normals.at(i)[0] << " " << normals.at(i)[1] << " " << normals.at(i)[2] << QString('\n');
  for (int i = 0; i < varray.size(); i += 3)
  {
    out << "f " << varray.at(i) + 1 << "//" << narray.at(i) + 1 << " "
      << varray.at(i + 1) + 1 << "//" << narray.at(i + 1) + 1 << " "
      << varray.at(i + 2) + 1 << "//" << narray.at(i + 2) + 1 << " "
      << "\n";
  }
  out.flush();
  data.close();
}

