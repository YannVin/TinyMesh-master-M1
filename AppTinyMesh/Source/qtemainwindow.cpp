#include "qte.h"
#include "implicits.h"
#include "ui_interface.h"

MainWindow::MainWindow() : QMainWindow(), uiw(new Ui::Assets)
{
	// Chargement de l'interface
    uiw->setupUi(this);

	// Chargement du GLWidget
	meshWidget = new MeshWidget;
	QGridLayout* GLlayout = new QGridLayout;
	GLlayout->addWidget(meshWidget, 0, 0);
	GLlayout->setContentsMargins(0, 0, 0, 0);
    uiw->widget_GL->setLayout(GLlayout);

	// Creation des connect
	CreateActions();

	meshWidget->SetCamera(Camera(Vector(10, 0, 0), Vector(0.0, 0.0, 0.0)));
}

MainWindow::~MainWindow()
{
	delete meshWidget;
}

void MainWindow::CreateActions()
{
	// Buttons
    connect(uiw->boxMesh, SIGNAL(clicked()), this, SLOT(BoxMeshExample()));
    connect(uiw->disqueMeshButton, SIGNAL(clicked()), this, SLOT(DisqueMeshExample()));
    connect(uiw->cylindreMeshButton, SIGNAL(clicked()), this, SLOT(CylindreMeshExample()));
    connect(uiw->sphereMeshButton, SIGNAL(clicked()), this, SLOT(SphereMeshExample()));
    connect(uiw->capsuleMeshButton, SIGNAL(clicked()), this, SLOT(CapsuleMeshExample()));
    connect(uiw->toreMeshButton, SIGNAL(clicked()), this, SLOT(ToreMeshExample()));
    connect(uiw->heightFieldMeshButton, SIGNAL(clicked()), this, SLOT(HeightFieldMeshExample()));
    connect(uiw->complexMeshButton, SIGNAL(clicked()), this, SLOT(Complex()));

    connect(uiw->sphereImplicit, SIGNAL(clicked()), this, SLOT(SphereImplicitExample()));
    connect(uiw->resetcameraButton, SIGNAL(clicked()), this, SLOT(ResetCamera()));
    connect(uiw->wireframe, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
    connect(uiw->radioShadingButton_1, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
    connect(uiw->radioShadingButton_2, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));

	// Widget edition
	connect(meshWidget, SIGNAL(_signalEditSceneLeft(const Ray&)), this, SLOT(editingSceneLeft(const Ray&)));
	connect(meshWidget, SIGNAL(_signalEditSceneRight(const Ray&)), this, SLOT(editingSceneRight(const Ray&)));
}

void MainWindow::editingSceneLeft(const Ray&)
{
}

void MainWindow::editingSceneRight(const Ray&)
{
}

void MainWindow::BoxMeshExample()
{
	Mesh boxMesh = Mesh(Box(1.0));

	std::vector<Color> cols;
	cols.resize(boxMesh.Vertexes());
    for (size_t i = 0; i < cols.size(); i++)
		cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

	meshColor = MeshColor(boxMesh, cols, boxMesh.VertexIndexes());
	UpdateGeometry();
}

void MainWindow::DisqueMeshExample()
{
    Mesh disqueMesh = Mesh(Disque(Vector(0,0,0),1),200);

    std::vector<Color> cols;
    cols.resize(disqueMesh.Vertexes());
    for (size_t i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    meshColor = MeshColor(disqueMesh, cols, disqueMesh.VertexIndexes());
    UpdateGeometry();
}

void MainWindow::CylindreMeshExample()
{
    Mesh cylindreMesh = Mesh(Cylindre(Vector(0,0,0),Vector(0,0,10),5),100);

    std::vector<Color> cols;
    cols.resize(cylindreMesh.Vertexes());
    for (size_t i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    meshColor = MeshColor(cylindreMesh, cols, cylindreMesh.VertexIndexes());
    UpdateGeometry();
}

void MainWindow::SphereMeshExample()
{
    Mesh sphereMesh = Mesh(Sphere(Vector(0,0,0),3,3),100);

    std::vector<Color> cols;
    cols.resize(sphereMesh.Vertexes());
    for (size_t i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    meshColor = MeshColor(sphereMesh, cols, sphereMesh.VertexIndexes());
    UpdateGeometry();
}

void MainWindow::CapsuleMeshExample()
{
    Mesh CapsuleMesh = Mesh(Capsule(Vector(0,0,0),3,1),100);

    std::vector<Color> cols;
    cols.resize(CapsuleMesh.Vertexes());
    for (size_t i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    meshColor = MeshColor(CapsuleMesh, cols, CapsuleMesh.VertexIndexes());
    UpdateGeometry();
}

void MainWindow::ToreMeshExample()
{
    Mesh ToreMesh = Mesh(Tore(),100);

    std::vector<Color> cols;
    cols.resize(ToreMesh.Vertexes());
    for (size_t i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    meshColor = MeshColor(ToreMesh, cols, ToreMesh.VertexIndexes());
    UpdateGeometry();
}

void MainWindow::HeightFieldMeshExample()
{
    int echelle = 5;
    QString s = QCoreApplication::applicationDirPath();
    QImage a = QImage(s+":/../../TinyMesh-master/AppTinyMesh/data/Texture_terrain.png","png");

    Mesh HeightFieldMesh = Mesh(HeightField(a,echelle));

    std::vector<Color> cols;
    cols.resize(HeightFieldMesh.Vertexes());
    for (size_t i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    meshColor = MeshColor(HeightFieldMesh, cols, HeightFieldMesh.VertexIndexes());
    UpdateGeometry();
}

// ensemble de meshes pour avoir d'un bonhomme de neige
void MainWindow::Complex()
{
    int nbpoint = 50;
    Mesh globale = Mesh(Sphere(Vector(0,0,0),5,5),nbpoint);
    Mesh tronc = Mesh(Sphere(Vector(0,0,7),4,4),nbpoint);
    Mesh tete = Mesh(Sphere(Vector(0,0,13),3,3),nbpoint);
    Mesh oeuil = Mesh(Tore(Vector(0,0,20),1,3),nbpoint);
    Mesh oeuil2 = Mesh(Tore(Vector(0,0,20),1,3),nbpoint);
    Mesh chapeau = Mesh(Cylindre(Vector(0,0,30),Vector(0,0,1),8),nbpoint);
    Mesh chapeau2 = Mesh(Cylindre(Vector(0,0,30),Vector(0,0,10),5),nbpoint);
    Mesh bouche = Mesh(Tore(Vector(0,0,20),1,3),nbpoint);
    Mesh bras1 = Mesh(Capsule(Vector(0,0,0),6,0.5),nbpoint);
    Mesh bras2 = Mesh(Capsule(Vector(0,0,0),6,0.5),nbpoint);
    Mesh baton = Mesh(Capsule(Vector(0,0,0),20,0.5),nbpoint);
    Mesh baton2 = Mesh(Capsule(Vector(0,0,0),4,0.5),nbpoint);
    Mesh baton3 = Mesh(Capsule(Vector(0,0,0),4,0.5),nbpoint);

    //oreilles
    tete.Warp(Sphere(Vector(-3,0,13),1,1),Vector(-1,0,0.5));
    tete.Warp(Sphere(Vector(3,0,13),1,1),Vector(1,0,0.5));

    //nez
    tete.Warp(Sphere(Vector(0,-3,12.8),0.7,0.7),Vector(0,-2,0));

    baton3.Homothetie(Vector(0.5,1,1));
    baton3.RotationZ(90);
    baton3.RotationY(135);
    baton3.Translation(Vector(8,0,17));

    baton2.Homothetie(Vector(0.5,1,1));
    baton2.RotationZ(90);
    baton2.RotationY(45);
    baton2.Translation(Vector(8,0,17));

    baton.Homothetie(Vector(1,1,0.5));
    baton.RotationX(90);
    baton.Translation(Vector(8,0,20));

    bras2.RotationZ(90);
    bras2.RotationY(145);
    bras2.Translation(Vector(-3,0,8));

    bras1.RotationZ(90);
    bras1.RotationY(35);

    bras1.Translation(Vector(3,0,8));

    bouche.Homothetie(Vector(0.5,0.3,0.2));
    bouche.Translation(Vector(0,-2.4,8));

    chapeau.Homothetie(Vector(0.5,0.5,0.5));
    chapeau2.Homothetie(Vector(0.5,0.5,0.5));

    oeuil.Homothetie(Vector(0.2,0.2,0.2));
    oeuil.Translation(Vector(1,-2.75,10));

    oeuil2.Homothetie(Vector(0.2,0.2,0.2));
    oeuil2.Translation(Vector(-1,-2.75,10));

    globale.Merge(tronc);
    globale.Merge(tete);
    globale.Merge(oeuil);
    globale.Merge(oeuil2);
    globale.Merge(chapeau);
    globale.Merge(chapeau2);
    globale.Merge(bouche);
    globale.Merge(bras1);
    globale.Merge(bras2);
    globale.Merge(baton);
    globale.Merge(baton2);
    globale.Merge(baton3);

    std::vector<Color> cols;
    cols.resize(globale.Vertexes());
    for (size_t i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    meshColor = MeshColor(globale, cols, globale.VertexIndexes());
    UpdateGeometry();
}



void MainWindow::SphereImplicitExample()
{
  AnalyticScalarField implicit;

  Mesh implicitMesh;
  implicit.Polygonize(31, implicitMesh, Box(2.0));

  std::vector<Color> cols;
  cols.resize(implicitMesh.Vertexes());
  for (size_t i = 0; i < cols.size(); i++)
    cols[i] = Color(0.8, 0.8, 0.8);

  meshColor = MeshColor(implicitMesh, cols, implicitMesh.VertexIndexes());
  UpdateGeometry();
}

void MainWindow::UpdateGeometry()
{
	meshWidget->ClearAll();
	meshWidget->AddMesh("BoxMesh", meshColor);

    uiw->lineEdit->setText(QString::number(meshColor.Vertexes()));
    uiw->lineEdit_2->setText(QString::number(meshColor.Triangles()));

	UpdateMaterial();
}

void MainWindow::UpdateMaterial()
{
    meshWidget->UseWireframeGlobal(uiw->wireframe->isChecked());

    if (uiw->radioShadingButton_1->isChecked())
		meshWidget->SetMaterialGlobal(MeshMaterial::Normal);
	else
		meshWidget->SetMaterialGlobal(MeshMaterial::Color);
}

void MainWindow::ResetCamera()
{
	meshWidget->SetCamera(Camera(Vector(-10.0), Vector(0.0)));
}
