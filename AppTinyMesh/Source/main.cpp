#include "qte.h"
#include <QtWidgets/qapplication.h>
#include <iostream>

int main(int argc, char *argv[])
{
	QApplication app(argc, argv);

	MainWindow mainWin;
	mainWin.showMaximized();

	return app.exec();
}
