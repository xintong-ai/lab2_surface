/****************************************************************************
**
** Copyright (C) 2012 Nokia Corporation and/or its subsidiary(-ies).
** All rights reserved.
** Contact: Nokia Corporation (qt-info@nokia.com)
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** You may use this file under the terms of the BSD license as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of Nokia Corporation and its Subsidiary(-ies) nor
**     the names of its contributors may be used to endorse or promote
**     products derived from this software without specific prior written
**     permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
** $QT_END_LICENSE$
**
****************************************************************************/



#include "window.h"

//! [0]
Window::Window()
    : QWidget()
{
    //Widget *native = new Widget(&helper, this);
    openGL = new GLWidget(&helper, this);
    //QLabel *nativeLabel = new QLabel(tr("Native"));
   // nativeLabel->setAlignment(Qt::AlignHCenter);
    QLabel *openGLLabel = new QLabel(tr("OpenGL"));
    openGLLabel->setAlignment(Qt::AlignHCenter);

    // radio buttons
    QGroupBox *groupBox = new QGroupBox(tr("Spline Types:"));
    radio0 = new QRadioButton(tr("Bezier curve"));
    radio1 = new QRadioButton(tr("Cubic B-spline"));
    radio2 = new QRadioButton(tr("Subdivision curves"));
    radio3 = new QRadioButton(tr("Closed Bezier curve"));
    radio4 = new QRadioButton(tr("Rational Bezier curve"));
	textX = new QLineEdit(tr("0"));
	textY = new QLineEdit(tr("0"));
	textZ = new QLineEdit(tr("600"));
    textNumSlice = new QLineEdit(tr("16"));
    textSampleNumX = new QLineEdit(tr("64"));
    textSampleNumY = new QLineEdit(tr("64"));
    radio0->setChecked(true);
    QPushButton *buttonExtrusion = new QPushButton(tr("Extrusion"));
    QPushButton *buttonRevolution = new QPushButton(tr("Revolution Bezier"));
    QPushButton *buttonRevolution2 = new QPushButton(tr("Revolution B-spline"));
    QPushButton *buttonClearPoints = new QPushButton(tr("Clear Points"));
    QPushButton *buttonRecordWire = new QPushButton(tr("Record Generator"));
    QPushButton *buttonSweep = new QPushButton(tr("Sweep"));
    QPushButton *buttonLoadMesh = new QPushButton(tr("Load Mesh"));
    QPushButton *buttonDooSabin = new QPushButton(tr("Doo Sabin"));
	QPushButton *buttonCatmullClark = new QPushButton(tr("Catmull-Clark"));
	QPushButton *buttonLoop = new QPushButton(tr("Loop"));
	QPushButton *buttonNNCrust = new QPushButton(tr("NNCrust"));
	QPushButton *buttonCrust = new QPushButton(tr("Crust"));
    QVBoxLayout *vbox = new QVBoxLayout;
    vbox->addWidget(radio0);
    vbox->addWidget(radio1);
    vbox->addWidget(radio2);
    vbox->addWidget(radio3);
    vbox->addWidget(radio4);
	vbox->addWidget(textX);
	vbox->addWidget(textY);
	vbox->addWidget(textZ);
	vbox->addWidget(textNumSlice);
    vbox->addWidget(textSampleNumX);
    vbox->addWidget(textSampleNumY);

    vbox->addWidget(buttonExtrusion);
	vbox->addWidget(buttonRevolution);
	vbox->addWidget(buttonRevolution2);
	vbox->addWidget(buttonClearPoints);
	vbox->addWidget(buttonRecordWire);
	vbox->addWidget(buttonSweep);
    vbox->addWidget(buttonLoadMesh);
    vbox->addWidget(buttonDooSabin);
    vbox->addWidget(buttonCatmullClark);
    vbox->addWidget(buttonLoop);
    vbox->addWidget(buttonNNCrust);
    vbox->addWidget(buttonCrust);

    vbox->addStretch(1);
    groupBox->setLayout(vbox);




    QGridLayout *layout = new QGridLayout;
    //layout->addWidget(native, 0, 0);
    layout->addWidget(openGL, 0, 0);
    //layout->addWidget(nativeLabel, 1, 0);
    layout->addWidget(openGLLabel, 1, 0);
    layout->addWidget(groupBox, 0,1);
    setLayout(layout);



    connect(radio0, SIGNAL(clicked()), this, SLOT(SplineTypeSelected()));
    connect(radio1, SIGNAL(clicked()), this, SLOT(SplineTypeSelected()));
    connect(radio2, SIGNAL(clicked()), this, SLOT(SplineTypeSelected()));
    connect(radio3, SIGNAL(clicked()), this, SLOT(SplineTypeSelected()));
    connect(radio4, SIGNAL(clicked()), this, SLOT(SplineTypeSelected()));
    connect(buttonExtrusion, SIGNAL(clicked()), this, SLOT(GenerateExtrusionSurface()));
    connect(buttonRevolution, SIGNAL(clicked()), this, SLOT(GenerateRevolutionSurface()));
    connect(buttonRevolution2, SIGNAL(clicked()), this, SLOT(GenerateRevolutionSurface2()));
    connect(buttonClearPoints, SIGNAL(clicked()), this, SLOT(ClearPoints()));
    connect(buttonRecordWire, SIGNAL(clicked()), this, SLOT(RecordWire()));
    connect(buttonSweep, SIGNAL(clicked()), this, SLOT(GenerateSweep()));
    connect(buttonLoadMesh, SIGNAL(clicked()), this, SLOT(LoadMesh()));
    connect(buttonDooSabin, SIGNAL(clicked()), this, SLOT(DooSabin()));
    connect(buttonCatmullClark, SIGNAL(clicked()), this, SLOT(CatmullClark()));
    connect(buttonLoop, SIGNAL(clicked()), this, SLOT(Loop()));
    connect(buttonNNCrust, SIGNAL(clicked()), this, SLOT(NNCrust()));
    connect(buttonCrust, SIGNAL(clicked()), this, SLOT(Crust()));

 //   connect(timer, SIGNAL(timeout()), openGL, SLOT(animate()));


    setWindowTitle(tr("Lab2: Surface Generation"));
}

void Window::SplineTypeSelected()
{
    int type = -1;
    if(radio0->isChecked())
        type = 0;
    else if(radio1->isChecked())
        type = 1;
    else if(radio2->isChecked())
        type = 2;
    else if(radio3->isChecked())
        type = 3;
    else if(radio4->isChecked())
        type = 4;
    helper.SetSplineType(type);
    openGL->repaint();
}

void Window::GenerateExtrusionSurface()
{
	helper.GenerateExtrusion(textNumSlice->text().toInt(), 
		QVector3D(textX->text().toFloat(), textY->text().toFloat(), textZ->text().toFloat()));
}


void Window::GenerateRevolutionSurface()
{
	int option = 0;
	//if(radio0->isChecked())
	//	option = 0;
	//else //if(radio1->isChecked())
	//	option = 1;
    helper.GenerateRevolution(textNumSlice->text().toInt(), textSampleNumX->text().toInt(), textSampleNumY->text().toInt(), option);
}

void Window::GenerateRevolutionSurface2()
{
	int option = 1;
	//if(radio0->isChecked())
	//	option = 0;
	//else //if(radio1->isChecked())
	//	option = 1;
    helper.GenerateRevolution(textNumSlice->text().toInt(), textSampleNumX->text().toInt(), textSampleNumY->text().toInt(), option);
}

void Window::GenerateSweep()
{
	helper.GenerateSweep();
}

void Window::LoadMesh()
{
	QString	filename = QFileDialog::getOpenFileName(this,
			 tr("Open OFF File"), "data/", tr("OFF Files (*.off)"));
	//char* filename = setText(file1Name);

	helper.LoadMesh(filename.toStdString().c_str());
}

void Window::DooSabin()
{
	helper.GenerateDooSabin(textNumSlice->text().toInt());
}

void Window::CatmullClark()
{
	helper.GenerateCatmullClark(textNumSlice->text().toInt());
}

void Window::Loop()
{
	helper.GenerateLoop(textNumSlice->text().toInt());
}

void Window::NNCrust()
{
	helper.GenNNCrust();
    openGL->repaint();
}

void Window::Crust()
{
	helper.GenCrust();
    openGL->repaint();
}

void Window::ClearPoints()
{
	helper.ClearPoints();
	openGL->repaint();
}

void Window::RecordWire()
{
	helper.RecordWire();
}
//! [0]
