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

#ifndef HELPER_H
#define HELPER_H

#include <QBrush>
#include <QFont>
#include <QPen>
#include <QVector>
#include <QVector3D>
#include "mesh.h"
//#include <vector>

QT_BEGIN_NAMESPACE
class QPainter;
class QPaintEvent;
QT_END_NAMESPACE

//using namespace std;

//struct pt2
//{
//    float x;
//    float y;
//    pt2(float _x, float _y) {
//        x = _x;
//        y = _y;
//    }
//};

//! [0]
class Helper
{
public:
    Helper();

public:
    void paint(QPainter *painter, QPaintEvent *event, int elapsed);
    void AddCtrlPt(QPoint pt);
    void SetSplineType(int i);
    void FindAndSetPoint(QPoint p_find, QPoint p_set);
    void FindAndDuplicate(QPoint p);
   // Mesh Extrusion(QVector<QPoint> curve, int ns, QVector3D v, bool close);
    void GenerateExtrusion(int ns, QVector3D v);
//	void GenerateRevolution(int ns);
	void ClearPoints();
	void RecordWire();
    void GenerateSweep();

    void GenerateRevolution(int ns, int nu, int nv);
    void LoadMesh();

	void GenerateDooSabin();
	void GenerateCatmullClark();
	void GenerateLoop();
private:
    QBrush background;
    QBrush circleBrush;
    QFont textFont;
    QPen circlePen;
    QPen textPen;
    int _splineType;
    Mesh _inputMesh;

    QVector<QPoint> ctrlPts;
	QVector<QPointF> GetCurve();
	QVector<QPointF> _generator;

};



inline bool farenough(QPoint p0, QPoint p1)
{
    QPoint d = p0 - p1;
    return (d.x()*d.x() + d.y()*d.y() > 16)?true:false;

}
//! [0]

#endif
