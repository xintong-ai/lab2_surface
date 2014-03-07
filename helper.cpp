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

#include <QtGui>
#include "helper.h"

//! [0]
Helper::Helper()
{
    QLinearGradient gradient(QPointF(50, -20), QPointF(80, 20));
    gradient.setColorAt(0.0, Qt::white);
    gradient.setColorAt(1.0, QColor(0xa6, 0xce, 0x39));

    background = QBrush(QColor(255, 255, 255));
    circleBrush = QBrush(gradient);
    circlePen = QPen(Qt::black);
    circlePen.setWidth(1);
    textPen = QPen(Qt::white);
    textFont.setPixelSize(50);
    _splineType = 0;
 //   ctrlPts = new vector<QPoint>();
}

unsigned int Combinations(unsigned int n, unsigned int k)
{
     if (k > n)
         return 0;
     unsigned int r = 1;
     for (unsigned int d = 1; d <= k; ++d)
     {
         r *= n--;
         r /= d;
     }
     return r;
}

//! [0]
inline float Bernstein(float u, int i, int n)   {
    return Combinations(n, i) * pow(u, i) * pow(1.0f - u, n - i);
}

inline QVector<QPoint> Bezier(QVector<QPoint> cp, int numSample)
{
    QVector<QPoint> ret;
    int n = cp.size() - 1;
    for(int iu = 0; iu < numSample; iu++)    {
        float u = 1.0 * iu / (float)(numSample - 1);
        QPointF p(0, 0);
        for(int i = 0; i < (n + 1); i++) {
            p += (Bernstein(u, i, n) * (QPointF)cp.at(i));
        }
        ret.push_back(p.toPoint());
    }
    return ret;
}

inline QVector<QPointF> BezierF(QVector<QPoint> cp, int numSample)
{
    QVector<QPointF> ret;
    int n = cp.size() - 1;
    for(int iu = 0; iu < numSample; iu++)    {
        float u = 1.0 * iu / (float)(numSample - 1);
        QPointF p(0, 0);
        for(int i = 0; i < (n + 1); i++) {
            p += (Bernstein(u, i, n) * (QPointF)cp.at(i));
        }
        ret.push_back(p);
    }
    return ret;
}

inline QVector<QPoint> RationalBezier(QVector<QPoint> cp, int numSample)
{
    QVector<QPoint> ret;
    int n = cp.size() - 1;
    int n_mid = n * 0.5;
    float w = 0;
    for(int iu = 0; iu < numSample; iu++)    {
        float u = 1.0 * iu / (float)(numSample - 1);
        QPointF p(0, 0);
        float denom = 0;
        for(int i = 0; i < (n + 1); i++) {
            w = (i == n_mid)? 5:1;
            p += (Bernstein(u, i, n) * (QPointF)cp.at(i)) * w;
            denom += Bernstein(u, i, n) * w;
        }
        p /= denom;
        ret.push_back(p.toPoint());
    }
    return ret;
}

inline QVector<QPointF> RationalBezierF(QVector<QPoint> cp, int numSample)
{
    QVector<QPointF> ret;
    int n = cp.size() - 1;
    int n_mid = n * 0.5;
    float w = 0;
    for(int iu = 0; iu < numSample; iu++)    {
        float u = 1.0 * iu / (float)(numSample - 1);
        QPointF p(0, 0);
        float denom = 0;
        for(int i = 0; i < (n + 1); i++) {
            w = (i == n_mid)? 5:1;
            p += (Bernstein(u, i, n) * (QPointF)cp.at(i)) * w;
            denom += Bernstein(u, i, n) * w;
        }
        p /= denom;
        ret.push_back(p);
    }
    return ret;
}

inline QVector<QPoint> ClosedBezier(QVector<QPoint> cp, int numSample)
{
    QVector<QPoint> ret;
	if(cp.size() < 3)
		return ret;
    QPointF dir = (cp[1] - cp.back()) * 0.2;
    QPointF newEnd = cp.front() - dir;
    QPointF newFront = cp.front() + dir;
    cp.insert(1, newFront.toPoint());
    cp.push_back(newEnd.toPoint());
    cp.push_back(cp.front());
    int n = cp.size() - 1;
    for(int iu = 0; iu < numSample; iu++)    {
        float u = 1.0 * iu / (float)(numSample - 1);
        QPointF p(0, 0);
        for(int i = 0; i < (n + 1); i++) {
            p += (Bernstein(u, i, n) * (QPointF)cp.at(i));
        }
        ret.push_back(p.toPoint());
    }
    return ret;
}

inline QVector<QPointF> ClosedBezierF(QVector<QPoint> cp, int numSample)
{
    QVector<QPointF> ret;
	if(cp.size() < 3)
		return ret;
    QPointF dir = (cp[1] - cp.back()) * 0.2;
    QPointF newEnd = cp.front() - dir;
    QPointF newFront = cp.front() + dir;
    cp.insert(1, newFront.toPoint());
    cp.push_back(newEnd.toPoint());
    cp.push_back(cp.front());
    int n = cp.size() - 1;
    for(int iu = 0; iu < numSample; iu++)    {
        float u = 1.0 * iu / (float)(numSample - 1);
        QPointF p(0, 0);
        for(int i = 0; i < (n + 1); i++) {
            p += (Bernstein(u, i, n) * (QPointF)cp.at(i));
        }
        ret.push_back(p);
    }
    return ret;
}


inline QVector<QPointF> BezierOneSubdivide(QVector<QPointF> cp, QVector<QPointF> poly1, QVector<QPointF> poly2, float u)
{
    if(cp.size() == 1)  {
        poly1.push_back(cp.front());
        return poly1 + poly2;
    }   else    {
        poly1.push_back(cp.front());

        QVector<QPointF> tmp = poly2;
        poly2.clear();
        poly2.push_back(cp.back());
        poly2 += tmp;

        QVector<QPointF> cp2;
        for(int i = 0; i < cp.size() - 1; i++)  {
            cp2.push_back(cp[i] + u * (cp[i + 1] - cp[i]));
        }
        return BezierOneSubdivide(cp2, poly1, poly2, u);
    }
}

inline QVector<QPointF> BezierSubdivideF(QVector<QPointF> cp, int m, float u)  {
    QVector<QPointF> leftSet, rightSet;
    if(cp.size() < 3)
        return leftSet;
    if(m==1)    {
        return BezierOneSubdivide(cp, leftSet, rightSet, u);
    }
    else    {
        QVector<QPointF> tmp = BezierOneSubdivide(cp, leftSet, rightSet, u);

//        leftSet.clear();
//        rightSet.clear();
        int mid = (tmp.size() - 1) / 2;
        for(int i = 0; i <= mid; i++)
            leftSet.push_back(tmp[i]);
        for(int i = mid; i < tmp.size(); i++)
            rightSet.push_back(tmp[i]);

        leftSet = BezierSubdivideF(leftSet, m - 1, u);
        rightSet = BezierSubdivideF(rightSet, m - 1, u);
        return leftSet + rightSet;
    }
}

//inline float N_i_1(float u, QVector<int> t, int i)  {
//    return (t[i] <= u && u < t[i+1])? 1:0;
//}

inline float N(float u, int i, int d, QVector<int> t)   {
    float ret;
    if(1 == d)  {
        ret = (t[i] <= u && u < t[i+1])? 1:0;
        if(t[i] == t[i+1] && t[i] == u)
            ret = 1;
    }
    else    {
        ret =   (u - t[i]) * N(u, i, d - 1, t) / (t[i + d - 1] - t[i])
              + (t[i + d] - u) * N(u, i + 1, d - 1, t) / (t[i + d] - t[i + 1]);
    }
    return ret;
}

inline QVector<int> FixedEndKnot(int n, int D)
{
    QVector<int> t;
    for(int j = 0; j <= (n + D); j++)  {
        if(j < D)
            t.push_back(0);
        else if(D <= j && j <= n)
            t.push_back(j - D + 1);
        else
            t.push_back(n - D + 2);
    }
    return t;
}

inline QVector<int> UniformKnot(int n, int D)
{
    QVector<int> t;
    for(int j = 0; j < (n + D + 1); j++)  {
        t.push_back(j);
    }
    return t;
}

inline QVector<QPoint> BSpline(QVector<QPoint> cp, int D, int numSample)
{

    QVector<QPoint> ret;
    if(cp.size() < D)
        return ret;
    int n = cp.size() - 1;

    QVector<int> t = UniformKnot(n, D);//FixedEndKnot(n, D);


    for(int iu = 0; iu < numSample; iu++)    {
        float u = t.back() * iu / (float)(numSample - 1);
        QPointF p(0, 0);
        for(int i = 0; i <= n; i++) {
            p += (N(u, i, D, t) * (QPointF)cp.at(i));
        }
        ret.push_back(p.toPoint());
    }
    return ret;
}

inline QPointF CubicMatrix(float u, QPointF p0, QPointF p1, QPointF p2, QPointF p3)
{
    QPointF ret =(    pow(u, 3)                                           * p3
                   + (pow(u, 3) * (-3) + pow(u, 2) * 3    + u * 3    + 1) * p2
                   + (pow(u, 3) * 3    + pow(u, 2) * (-6)            + 4) * p1
                   + (pow(u, 3) * (-1) + pow(u, 2) * 3    + u * (-3) + 1) * p0
                     ) / 6.0;
    return ret;
}

inline QVector<QPoint> CubicBSpline(QVector<QPoint> cp, int numSample)
{
    QVector<QPoint> ret;
    int D = 4;
    if(cp.size() < D)
        return ret;
    //int n = cp.size() - 1;

   // QVector<int> t = UniformKnot(n, D);//FixedEndKnot(n, D);

    for(int i = 0; i <= (cp.size() - D); i++)   {
        QPointF p(0, 0);
        for(int iu = 0; iu < numSample; iu++)    {
            float u = (float)iu / (numSample - 1);
            ret.push_back(CubicMatrix(u, cp[i], cp[i + 1], cp[i + 2], cp[i + 3]).toPoint());
        }
    }
    return ret;
}

inline QVector<QPointF> CubicBSplineF(QVector<QPoint> cp, int numSample)
{
    QVector<QPointF> ret;
    int D = 4;
    if(cp.size() < D)
        return ret;
    //int n = cp.size() - 1;

   // QVector<int> t = UniformKnot(n, D);//FixedEndKnot(n, D);

    for(int i = 0; i <= (cp.size() - D); i++)   {
        QPointF p(0, 0);
        for(int iu = 0; iu < numSample; iu++)    {
            float u = (float)iu / (numSample - 1);
            ret.push_back(CubicMatrix(u, cp[i], cp[i + 1], cp[i + 2], cp[i + 3]));
        }
    }
    return ret;
}

inline Mesh Extrusion(QVector<QPointF> curve, int ns, QVector3D v, bool close)
{
    Mesh mMesh;
    QVector<QVector3D> c0, c1;
     //QVector<QPoint> curve = CubicBSpline(ctrlPts, 16);
    for(int i = 0; i < curve.size(); i++) {
        c0.push_back(QVector3D(curve[i].x(), curve[i].y(), 0));
    }
    QVector3D vstep = v / (ns - 1);
    for(int i = 1; i < ns; i++)   {
        c1.clear();
        for(int j = 0; j < curve.size(); j++)
            c1.push_back(QVector3D(curve[j].x(), curve[j].y(), 0) + vstep * i);
        for(int j = 0; j < (c0.size() - 1); j++) {
            mMesh.AddFacet(c0[j].x(), c0[j].y(), c0[j].z(),
                           c0[j + 1].x(), c0[j + 1].y(), c0[j + 1].z(),
                           c1[j].x(), c1[j].y(), c1[j].z());
            mMesh.AddFacet(c0[j + 1].x(), c0[j + 1].y(), c0[j + 1].z(),
                           c1[j + 1].x(), c1[j + 1].y(), c1[j + 1].z(),
                           c1[j].x(), c1[j].y(), c1[j].z());
        }
        c0 = c1;
    }
    return mMesh;
}

inline Mesh Revolution(QVector<QPointF> curve, int ns, bool close)
{
    Mesh mMesh;
    QVector<QVector3D> c0, c1, cc;
    for(int i = 0; i < curve.size(); i++) {
        c0.push_back(QVector3D(curve[i].x(), curve[i].y(), 0));
    }
	cc = c0;
	float aStep = 2 * M_PI / ns;
    for(int i = 1; i <= ns; i++)   {
        c1.clear();
		float angle = (i % ns) * aStep;
        for(int j = 0; j < curve.size(); j++)
            c1.push_back(QVector3D(cc[j].x(), cc[j].y() * cos(angle), cc[j].y() * sin(angle)));
        for(int j = 0; j < (c0.size() - 1); j++) {
            mMesh.AddFacet(c0[j].x(), c0[j].y(), c0[j].z(),
                           c0[j + 1].x(), c0[j + 1].y(), c0[j + 1].z(),
                           c1[j].x(), c1[j].y(), c1[j].z());
            mMesh.AddFacet(c0[j + 1].x(), c0[j + 1].y(), c0[j + 1].z(),
                           c1[j + 1].x(), c1[j + 1].y(), c1[j + 1].z(),
                           c1[j].x(), c1[j].y(), c1[j].z());
        }
        c0 = c1;
    }
    return mMesh;
}

inline Mesh Sweep(QVector<QPointF> generator, QVector<QPointF> trajectory)
{
    Mesh mMesh;
    QVector<QVector3D> c0, c1;
     //QVector<QPoint> curve = CubicBSpline(ctrlPts, 16);
    for(int i = 0; i < generator.size(); i++) {
        c0.push_back(QVector3D(generator[i].x(), generator[i].y(), 0));
    }
//    QVector3D vstep = v / (ns - 1);
	for(int i = 1; i < trajectory.size(); i++)   {
        c1.clear();
        for(int j = 0; j < generator.size(); j++)
            c1.push_back(
				QVector3D(generator[j].x(), 
				generator[j].y() + trajectory[i].x() - trajectory[0].x(), 
				trajectory[i].y() - trajectory[0].y()));
        for(int j = 0; j < (c0.size() - 1); j++) {
            mMesh.AddFacet(c0[j].x(), c0[j].y(), c0[j].z(),
                           c0[j + 1].x(), c0[j + 1].y(), c0[j + 1].z(),
                           c1[j].x(), c1[j].y(), c1[j].z());
            mMesh.AddFacet(c0[j + 1].x(), c0[j + 1].y(), c0[j + 1].z(),
                           c1[j + 1].x(), c1[j + 1].y(), c1[j + 1].z(),
                           c1[j].x(), c1[j].y(), c1[j].z());
        }
        c0 = c1;
    }
    return mMesh;
}
	
void Helper::GenerateRevolution(int ns)
{
    QVector<QPointF> curveF = GetCurve();
	Mesh mMesh;
    bool closed = (_splineType == 3);
    mMesh = Revolution(curveF, ns, closed);
    mMesh.WritePLY("data/revolution.ply");
}

void Helper::GenerateSweep()
{
	QVector<QPointF> trajectory = GetCurve();
	Mesh mMesh;
    mMesh = Sweep(_generator, trajectory);
    mMesh.WritePLY("data/sweep.ply");
}

void Helper::ClearPoints()
{
	ctrlPts.clear();
}

void Helper::RecordWire()
{
	_generator = GetCurve();
}

QVector<QPointF> Helper::GetCurve()
{
    QVector<QPointF> curveF;
    QVector<QPointF> ctrlPtsF;
    for(int i = 0; i < ctrlPts.size(); i++)
        ctrlPtsF.push_back(ctrlPts[i]);

    switch(_splineType)
    {
    case 0:
        curveF = BezierF(ctrlPts, ctrlPts.size() * 32);
        break;
    case 1:
        curveF = CubicBSplineF(ctrlPts, 16);
        break;
    case 2:
        curveF = BezierSubdivideF(ctrlPtsF, 4, 0.5);
        break;
    case 3:
        curveF = ClosedBezierF(ctrlPts, ctrlPts.size() * 32);
        break;
    case 4:
        curveF = RationalBezierF(ctrlPts, ctrlPts.size() * 32);
        break;
    }
	return curveF;
}

void Helper::GenerateExtrusion(int ns, QVector3D v)
{
    QVector<QPointF> curveF = GetCurve();
    Mesh mMesh;
    bool closed = (_splineType == 3);
    mMesh = Extrusion(curveF, ns, v, closed);
    //mMesh = Extrusion(ctrlPtsF, ns, v, closed);
    mMesh.WritePLY("data/extrusion.ply");

}

//! [1]
void Helper::paint(QPainter *painter, QPaintEvent *event, int elapsed)
{

    painter->fillRect(event->rect(), background);
//    painter->translate(100, 100);
//! [1]

//! [2]
    painter->save();
    painter->setBrush(circleBrush);
    painter->setPen(circlePen);
//    painter->rotate(elapsed * 0.030);

	if(ctrlPts.size() < 1)
		return;
    for(int i = 0; i < ctrlPts.size(); i++)    {
        painter->drawEllipse(ctrlPts.at(i), 4, 4);
    }
//
//    painter->drawPolyline(ctrlPts);
//    painter->drawPolyline(BSpline(ctrlPts, 1, ctrlPts.size() * 32));
//
    QVector<QPointF> ctrlPtsF;
    switch(_splineType)
    {
    case 0:
        painter->drawPolyline(Bezier(ctrlPts, ctrlPts.size() * 32));
        break;
    case 1:
        painter->drawPolyline(CubicBSpline(ctrlPts, 16));
        break;
    case 2:
        for(int i = 0; i < ctrlPts.size(); i++)
            ctrlPtsF.push_back(ctrlPts[i]);
        painter->drawPolyline(BezierSubdivideF(ctrlPtsF, 4, 0.5));
        break;
    case 3:
        painter->drawPolyline(ClosedBezier(ctrlPts, ctrlPts.size() * 32));
        break;
    case 4:
        painter->drawPolyline(RationalBezier(ctrlPts, ctrlPts.size() * 32));
        break;
    }





//    qreal r = elapsed/1000.0;
//    int n = 30;
//    for (int i = 0; i < n; ++i) {
//        painter->rotate(30);
//        qreal radius = 0 + 120.0*((i+r)/n);
//        qreal circleRadius = 1 + ((i+r)/n)*20;
//        painter->drawEllipse(QRectF(radius, -circleRadius,
//                                    circleRadius*2, circleRadius*2));
//    }
//    painter->restore();
////! [2]

////! [3]
//    painter->setPen(textPen);
//    painter->setFont(textFont);
//    painter->drawText(QRect(-50, -50, 100, 100), Qt::AlignCenter, "Qt");
}
//! [3]

void Helper::AddCtrlPt(QPoint pt)
{
    ctrlPts.push_back(pt);
}

void Helper::SetSplineType(int i)
{
    _splineType = i;

}

void Helper::FindAndSetPoint(QPoint p_find, QPoint p_set)
{
    for(int i = 0; i < ctrlPts.size(); i++) {
        if(!farenough(p_find, ctrlPts[i]))   {
            ctrlPts[i] = p_set;
            break;
        }
    }
}

void Helper::FindAndDuplicate(QPoint p)
{
    for(int i = 0; i < ctrlPts.size(); i++) {
        if(!farenough(p, ctrlPts[i]))   {
            ctrlPts.insert(i + 1, ctrlPts[i]);
            break;
        }
    }
}

