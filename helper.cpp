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
#include <CGAL/Vector_24.h>

struct quad
{
	QVector3D v[4];
};

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

    
inline Mesh Grid2Mesh(QVector<QVector<QVector3D> > grid)
{
    Mesh mMesh;
	for(int i = 0; i < (grid.size() - 1); i ++)
		for(int j = 0; j < (grid.at(0).size() - 1); j++)
		{
			vector<GeomVert> facet;
			facet.push_back(GeomVert(grid[i][j].x(),	grid[i][j].y(),		grid[i][j].z()));
			facet.push_back(GeomVert(grid[i+1][j].x(),	grid[i+1][j].y(),	grid[i+1][j].z()));
			facet.push_back(GeomVert(grid[i+1][j+1].x(),grid[i+1][j+1].y(),	grid[i+1][j+1].z()));
			facet.push_back(GeomVert(grid[i][j+1].x(),	grid[i][j+1].y(),	grid[i][j+1].z()));
			mMesh.AddFacet(facet);
		}
    return mMesh;
}

inline void Grid2Mesh(Mesh &mMesh, QVector<QVector<QVector3D> > grid)
{
    
	for(int i = 0; i < (grid.size() - 1); i ++)
		for(int j = 0; j < (grid.at(0).size() - 1); j++)
		{
         //   mMesh.AddFacet(grid[i][j].x(),		grid[i][j].y(),		grid[i][j].z(),
         //                  grid[i+1][j].x(),	grid[i+1][j].y(),	grid[i+1][j].z(),
         //                  grid[i][j+1].x(),	grid[i][j+1].y(),	grid[i][j+1].z() );
         //   mMesh.AddFacet(grid[i][j+1].x(),	grid[i][j+1].y(),	grid[i][j+1].z(),
						   //grid[i+1][j].x(),	grid[i+1][j].y(),	grid[i+1][j].z(),
						   //grid[i+1][j+1].x(),	grid[i+1][j+1].y(),	grid[i+1][j+1].z() );
			vector<GeomVert> facet;
			facet.push_back(GeomVert(grid[i][j].x(),	grid[i][j].y(),		grid[i][j].z()));
			facet.push_back(GeomVert(grid[i+1][j].x(),	grid[i+1][j].y(),	grid[i+1][j].z()));
			facet.push_back(GeomVert(grid[i+1][j+1].x(),grid[i+1][j+1].y(),	grid[i+1][j+1].z()));
			facet.push_back(GeomVert(grid[i][j+1].x(),	grid[i][j+1].y(),	grid[i][j+1].z()));
			mMesh.AddFacet(facet);
		}
}

QVector3D CasteljauSurf(QVector<QVector<QVector3D> > cp, float w, float u)
{
    int m = cp.size() - 1;
    int n = cp.at(0).size() - 1;
    QVector<QVector3D> p_i;
	for(int i = 0; i <= m; i++)	{
        QVector<QVector3D> p0;
        QVector<QVector3D> p;
        p0 = cp[i];
		for(int j = 1; j <=n; j++)	{
			for(int k = 0; k <= (n - j); k++)	{
                p.push_back((1 - w) * p0[k] + w * p0[k+1]);
			}
            p0 = p;
            p.clear();
		}
        p_i.push_back(p0[0]);
	}
    QVector<QVector3D> p0 = p_i;
    QVector<QVector3D> p;
    for(int j = 1; j <= m; j++)  {
        for(int i = 0; i <= (m - j); i++)    {
            p.push_back((1 - u) * p0[i] + u * p0[i+1]);
        }
        p0 = p;
        p.clear();
    }
    return p0[0];
}

inline QVector<QVector<QVector3D> > ControlGrid2BeizerGrid(QVector<QVector<QVector3D> > cp, int ni, int nj)
{
    QVector<QVector<QVector3D> > grid;
	int m = cp.size() - 1;
	int n = cp.at(0).size() - 1;
    float stepi = 1.0 / ni;
    float stepj = 1.0 / nj;

    for(int i = 0; i <= ni; i++)   {
        QVector<QVector3D>  line;
        for(int j = 0; j <= nj; j++)    {
            line.push_back(CasteljauSurf(cp, i * stepi, j * stepj));
        }
        grid.push_back(line);
    }

    return grid;
}

inline QVector<QVector<QVector3D> > Patch2Grid(QVector<QVector3D> patch, int nx, int ny)
{
	//QVector<quad> quads;
	float stepx = 1.0 / nx;
	float stepy = 1.0 / ny;
	QVector<QVector<QVector3D> > grid;
	for(int i = 0; i <= nx; i++)
	{
		QVector<QVector3D> line;
		float u = stepx * i;
		for(int j = 0; j <= ny; j++)
		{
			float w = stepy * j;
			QVector4D uu(u*u*u, u*u, u, 1);
			QVector4D ww(w*w*w, w*w, w, 1);
			QMatrix4x4 M(	-1.0/6 ,	1.0/2,	-1.0/2,	1.0/6,
							1.0/2,		-1,		1.0/2,	0,
							-1.0/2,		0,		1.0/2,	0,
							1.0/6,		2.0/3,	1.0/6,	0);
			QVector4D left = uu * M;
			QVector4D right = M.transposed() * ww;
			QVector<QVector3D> tmp;
			tmp.push_back(left.x() * patch[0] + left.y() * patch[4] + left.z() * patch[8] + left.w() * patch[12]);
			tmp.push_back(left.x() * patch[1] + left.y() * patch[5] + left.z() * patch[9] + left.w() * patch[13]);
			tmp.push_back(left.x() * patch[2] + left.y() * patch[6] + left.z() * patch[10] + left.w() * patch[14]);
			tmp.push_back(left.x() * patch[3] + left.y() * patch[7] + left.z() * patch[11] + left.w() * patch[15]);

			line.push_back(tmp[0] * right.x() + tmp[1] * right.y() + tmp[2] * right.z() + tmp[3] * right.w());
		}
		grid.push_back(line);
	}
	return grid;
}

inline Mesh ControlGrid2CubicBSplineMesh(QVector<QVector<QVector3D> > cp, int ni, int nj)
{
	Mesh mMesh; 
	bool close_x = (cp.front().front() == cp.back().front());
	bool close_y = (cp.front().front() == cp.front().back());
	if(close_x )
		cp.pop_back();
	
	if(close_y)
	{
		for(int i = 0; i < cp.size(); i++)	{
			cp[i].pop_back();
		}
	}
	
	int nx = cp.size();
	int ny = cp.at(0).size();


	int nx_add = close_x? 3:0;
	int ny_add = close_y? 3:0;
	for(int i = 0; i < (nx + nx_add - 3); i++)
		for(int j = 0; j < (ny + ny_add - 3); j++)
		{
			QVector<QVector3D> patch;
			for(int k = 0; k < 4; k++)
				for(int l = 0; l < 4; l++)
					patch.push_back(cp[(i + k) % nx][(j + l) % ny]);
			Grid2Mesh(mMesh, Patch2Grid(patch, ni, nj));
		}
	return mMesh;
}

inline Mesh RevolutionBezier(QVector<QPointF> curve, int ns, int nu, int nv, int option)
{
    QVector<QVector3D> c1, cc;
    QVector<QVector<QVector3D> > cp;
    for(int i = 0; i < curve.size(); i++) {
        cc.push_back(QVector3D(curve[i].x(), curve[i].y(), 0));
    }
	float aStep = 2 * M_PI / ns;
    for(int i = 0; i <= ns; i++)   {
		float angle = (i % ns) * aStep;
        for(int j = 0; j < curve.size(); j++)
        {
            QVector3D v;
            if(i == (ns - 1))
                v = 2 * cp[0][j] - cp[1][j];
            else
                v = QVector3D(cc[j].x(), cc[j].y() * cos(angle), cc[j].y() * sin(angle));
            c1.push_back(v);
        }
		cp.push_back(c1);
        c1.clear();
    }
	if(0 == option)
		return Grid2Mesh(ControlGrid2BeizerGrid(cp, nu, nv));
	else
	    return ControlGrid2CubicBSplineMesh(cp, nu, nv);
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

inline GeomVert GetFaceVertex(Mesh mMesh, GeomVert v, TopoFacet f)
{
	int j0;
	int nv = f.GetNumberVertices();
	GeomVert x(0, 0, 0);
	for(int j = 0; j < nv; j++)    {
		int vid = f.GetVertexInd(j);
		GeomVert jv = mMesh.GetGeomVertex(vid);
		if(v == jv)	
			j0 = j;
		x += jv;
	}
	x *= (1.0 / nv);
	GeomVert v1 = mMesh.GetGeomVertex(f.GetVertexInd((j0 + nv - 1) % nv));
	GeomVert v2 = mMesh.GetGeomVertex(f.GetVertexInd((j0 + nv + 1) % nv));

	GeomVert ev1 = (v + v1) * 0.5;
	GeomVert ev2 = (v + v2) * 0.5;


	return (x + ev1 + ev2 + v) * 0.25;

}

inline Mesh DooSabin(Mesh inMesh)
{
    Mesh outMesh;

    /////////////////for face face
    for(int i = 0; i < inMesh.GetNumberFacets(); i++) {
		vector<GeomVert> facet;
        TopoFacet f = inMesh.GetFacet(i);
		for(int j = 0; j < f.GetNumberVertices(); j++)  {
			GeomVert v = inMesh.GetGeomVertex( f.GetVertexInd(j));
            facet.push_back(GetFaceVertex(inMesh, v, f));
        }
        outMesh.AddFacet(facet);
    }

    for(int i = 0 ; i < inMesh.GetNumberEdges(); i++){
        ////////////////for edge face
		vector<GeomVert> facet;
        TopoEdge e = inMesh.GetEdge(i);
		TopoFacet f0 = inMesh.GetFacet(e.GetIncFacet(0));
        TopoFacet f1 = inMesh.GetFacet(e.GetIncFacet(1));
		GeomVert v0 = inMesh.GetGeomVertex(e.GetVertex(0));
        GeomVert v1 = inMesh.GetGeomVertex(e.GetVertex(1));
		facet.push_back(GetFaceVertex(inMesh, v0, f1));        
		facet.push_back(GetFaceVertex(inMesh, v1, f1));        
		facet.push_back(GetFaceVertex(inMesh, v1, f0));        
		facet.push_back(GetFaceVertex(inMesh, v0, f0));        
        outMesh.AddFacet(facet);
    }

	//vertex face
	for(int i = 0; i < inMesh.GetNumberVertices(); i++)	{
		vector<GeomVert> facet;
		TopoVert v = inMesh.GetVertex(i);
		GeomVert gv = inMesh.GetGeomVertex(i);
		for(int j = v.GetNumberIncFacets() - 1; j >= 0; j--)	{
			TopoFacet f = inMesh.GetFacet(v.GetIncFacet(j));
			facet.push_back(GetFaceVertex(inMesh, gv, f));
		}
        outMesh.AddFacet(facet);
	}

	return outMesh;
}

inline GeomVert FacePoint(TopoFacet facet, Mesh mMesh)
{
	GeomVert vf(0, 0, 0);
	int nv = facet.GetNumberVertices();
	for(int i = 0; i < nv; i++)	{
		vf += mMesh.GetGeomVertex(facet.GetVertexInd(i));
	}
	vf *= (1.0 / nv);
	return vf;
}

inline Mesh CatmullClark(Mesh inMesh)
{
    Mesh outMesh;

    /////////////////for face face
    for(int i = 0; i < inMesh.GetNumberFacets(); i++) {
		TopoFacet f = inMesh.GetFacet(i);
		GeomVert vf = FacePoint(f, inMesh);
		vector<GeomVert> ve;
		for(int j = 0; j < f.GetNumberEdges(); j++)	{
			TopoEdge e = inMesh.GetEdge(f.GetIncEdge(j));
			GeomVert v = inMesh.GetGeomVertex(e.GetVertex(0));
			GeomVert w = inMesh.GetGeomVertex(e.GetVertex(1));
			TopoFacet f1 = inMesh.GetFacet(e.GetIncFacet(0));
			TopoFacet f2 = inMesh.GetFacet(e.GetIncFacet(1));
			GeomVert vf1 = FacePoint(f1, inMesh);
			GeomVert vf2 = FacePoint(f2, inMesh);
			ve.push_back((v + w + vf1 + vf2) * 0.25);
		}
		vector<GeomVert> vv;
		for(int j = 0; j <f.GetNumberVertices(); j++)	{
			TopoVert v = inMesh.GetVertex(f.GetVertexInd(j));
			GeomVert vg = inMesh.GetGeomVertex(f.GetVertexInd(j));
			GeomVert Q(0, 0, 0);
			for(int k = 0; k < v.GetNumberIncFacets(); k++)	{
				TopoFacet fa = inMesh.GetFacet(v.GetIncFacet(k));
				GeomVert vfa = FacePoint(fa, inMesh);
				Q += vfa;
			}
			Q *= (1.0/v.GetNumberIncFacets());

			GeomVert R(0, 0, 0);
			int n = v.GetNumberIncEdges();
			for(int k = 0; k < v.GetNumberIncEdges(); k++)	{
				TopoEdge e = inMesh.GetEdge(v.GetIncEdge(k));
				GeomVert v1 = inMesh.GetGeomVertex(e.GetVertex(0));
				GeomVert v2 = inMesh.GetGeomVertex(e.GetVertex(1));
				R += ((v1 + v2) * 0.5);
			}
			R *= (1.0 / v.GetNumberIncEdges());

			vv.push_back((Q + R * 2 + vg * (n - 3) ) * (1.0 / n));
		}
		int ne = f.GetNumberEdges();
		for(int j = 0; j < f.GetNumberEdges(); j++)	{
			vector<GeomVert> facet;
			facet.push_back(ve[j]);
			facet.push_back(vv[j]);
			facet.push_back(ve[(j + 1) % ne]);
			facet.push_back(vf);
			outMesh.AddFacet(facet);
		}
    }

	return outMesh;
}



inline Mesh Loop(Mesh inMesh)
{
    Mesh outMesh;
	Mesh triMesh;
	triMesh = inMesh;
	//for(int i = 0; i < inMesh.GetNumberFacets(); i++)	{
	//	TopoFacet f = inMesh.GetFacet(i);
	//	GeomVert vf = FacePoint(f, inMesh);
	//	int nv = f.GetNumberVertices();
	//	for(int j = 0; j < nv; j++)	{
	//		vector<GeomVert> facet;
	//		facet.push_back(inMesh.GetGeomVertex(f.GetVertexInd(j)));
	//		facet.push_back(inMesh.GetGeomVertex(f.GetVertexInd((j + 1 + nv) % nv)));
	//		facet.push_back(vf);
	//		triMesh.AddFacet(facet);
	//	}
	//}

	for(int i = 0; i < triMesh.GetNumberFacets(); i++)	{
		TopoFacet f = inMesh.GetFacet(i);
		vector<GeomVert> ers;
		vector<GeomVert> vv;
		for(int j = 0; j < f.GetNumberEdges(); j++)	{
			TopoEdge e = inMesh.GetEdge(f.GetIncEdge(j));
			GeomVert r = inMesh.GetGeomVertex(e.GetVertex(0));
			GeomVert s = inMesh.GetGeomVertex(e.GetVertex(1));
			TopoFacet f1 = inMesh.GetFacet(e.GetIncFacet(0));
			TopoFacet f2 = inMesh.GetFacet(e.GetIncFacet(1));
			GeomVert p(0,0,0), q(0,0,0);
			for(int k = 0; k < f1.GetNumberVertices(); k++)	{
				GeomVert v0 = inMesh.GetGeomVertex(f1.GetVertexInd(k));
				if(!((v0 == r) || (v0 == s)))
					p = v0;
			}
			for(int k = 0; k < f2.GetNumberVertices(); k++)	{
				GeomVert v0 = inMesh.GetGeomVertex(f2.GetVertexInd(k));
				if(!((v0 == r) || (v0 == s)))
					q = v0;
			}
			ers.push_back( (p + r * 3 + s * 3 + q) * (1.0 / 8));
		}
		for(int j = 0; j < f.GetNumberVertices(); j++)	{
			TopoVert v = inMesh.GetVertex(f.GetVertexInd(j));
			GeomVert vg = inMesh.GetGeomVertex(f.GetVertexInd(j));
			GeomVert p(0, 0, 0);
			int n = v.GetNumberIncVertices();
			for(int k = 0; k < n; k++)	{
				p += inMesh.GetGeomVertex(v.GetIncVertex(k));
			}
			p *= (1.0 / n);
			float alpha = 0;
			if(n > 3)
				alpha = 3.0 / (8 * n);
			else
				alpha = 3.0 / 16;
			vv.push_back(p * alpha + vg * (1 - alpha));
		}
		vector<GeomVert> facet0;
		for(int j = 0; j < f.GetNumberEdges(); j++)	{
			vector<GeomVert> facet;
			facet.push_back(vv[j]);
			facet.push_back(ers[(j + 1) % 3]);
			facet.push_back(ers[j]);
			//facet.push_back(vv[(j + 1) % 3]);
			outMesh.AddFacet(facet);

			facet0.push_back(ers[j]);
		}
		outMesh.AddFacet(facet0);
	}

	return outMesh;
}


void Helper::GenerateDooSabin(int n)
{
	Mesh outMesh = _inputMesh;
	for(int i = 0; i < n; i++)	{
		outMesh = DooSabin(outMesh);
	}
	outMesh.WritePLY("data/DooSabin.ply");;
}

void Helper::GenerateCatmullClark(int n)
{
	Mesh outMesh = _inputMesh;
	for(int i = 0; i < n; i++)	{
		outMesh = CatmullClark(outMesh);
	}
	outMesh.WritePLY("data/CatmullClark.ply");;
}

void Helper::GenerateLoop(int n)
{
	Mesh outMesh = _inputMesh;
	for(int i = 0; i < n; i++)	{
		outMesh = Loop(outMesh);
	}
	outMesh.WritePLY("data/Loop.ply");;
}
	
void Helper::GenerateRevolution(int ns, int nu, int nv, int option)
{
    QVector<QPointF> curveF = GetCurve();
	Mesh mMesh;
    bool closed = (_splineType == 3);
    //mMesh = Revolution(curveF, ns, closed);

    mMesh = RevolutionBezier(curveF, ns, nu, nv, option);
    mMesh.WritePLY("data/revolution.ply");
}

void Helper::GenerateSweep()
{
	QVector<QPointF> trajectory = GetCurve();
	Mesh mMesh;
    mMesh = Sweep(_generator, trajectory);
    mMesh.WritePLY("data/sweep.ply");
}

void Helper::LoadMesh(const char* filename)
{
	Mesh newMesh;
	_inputMesh = newMesh;
    _inputMesh.ReadOFF(filename);//"data/octa.off");//"data/dodec.off");
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
        curveF = BezierF(ctrlPts, ctrlPts.size() * 2);
        break;
    case 1:
        curveF = CubicBSplineF(ctrlPts, 16);
        break;
    case 2:
        curveF = BezierSubdivideF(ctrlPtsF, 4, 0.5);
        break;
    case 3:
        curveF = ClosedBezierF(ctrlPts, ctrlPts.size() * 2);
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
        painter->drawPolyline(Bezier(ctrlPts, ctrlPts.size() * 2));
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
        painter->drawPolyline(ClosedBezier(ctrlPts, ctrlPts.size() * 2));
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

