// ------------------------------------------------------------
// Mesh.cpp:  Implementation of mesh class
// ------------------------------------------------------------

#include "mesh.h"
#include <fstream>
#include <string>
#include <stdlib.h>
#include <iostream>

using namespace std;

// ------------------------------------------------------------
// AddFacet:  Adds a triangle to the mesh.
//            This is one of 2 functions that can be used to build a mesh
// ------------------------------------------------------------
void Mesh::AddFacet(datatype x1, datatype y1, datatype z1, datatype x2, datatype y2, datatype z2, datatype x3, datatype y3, datatype z3) {
	
	vector<GeomVert> geomfacet;
	geomfacet.push_back( GeomVert(x1, y1, z1) );
	geomfacet.push_back( GeomVert(x2, y2, z2) );
	geomfacet.push_back( GeomVert(x3, y3, z3) );

	AddFacet( geomfacet );
}
// ------------------------------------------------------------





// ------------------------------------------------------------
// AddFacet:  Adds a facet with arbitrary number of vertices to mesh
//            This is one of 2 functions that can be used to build a mesh
// ------------------------------------------------------------
void Mesh::AddFacet(vector<GeomVert> geomfacet) {
	int i;

	// --------------
	// Create topo facet (list of geom vertex indices)
	TopoFacet topofacet;
	// Look for facet vertices in mesh - if they don't already exist in mesh then add them
	for (i = 0; i < geomfacet.size(); i++) {
		int v_ind = FindGeomVertex( geomfacet[i] );
		if (v_ind == -1) {
			// New vertex:  add geomtric vertex
			v_ind = mGeomVerts.size();
			mGeomVerts.push_back( geomfacet[i] );

			// Add topo vertex
			TopoVert topovert;
			mTopoVerts.push_back( topovert );
		}

		// Add vertex indice to topo facet
		topofacet.AddIncVertex( v_ind );
	}

	// Add this new topo facet to mesh	
	int facet_ind = mTopoFacets.size();
	mTopoFacets.push_back( topofacet );


	// Add edges of facet to mesh, again checking if they already exist
	for (i = 0; i < topofacet.GetNumberVertices(); i++) {
		int prev = (i == 0) ? topofacet.GetNumberVertices() - 1  :  i - 1;
		
		// Create edge
		TopoEdge e;
		e.SetVertex(0, topofacet.GetVertexInd(prev) );
		e.SetVertex(1, topofacet.GetVertexInd(i) );

		// Check if exists
		int e_ind = FindTopoEdge( e );
		
		if (e_ind == -1) {
			// Didn't exist, add to mesh
			e_ind = mTopoEdges.size();
			mTopoVerts[ e.GetVertex(0) ].AddIncEdge( e_ind );
			mTopoVerts[ e.GetVertex(1) ].AddIncEdge( e_ind );
			mTopoEdges.push_back( e );			
		}

		// Point edge to this facet
		mTopoEdges[e_ind].AddIncFacet( facet_ind );

		// Point facet to this edge
		mTopoFacets[ facet_ind ].AddIncEdge( e_ind );
	}
	// --------------
		

	
	// Compute other connectivity
	for (i = 0; i < topofacet.GetNumberVertices(); i++) {
		// Add vertex-facet topology
		mTopoVerts[  topofacet.GetVertexInd(i) ].AddIncFacet( facet_ind );

		// Add vertex-vertex (edge) topology
		int prev = (i == 0) ? topofacet.GetNumberVertices() - 1  :  i - 1;
		int next = (i == topofacet.GetNumberVertices() - 1) ? 0 : i + 1;

		mTopoVerts[ topofacet.GetVertexInd(i) ].AddIncVert( topofacet.GetVertexInd( prev ) );
		mTopoVerts[ topofacet.GetVertexInd(i) ].AddIncVert( topofacet.GetVertexInd( next ) );
	}
	
	// Facet-facet adjacency...
	for (i = 0; i < mTopoFacets[ facet_ind ].GetNumberEdges(); i++) {		
		TopoEdge edge = mTopoEdges[ mTopoFacets[ facet_ind ].GetIncEdge(i) ];
		for (int j = 0; j < edge.GetNumberIncFacets(); j++) {
			if (edge.GetIncFacet(j) != facet_ind) {
				mTopoFacets[ facet_ind ].AddIncFacet( edge.GetIncFacet(j) );
				mTopoFacets[ edge.GetIncFacet(j) ].AddIncFacet( facet_ind );
			}
		}
	}
}
// ------------------------------------------------------------





// ------------------------------------------------------------
// Erase:  Releases all memory used by object
// ------------------------------------------------------------
void Mesh::Erase() {
	mGeomVerts.clear();
	mTopoVerts.clear();
	mTopoEdges.clear();
	mTopoFacets.clear();
}
// ------------------------------------------------------------





// ------------------------------------------------------------
// FindGeomVertex:  Searches for a geometric vertex in the mesh,
//                  returning its indice if found, -1 otherwise
// ------------------------------------------------------------
int Mesh::FindGeomVertex(GeomVert v) {
	for (int i = 0; i < mGeomVerts.size(); i++) {
		if (mGeomVerts[i] == v) return i;
	}
	return -1;
}
// ------------------------------------------------------------





// ------------------------------------------------------------
// FindTopoEdge:  Searches for an edge in the mesh, returing 
//                its indice if found, -1 otherwise
// ------------------------------------------------------------
int	Mesh::FindTopoEdge(TopoEdge e) {
	for (int i = 0; i < mTopoEdges.size(); i++) {
		if (mTopoEdges[i] == e) return i;
	}
	return -1;
}
// ------------------------------------------------------------


void Mesh::WritePLY(char* filename)
{
    ofstream ofs;
    ofs.open(filename, ofstream::out);
    ofs<<"ply\nformat ascii 1.0\nelement vertex ";
    ofs<<GetNumberVertices();
    ofs<<"\nproperty float32 x\nproperty float32 y\nproperty float32 z\nelement face ";
    ofs<<GetNumberFacets();
    ofs<<"\nproperty list uint8 int32 vertex_indices\nend_header\n";
    for(int i = 0; i < GetNumberVertices(); i++)    {
        ofs<<GetGeomVertex(i).GetCo(0)<<" "<<
             GetGeomVertex(i).GetCo(1)<<" "<<
             GetGeomVertex(i).GetCo(2)<<endl;
    }
    for(int i = 0; i < GetNumberFacets(); i++)  {
		int nv = GetFacet(i).GetNumberVertices();
		ofs<<nv << " ";
		for(int j = 0; j < nv; j++)	{
             ofs << GetFacet(i).GetVertexInd(j)<<" ";
		}
		ofs<<endl;
    }
}

//void Mesh::ReadOFF(char* filename)
//{
//	ifstream in;
//	in.open(filename);
//	string readLine;
//	// Check if file is in OFF format
//	getline(in,readLine);
//	if (readLine != "OFF")
//	{
//	cout << "The file to read is not in OFF format." << endl;
//	return;
//	}
//	int nv, nf;
//	int delimiterPos_1, delimiterPos_2, delimiterPos_3, delimiterPos_4;
//	// Read values for Nv and Nf
//	getline(in,readLine);
//	delimiterPos_1 = readLine.find(" ", 0);
//	nv = atoi(readLine.substr(0,delimiterPos_1+1).c_str());
//	delimiterPos_2 = readLine.find(" ", delimiterPos_1);
//	nf = atoi(readLine.substr(delimiterPos_1,delimiterPos_2 +1).c_str());
//
//	for (int n=0; n<nv; n++)
//	{
//		getline(in,readLine);
//		delimiterPos_1 = readLine.find(" ", 0);
//		vertices[n].x = atof(readLine.substr(0,delimiterPos_1).c_str());
//		delimiterPos_2 = readLine.find(" ", delimiterPos_1+1);
//		vertices[n].y =
//		atof(readLine.substr(delimiterPos_1,delimiterPos_2 ).c_str());
//		delimiterPos_3 = readLine.find(" ", delimiterPos_2+1);
//		vertices[n].z =
//		atof(readLine.substr(delimiterPos_2,delimiterPos_3 ).c_str());
//
//		cout << vertices[n].x << "\t" << vertices[n].y << "\t" <<
//		vertices[n].z << "\t" << endl;
//	}
//
//	// Read the facades
////	facades = new facade[nf];
//
//	for (int n=0; n<nf; n++)
//	{
//		getline(in,readLine);
//		delimiterPos_1 = readLine.find(" ", 0);
//		delimiterPos_2 = readLine.find(" ", delimiterPos_1+1);
//		facades[n].v1 =
//		atoi(readLine.substr(delimiterPos_1,delimiterPos_2 ).c_str());
//		delimiterPos_3 = readLine.find(" ", delimiterPos_2+1);
//		facades[n].v2 =
//		atoi(readLine.substr(delimiterPos_2,delimiterPos_3 ).c_str());
//		delimiterPos_4 = readLine.find(" ", delimiterPos_3+1);
//		facades[n].v3 =
//		atoi(readLine.substr(delimiterPos_3,delimiterPos_4 ).c_str());
//
//		cout << facades[n].v1 << "\t" << facades[n].v2 << "\t" <<
//		facades[n].v3 << "\t" << endl;
//	}
//}


void Mesh::WriteASCII(char* filename)
{
    ofstream ofs;
    ofs.open(filename, ofstream::out);
    ofs<<GetNumberVertices()<<" "<<GetNumberFacets()<<endl;
    for(int i = 0; i < GetNumberVertices(); i++)    {
        ofs<<GetGeomVertex(i).GetCo(0)<<" "<<
             GetGeomVertex(i).GetCo(1)<<" "<<
             GetGeomVertex(i).GetCo(2)<<endl;
    }
    for(int i = 0; i < GetNumberFacets(); i++)  {
        ofs<<"3 "<<
             GetFacet(i).GetVertexInd(0)<<" "<<
             GetFacet(i).GetVertexInd(1)<<" "<<
             GetFacet(i).GetVertexInd(2)<<endl;
    }
}

