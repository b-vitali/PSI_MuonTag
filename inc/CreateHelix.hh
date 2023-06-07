#include <iostream>
#include <vector>
#include "TVector3.h"

////////////////////////////////////////////////////////////////////!
/*
    Two functions for debug: printvertex and printvector(of verteces)
*/
void PrintVertex(TVector3 v){
    G4cout<<"vertex : "<<v[0]<<" "<<v[1]<<" "<<v[2]<<G4endl;
}

void PrintVector(std::vector<TVector3> v){
    G4cout<<G4endl<<"Vector length = "<<v.size()<<G4endl;
    for(int i = 0; i<v.size();i++){
        PrintVertex(v[i]);
    }
    G4cout<<G4endl;
}

////////////////////////////////////////////////////////////////////!
/*
    Auxiliary functions to create the helix shape
*/

// From the angle to the number of turns
double AngleToTurns(double angle, double length, double R);

// Given the center and the size of the fiber creates a squared base
// Change this to move from *squared* helix to other shape.
std::vector<TVector3>  CreateBase(TVector3 center, double size);

// Evaluate the path using the length t. 
// Change this to move from *helix* to other extrusion shape
TVector3 Path(double t, double turns, TVector3 center, double length);

// Function to flip the closing cap (to fix normals)
std::vector<TVector3> FlipLoop(std::vector<TVector3> cap);

// Function to create the triangulation for the endcaps
std::vector<std::tuple<TVector3, TVector3, TVector3>> TriangulateBase(std::vector<TVector3>v);

// Function to create the triangulation given two verteces loops to be bridged
std::vector<std::tuple<TVector3, TVector3, TVector3>> TriangulateSide(std::vector<TVector3>v, std::vector<TVector3>u);

// The square needs to be alligned to the angle of the fiber
std::vector<TVector3> Tilt(std::vector<TVector3>v,  TVector3 center, double angle);

// Function to add a triangulation to a G4TessellatedSolid.
void MyADD(G4TessellatedSolid* helix, std::vector<std::tuple<TVector3, TVector3, TVector3>> triang);

// Given a lopp create a second loop along the path. These two are going to be bridged with triangulation
std::vector<TVector3> Transform(std::vector<TVector3>v, double turns, int steps, TVector3 center, double length);

// Function to flip the closing cap (to fix normals)
std::vector<std::tuple<TVector3, TVector3, TVector3>> AddExtrusion(std::vector<TVector3> base, double extrusion);
