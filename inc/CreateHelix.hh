
#ifndef CreateHelix_h
#define CreateHelix_h 1

#include <iostream>
#include <vector>
#include "TVector3.h"

#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
// My additional file to create the G4TessellatedSolid

//!
#include "G4VUserDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4TwistedBox.hh"
#include "G4IntersectionSolid.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "VirtualDetectorSD.hh"
#include "ScintSD.hh"
#include "SiPMSD.hh"

#include "DetectorMessenger.hh"

#include "G4GlobalMagFieldMessenger.hh"

#include "G4Cache.hh"
#include "G4SDManager.hh"

#include "G4SubtractionSolid.hh"

namespace HelixMaker{
    ////////////////////////////////////////////////////////////////////!
    /*
        Two functions for debug: printvertex and printvector(of verteces)
    */
    void PrintVertex(TVector3 v);

    void PrintVector(std::vector<TVector3> v);

    void print_progress_bar(double percentage, std::string name);

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

    ////////////////////////////////////////////////////////////////////!
    /*
        Main function to create the helix. 
        This calls all the previous and returns a G4TessellatedSolid
    */

    // Create a triangulated helix using G4TessellatedSolid
    G4TessellatedSolid* CreateHelix(G4String name, TVector3 center, double size, double runningangle, double length, int steps, double extrusion); 
}

bool print = false
void start_print(G4String s);
void running_print(G4String s);
void finish_print(G4String s);

#endif