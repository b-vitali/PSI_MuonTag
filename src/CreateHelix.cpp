#include <iostream>
#include <vector>
#include "TVector3.h"
#include "CreateHelix.hh"

////////////////////////////////////////////////////////////////////!
/*
    Main function to create the helix. 
    This calls all the previous and returns a G4TessellatedSolid
*/

// Create a triangulated helix using G4TessellatedSolid
G4TessellatedSolid* CreateHelix(G4String name, TVector3 center, double size, double runningangle, double length, int steps, double extrusion)
{
    // Create a new G4TessellatedSolid
    G4TessellatedSolid* helix = new G4TessellatedSolid(name);

    // Calculate number of turns
    double turns = AngleToTurns(runningangle, length, center.x());
	
    // Create base starting from 'center' and 'size'; triangulate it and add it to the helix
	std::vector<TVector3> base = CreateBase(center, size);
    // Tilt to metch the running angle
    base = Tilt(base, center, runningangle);

    // If no extrusion just cap it, otherwise use the AddExtrusion function
    if(extrusion == 0) {
        std::vector<std::tuple<TVector3, TVector3, TVector3>> triang_base = TriangulateBase(base);
	    MyADD(helix,triang_base);
    }
    else{
        std::vector<std::tuple<TVector3, TVector3, TVector3>> triang_base_extrusion = AddExtrusion(base, extrusion);
        MyADD(helix,triang_base_extrusion);
    }

	std::vector<TVector3> first=base;
	std::vector<TVector3> second;

	// Create the side: duplicate and move the base; triangulate the side; add to helix and repeat. 
	for(int i=0; i<steps; i++){
		second=Transform(first, turns, steps, center, length);
		std::vector<std::tuple<TVector3, TVector3, TVector3>> triang_side = TriangulateSide(first,second);
		MyADD(helix,triang_side);
		first = second; 
	}

	// The last vertex loop is use to cap the helix (NB it needs to be flipped)
    std::vector<TVector3> cap = first;
    cap = FlipLoop(cap);

    // If no extrusion just cap it, otherwise use the AddExtrusion function
    if(extrusion==0){
        std::vector<std::tuple<TVector3, TVector3, TVector3>> triang_cap = TriangulateBase(cap);
	    MyADD(helix,triang_cap);
	
    }
    else{
        std::vector<std::tuple<TVector3, TVector3, TVector3>> triang_cap_extrusion = AddExtrusion(cap, extrusion);
        MyADD(helix,triang_cap_extrusion);
    }
	
    // Close the solid to ensure correct surface normals
    helix->SetSolidClosed(true);

    return helix;
}

////////////////////////////////////////////////////////////////////!
/*
    Auxiliary functions to create the helix shape
*/

// From the angle to the number of turns
double AngleToTurns(double angle, double length, double R){
    double turns = tan(angle)*length/(2*M_PI*R);
    return turns;
}

// Given the center and the size of the fiber creates a squared base
// Change this to move from *squared* helix to other shape.
std::vector<TVector3>  CreateBase(TVector3 center, double size){
    std::vector<TVector3> base = {
    center + TVector3(-size*0.5, -size*0.5, 0),
    center + TVector3(size*0.5, -size*0.5, 0),
    center + TVector3(size*0.5, size*0.5, 0),
    center + TVector3(-size*0.5, size*0.5, 0),};
    return base;
}

// Evaluate the path using the length t. 
// Change this to move from *helix* to *other extrusion shape*
TVector3 Path(double t, double turns, TVector3 center, double length){
    double R = center.x();
    double x,y,z;
    x = R * cos(t/length* 2*M_PI *turns);
    y = R * sin(t/length* 2*M_PI *turns);
    z = t;
    TVector3 offset = TVector3(x,y,z);
    return offset;
}

// Function to flip the closing cap (to fix normals)
std::vector<TVector3> FlipLoop(std::vector<TVector3> cap){
    std::vector<TVector3> fliploop;
    fliploop.push_back(cap[1]);
    fliploop.push_back(cap[0]);
    fliploop.push_back(cap[3]);
    fliploop.push_back(cap[2]);
    return fliploop;
}

// Function to create the triangulation for the endcaps
std::vector<std::tuple<TVector3, TVector3, TVector3>> TriangulateBase(std::vector<TVector3>v){
    std::vector<std::tuple<TVector3, TVector3, TVector3>> faces;
    faces.push_back(std::make_tuple(v[3],v[1],v[0]));
    faces.push_back(std::make_tuple(v[3],v[2],v[1]));
    return faces;
}

// Function to create the triangulation given two verteces loops to be bridged
std::vector<std::tuple<TVector3, TVector3, TVector3>> TriangulateSide(std::vector<TVector3>v, std::vector<TVector3>u){
    std::vector<std::tuple<TVector3, TVector3, TVector3>> faces;
    for(int i = 0; i<v.size()-1;i++){
        faces.push_back(std::make_tuple(v[i],v[i+1],u[i]));
        faces.push_back(std::make_tuple(v[i+1],u[i+1],u[i]));
    }
    // To close the loop
    faces.push_back(std::make_tuple(v[v.size()-1],v[0],u[v.size()-1]));
    faces.push_back(std::make_tuple(v[0],u[0],u[v.size()-1]));
    return faces;
}

// The square needs to be alligned to the angle of the fiber
std::vector<TVector3> Tilt(std::vector<TVector3>v,  TVector3 center, double angle){
    double x,y,z;
    std::vector<TVector3> u;

    for(int i = 0; i<v.size();i++){
        x = v[i].x();
        y = v[i].y()*cos(angle)+v[i].z()*sin(angle);
        z = -v[i].y()*sin(angle)+v[i].z()*cos(angle);
             
        u.push_back(TVector3(x,y,z));
    }
    return u;
}

// Function to add a triangulation to a G4TessellatedSolid.
void MyADD(G4TessellatedSolid* helix, std::vector<std::tuple<TVector3, TVector3, TVector3>> triang){
	TVector3 a,b,c;
	for(auto i : triang){
		a = std::get<0>(i);
		b = std::get<1>(i);
		c = std::get<2>(i);
		G4ThreeVector va(a.x(), a.y(), a.z());
        G4ThreeVector vb(b.x(), b.y(), b.z());
        G4ThreeVector vc(c.x(), c.y(), c.z());
        helix->AddFacet(new G4TriangularFacet(va, vb, vc, ABSOLUTE));
	}
}

// Given a lopp create a second loop along the path. These two are going to be bridged with triangulation
std::vector<TVector3> Transform(std::vector<TVector3>v, double turns, int steps, TVector3 center, double length){
    std::vector<TVector3> u;
    double angle_step = turns/steps*2*M_PI;
    double length_step = length/steps;
    for(int i = 0; i<v.size();i++){
        v[i].RotateZ(angle_step);
        u.push_back(v[i]+length_step*TVector3(0,0,1));
    }
    return u;
}

// Function to flip the closing cap (to fix normals)
std::vector<std::tuple<TVector3, TVector3, TVector3>> AddExtrusion(std::vector<TVector3> base, double extrusion){
    std::vector<std::tuple<TVector3, TVector3, TVector3>> faces;
    std::vector<TVector3> extrusion_loop;
    TVector3 direction =((base[2]-base[0])).Cross(base[1]-base[0]);
    for(int i = 0; i<base.size();i++){
        extrusion_loop.push_back(base[i] + extrusion * direction.Unit());
    }

    std::vector<std::tuple<TVector3, TVector3, TVector3>>  cap = TriangulateBase(extrusion_loop);

    base = FlipLoop(base);
    extrusion_loop = FlipLoop(extrusion_loop);
    faces = TriangulateSide(base, extrusion_loop);

    faces.insert(faces.end(), cap.begin(), cap.end());

    return faces;
}