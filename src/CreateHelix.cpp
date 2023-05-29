#include <iostream>
#include <vector>
#include "TVector3.h"

////////////////////////////////////////////////////////////////////!
/*
    Two functions for debug: printvertex and printvector(of verteces)
*/
void printvertex(TVector3 v){
    G4cout<<"vertex : "<<v[0]<<" "<<v[1]<<" "<<v[2]<<G4endl;
}

void printvector(std::vector<TVector3> v){
    G4cout<<G4endl<<"Vector length = "<<v.size()<<G4endl;
    for(int i = 0; i<v.size();i++){
        printvertex(v[i]);
    }
    G4cout<<G4endl;
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
std::vector<TVector3>  createbase(TVector3 center, double size){
    std::vector<TVector3> base = {
    center + TVector3(-size*0.5, -size*0.5, 0),
    center + TVector3(size*0.5, -size*0.5, 0),
    center + TVector3(size*0.5, size*0.5, 0),
    center + TVector3(-size*0.5, size*0.5, 0),};
    return base;
}

// Evaluate the path using the length t. 
// Change this to move from *helix* to *other extrusion shape*
TVector3 path(double t, double turns, TVector3 center, double length){
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
std::vector<std::tuple<TVector3, TVector3, TVector3>> triangulatebaseG4(std::vector<TVector3>v){
    std::vector<std::tuple<TVector3, TVector3, TVector3>> faces;
    faces.push_back(std::make_tuple(v[3],v[1],v[0]));
    faces.push_back(std::make_tuple(v[3],v[2],v[1]));
    return faces;
}

// Function to create the triangulation given two verteces loops to be bridged
std::vector<std::tuple<TVector3, TVector3, TVector3>> triangulatesideG4(std::vector<TVector3>v, std::vector<TVector3>u){
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


// Function to flip the closing cap (to fix normals)
std::vector<std::tuple<TVector3, TVector3, TVector3>> AddExtrusion(std::vector<TVector3> base, double extrusion){
    std::vector<std::tuple<TVector3, TVector3, TVector3>> faces;
    std::vector<TVector3> extrusion_loop;
    for(int i = 0; i<base.size();i++){
        extrusion_loop.push_back(base[i] + TVector3(0,0,extrusion));
    }
    printvector(base);
    printvector(extrusion_loop);
   
    std::vector<std::tuple<TVector3, TVector3, TVector3>>  cap = triangulatebaseG4(extrusion_loop);

    base = FlipLoop(base);
    extrusion_loop = FlipLoop(extrusion_loop);
    faces = triangulatesideG4(base, extrusion_loop);

    faces.insert(faces.end(), cap.begin(), cap.end());

    return faces;
}

// Extrude a loop of verteces around the path
std::vector<TVector3> translate(std::vector<TVector3>v, double turns, int steps, TVector3 center, double length){
    std::vector<TVector3> u;
    for(int i = 0; i<v.size();i++){
        TVector3 offset = path(v[i].z()+length/steps, turns, center, length) - path(v[i].z(), turns, center, length);
        u.push_back(v[i] + offset);
    }
    return u;
}

// The square needs to be translated and rotated around the path
std::vector<TVector3> rotate(std::vector<TVector3>v, double turns, int steps, TVector3 center, double length){
    double x,y,z;
    std::vector<TVector3> u;
    TVector3 rotationcenter = path(v[0].z(), turns, center, length);
    
    double angle = 2*M_PI / steps * turns;

    for(int i = 0; i<v.size();i++){
        x = v[i].x();
        y = v[i].y();
        z = v[i].z();

        // Translate the point to the origin
        double translatedX = x - rotationcenter.x();
        double translatedY = y - rotationcenter.y();

        // Apply the rotation around the z-axis
        double rotatedX = translatedX * cos(angle) - translatedY * sin(angle);
        double rotatedY = translatedX * sin(angle) + translatedY * cos(angle);

        // Translate the point back to the original position
        x = rotatedX + rotationcenter.x();
        y = rotatedY + rotationcenter.y();
        
        u.push_back(TVector3(x,y,z));
    }
    return u;
}

// Given a lopp create a second loop along the path. These two are going to be bridged with triangulation
std::vector<TVector3> transform(std::vector<TVector3>v, double turns, int steps, TVector3 center, double length){
    std::vector<TVector3> u = translate(v, turns, steps, center, length);
    u = rotate(u, turns, steps, center, length);
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

    double turns = AngleToTurns(runningangle, length, center.x());
    G4cout<<"number of turns for "<<name<<" : "<<turns<<G4endl;
	
    // Create base starting from 'center' and 'size'; triangulate it and add it to the helix
	std::vector<TVector3> base = createbase(center, size);
	
    // If no extrusion just cap it, otherwise use the AddExtrusion function
    if(extrusion == 0) {
        std::vector<std::tuple<TVector3, TVector3, TVector3>> triang_base = triangulatebaseG4(base);
	    MyADD(helix,triang_base);
    }
    else{
        std::vector<std::tuple<TVector3, TVector3, TVector3>> triang_base_extrusion = AddExtrusion(base, -extrusion);
        MyADD(helix,triang_base_extrusion);
    }

	std::vector<TVector3> first=base;
	std::vector<TVector3> second;

	// Create the side: duplicate and move the base; triangulate the side; add to helix and repeat. 
	for(int i=0; i<steps; i++){
		second=transform(first, turns, steps, center, length);
		//printvector(first);
		//printvector(second);
		std::vector<std::tuple<TVector3, TVector3, TVector3>> triang_side = triangulatesideG4(first,second);
		MyADD(helix,triang_side);
		//G4cout<<"done"<<G4endl;
		first = second; 
	}

	// The last vertex loop is use to cap the helix (NB it needs to be flipped)
    std::vector<TVector3> cap = first;
    cap = FlipLoop(cap);

    // If no extrusion just cap it, otherwise use the AddExtrusion function
    if(extrusion==0){
        std::vector<std::tuple<TVector3, TVector3, TVector3>> triang_cap = triangulatebaseG4(cap);
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