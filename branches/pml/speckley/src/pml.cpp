
/*****************************************************************************
*
* Copyright (c) 2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development from 2016 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

/*
This file contains the additional information used to implement perfectly matched
layers in the code.
*/

#include <vector>
#include <string>
#include <iostream>
#include <speckley/pml.h>

int autoWidth(int n0, int d0);
int autoMax(int n0, int d0, int w);
float sigma(int n, int n0, int d0, int w, int maxSigma);

pml_class::pml_class(){
	North_PML = false;
	South_PML = false;
	East_PML = false;
	West_PML = false;
	width = 0;
}

pml_class::~pml_class(){

}

void pml_class::set_width(int w, int n0, int d0){
	if(w <= 1){
		w = autoWidth(n0, d0);
	} else {
		width = w;	
	}
	
}

void pml_class::set_on(std::vector<bool> pml){

	if(pml.size() != 4){
		std::cout << "Warning: pml vector is the wrong size" << std::endl;
	}

	North_PML = pml.at(0);
	South_PML = pml.at(1);
	East_PML = pml.at(3);
	West_PML = pml.at(2);

}

bool pml_class::check(std::string direction){

	if(direction == "North"){ return North_PML;}
	else if(direction == "South"){ return South_PML; } 
	else if(direction == "East"){ return East_PML; } 
	else if(direction == "West"){ return West_PML; }
	else {
		std::cout << "Warning: check the direction." << std::endl;
		return false;
	}
}

std::string pml_class::print_info(){
	std::string answer = "";
	answer += "Current Width = " + std::to_string(width) + ". ";

	answer += "North: " + std::to_string(North_PML) + ", ";
	answer += "South: " + std::to_string(South_PML) + ", ";
	answer += "East: " + std::to_string(East_PML) + ", ";
	answer += "West: " + std::to_string(West_PML) + ". ";

	return answer;
}

int autoWidth(int n0, int d0){
	// Sets the width of the pml
	// n0 is the number of elements in this dimension
	// d0 is the number of subdivisions in this dimension

	int w = (int) 0.10 * n0 * d0;

	return w; // AEAE Improve on this
}

int autoMax(int n0, int d0, int w){
	// Sets the max value of sigma on the boundary
	// n0 is the number of elements in this dimension
	// d0 is the number of subdivisions in this dimension
	// w is the width of the pml

	return 100; // AEAE Improve on this
}

// Sigma function
float sigma(int n, int n0, int d0, int w, int maxSigma){
	// n is the node we are examining
	// n0 is the number of elements in this dimension
	// d0 is the number of subdivisions in this dimension
	// w is the width of the pml
	// maxSigma is the maximum value of sigma on the boundary
	
	long totalNodes = n0 * d0;
	
	// If we're not in the pml then return zero
	if(n > w && n < (totalNodes - w))
		return 0;

	// If necessary, set w and d0
	if( w == 0 ){
		w = autoWidth(n0, d0);
	}
	if( maxSigma == 0 ){
		maxSigma = autoMax(n0, d0, w);
	}

	// Calculate how far n is from the boundary
	float R;
	if (n <= w) {
		R = (float) w - n;
	} else {
		R = (float) n - totalNodes + w + 1.0; // AETODO - Double check this
	}

	// Return the value
	return (float) maxSigma * (R / w) * (R / w);
}