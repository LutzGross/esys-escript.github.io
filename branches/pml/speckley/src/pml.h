
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

// class pml_info;

// // public:
// // 	pml_info();
// // 	~pml_info();
// // };

//A class for keeping track of the pml information

#ifndef pml_class
// #define pml_class

class pml_class {
	private:
		bool North_PML, South_PML, East_PML, West_PML;
		int width;
		int dimensions;

	public:
		pml_class();
		~pml_class();
		void set_width(int w, int n0, int d0);
		void set_on(std::vector<bool> pml);
		bool check(std::string direction);
		std::string print_info();
};

#endif

// int autoWidth(int n0, int d0)

// int autoMax(int n0, int d0, int w)

// float sigma(int n, int n0, int d0, int w, int maxSigma)