#ifndef __FORTRANIO_H__
#define __FORTRANIO_H__

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <complex>
#include <limits>

using namespace std;


class FortranIO {

public:

	static string& read_string(ifstream& infile);
	static int read_int(ifstream& infile);
	static double read_double(ifstream& infile);
	static void skip_line(ifstream& file);

};

#endif
