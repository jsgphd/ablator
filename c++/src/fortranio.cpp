#include "fortranio.h"

string& FortranIO::read_string(ifstream& infile) {
	string* line = new string();
	//infile.getline(line->c_str());
	getline(infile, *line);
	return *line;
}

int FortranIO::read_int(ifstream& infile) {
	int num;
	infile >> num;
	infile.ignore(numeric_limits<streamsize>::max(),'\n');
	return num;
}

double FortranIO::read_double(ifstream& infile) {
  double num;  
  infile >> num;  
  infile.ignore(numeric_limits<streamsize>::max(),'\n');
  return num;  
}

void FortranIO::skip_line(ifstream& infile) {
  infile.ignore(numeric_limits<streamsize>::max(),'\n');
}
