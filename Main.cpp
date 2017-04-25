#include <math.h>
#include <vector>
#include <iostream>
#include <list>
#include <string>
#include "Parameterization.hpp"

void subdivide_loop(std::string meshfile, std::string outmeshfile);
void fillHoles(std::string meshfile, std::string outmeshfile);
void Parameterization(std::string meshfile, std::string outmeshfile);

std::string const meshfile = "smile.m";
std::string const outmeshfile = "smileuv.m";

int main(int argc, char * argv[])
{
	if (std::string(argv[1]) == "-subdivision") {
		subdivide_loop(argv[2], argv[3]);
		//subdivide_loop(meshfile, outmeshfile);
	}
	else if (std::string(argv[1]) == "-fillhole")
	{
		//fillHoles(meshfile, outmeshfile);
		fillHoles(argv[2], argv[3]);

	}
	else if (std::string(argv[1]) == "-parameterization") {
		Parameterization(argv[2], argv[3]);
	}
	return 0;
}
