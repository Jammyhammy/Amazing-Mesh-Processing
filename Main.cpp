#include <math.h>
#include <vector>
#include <iostream>
#include <list>
#include <string>

void subdivide_loop(std::string meshfile, std::string outmeshfile);
void fillHoles(std::string meshfile, std::string outmeshfile);
//std::string const meshfile = "camel.m";
//std::string const outmeshfile = "camel2.m";
//std::string const meshfile = "mannequin.m";
//std::string const outmeshfile = "mannequin2.m";

std::string const meshfile = "cube.m";
//std::string const outmeshfile = "tetra2.m";
std::string const outmeshfile = "cube2.m";


int main(int argc, char * argv[])
{
	if (std::string(argv[1]) == "-subdivision") {
		subdivide_loop(meshfile, outmeshfile);
	}
	else if (std::string(argv[1]) == "-fillhole")
	{
		fillHoles(meshfile, outmeshfile);
	}
	return 0;
}

