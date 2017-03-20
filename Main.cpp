#include <math.h>
#include <vector>
#include <iostream>
#include <list>
#include <string>

void subdivide_loop(std::string meshfile, std::string outmeshfile);
void fillHoles(std::string meshfile, std::string outmeshfile);

std::string const meshfile = "camel.m";
std::string const outmeshfile = "camel2.m";

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
	return 0;
}
