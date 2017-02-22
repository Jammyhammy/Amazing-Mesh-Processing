#include <math.h>
#include <iostream>

#include "./lib/Mesh_Library/Mesh/mesh.h"
#include "./lib/Mesh_Library/Mesh/iterators.h"

using namespace MeshLib;

int main(int argc, char * argv[])
{

	// check http://www.ltcconline.net/greenl/courses/107/vectors/dotcros.htm for reference of vector operations
	if (strcmp(argv[1], "-example") == 0)
	{
		CMesh mesh;
		mesh.read_m(argv[2]);		//read into mesh (triangle mesh with format .m)

		double area = 0;
		for (MeshFaceIterator fiter(&mesh); !fiter.end(); fiter++) //face iterator; go through each face of a mesh
		{
			CFace * pF = *fiter;
			CPoint p0 = pF->halfedge()->source()->point();
			CPoint p1 = pF->halfedge()->target()->point();
			CPoint p2 = pF->halfedge()->he_prev()->source()->point();
			area += ((p1 - p0) ^ (p2 - p1)).norm() * 0.5;  //area of a triangle = half of the length of the cross product of the two vectors = half of the area of the parallelogram determined by the two vectors 
		}
		area = area / mesh.numFaces(); //avearge face area

		double dihedral_angle = 0;
		for (MeshEdgeIterator eiter(&mesh); !eiter.end(); eiter++) //edge iterator; go through each edge of a mesh
		{
			CEdge * pE = *eiter;
			CFace * pF0 = pE->halfedge(0)->face(); // a triangle face neighboring edge pE
			CPoint p0 = pF0->halfedge()->source()->point();
			CPoint p1 = pF0->halfedge()->target()->point();
			CPoint p2 = pF0->halfedge()->he_prev()->source()->point();
			CPoint normal0 = (p1 - p0) ^ (p2 - p1); // compute normal of the triangle face

			if (pE->boundary())
				continue;
			CFace * pF1 = pE->halfedge(1)->face(); // the other triangle face neighboring edge pE if pE is not boundary edge
			p0 = pF1->halfedge()->source()->point();
			p1 = pF1->halfedge()->target()->point();
			p2 = pF1->halfedge()->he_prev()->source()->point();
			CPoint normal1 = (p1 - p0) ^ (p2 - p1); // compute normal of the trianlge face

													//the angle between two face normals = the angle of their dot product
			dihedral_angle += std::acos((normal0*normal1) / (normal0.norm() * normal1.norm()));
		}
		dihedral_angle = dihedral_angle / mesh.numEdges(); //avearge dihedral_angle

		printf("the avearge face area %f the average dihedral_angle %f ", area, dihedral_angle);

		mesh.write_m(argv[3]); //write the mesh
	}

	return 0;
}

