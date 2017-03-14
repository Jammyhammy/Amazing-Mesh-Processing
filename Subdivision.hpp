#pragma once
#ifndef SUBDIVISION_HPP
#define SUBDIVISION_HPP

#include <math.h>
#include <vector>
#include <iostream>
#include <list>
#include <string>


#include "./lib/Mesh_Library/Mesh/mesh.h"
#include "./lib/Mesh_Library/Mesh/iterators.h"
#include "./lib/mesh_library/Geometry/Point.h"

using namespace MeshLib;

// check http://www.ltcconline.net/greenl/courses/107/vectors/dotcros.htm for reference of vector operations

double inline roundDouble(double v) {
	v *= 10000.0;
	if (v < 0) return ceil(v - 0.5) / 10000.0;
	return floor(v + 0.5) / 10000.0;
}

CPoint* roundPoint(CPoint* cp) {

	auto npx = roundDouble(cp->get_x());
	auto npy = roundDouble(cp->get_y());
	auto npz = roundDouble(cp->get_z());
}

void inline getAdjacentVertices(CVertex *v, std::vector<CVertex*> *adj_v)
{
	if (!v->boundary())
	{
		auto h1 = v->halfedge()->he_next();
		auto cv = h1->vertex();
		adj_v->push_back(cv);

		auto vid = h1->vertex()->id();

		do {
			h1 = h1->he_sym()->he_next();
			cv = h1->vertex();
			adj_v->push_back(cv);
		} while (h1->he_sym()->he_next()->vertex()->id() != vid);
		auto t = "test";
	}
	else {

		auto h1 = v->halfedge()->he_next();
		auto cv = h1->vertex();
		adj_v->push_back(cv);

		auto vid = h1->vertex()->id();
		//TODO: Find some other way other than null checking in the while loop because this won't work for some.
		do {
			h1 = h1->he_sym()->he_next();
			cv = h1->vertex();
			adj_v->push_back(cv);
		} while (h1->he_sym() && !h1->he_sym()->he_next()->vertex()->id() != vid);
		auto t = "test";
	}
}

int inline getNewVertexIndexByOldVertexes(std::vector<std::pair<CVertex*, CVertex*>>& vertexes_pair, CVertex* v0, CVertex* v1, int old_vtx_num) {

	for (int i = 0; i <= vertexes_pair.size(); i++) {

		std::pair<CVertex*, CVertex*> pr = vertexes_pair.at(i);
		if ((pr.first->id() == v0->id() && pr.second->id() == v1->id())
			|| (pr.first->id() == v1->id() && pr.second->id() == v0->id())) {
			return old_vtx_num + i + 1;
		}
	}

	return -1;
}

void subdivide_loop(std::string meshfile, std::string outmeshfile)
{

	CMesh mesh;
	CMesh mesh2;
	mesh2.write_m("testempty.m");

	
	mesh.read_m(meshfile.c_str());		//read into mesh (triangle mesh with format .m)
	std::list<CVertex*> mvs = mesh.vertices();
	std::list<CFace*> mfs = mesh.faces();
	std::list<CEdge*> mes = mesh.edges();


	if (mvs.size() == 0) {
		throw std::invalid_argument("Mesh doesn't exist.");
	}
	std::vector<CVertex*> verts{ std::begin(mvs), std::end(mvs) };
	std::vector<CFace*> faces{ std::begin(mfs), std::end(mfs) };
	std::vector<CEdge*> edges{ std::begin(mes), std::end(mes) };

	std::vector<CVertex*> newvtx;
	std::vector<CVertex*> newnewvtx;

	std::vector<std::pair<CVertex*, CVertex*>> vtxpair;

	const int vertexNum = mesh.max_vertex_id();

	for (std::vector<CVertex*>::iterator it = verts.begin(); it != verts.end(); ++it) {
		CVertex* item = *it;
		newvtx.push_back(item);
		//CVertex* nv = mesh2.createVertex(item->id());
		//CPoint* np = &(item->point());
		//nv->set_point(*np);
	}

	//TODO: Maybe need to set updated on edge rather than half-edge.


	for (std::vector<CFace*>::iterator it = faces.begin(); it != faces.end(); ++it) {
		CFace* item = *it;
		auto he0 = item->halfedge();
		auto v0 = he0->source();

		auto he1 = he0->he_next();
		auto v1 = he1->source();

		auto he2 = he1->he_next();
		auto v2 = he2->source();

		if (mesh.isBoundary(he0) && !he0->updated()) {
			CVertex* nv = mesh.createVertex(mesh.max_vertex_id() + 1);
			//CVertex* nv2 = mesh2.createVertex(mesh2.max_vertex_id() + 1);

			CPoint* np = new CPoint(
				(v0->point()[0] + v1->point()[0]) / 2.0,
				(v0->point()[1] + v1->point()[1]) / 2.0,
				(v0->point()[2] + v1->point()[2]) / 2.0
			);
			nv->set_point(*np);
			//nv2->set_point(*np);
			newvtx.push_back(nv);
			//newnewvtx.push_back(nv);


			std::pair<CVertex*, CVertex*> pair(v0, v1);
			vtxpair.push_back(pair);

			he0->set_updated(true);
		}

		else if (!he0->updated()) {
			CHalfEdge* he0_prev = he0->he_prev();
			CHalfEdge* he0_prev_opp = he0->he_sym()->he_prev();

			CVertex* vtxtem1 = he0_prev->source();
			CVertex* vtxtem2 = he0_prev_opp->source();

			CVertex* nv = mesh.createVertex(mesh.max_vertex_id() + 1);
			//CVertex* nv2 = mesh2.createVertex(mesh2.max_vertex_id() + 1);

			CPoint* np = new CPoint(
				3.0 * (v0->point()[0] + v1->point()[0]) / 8.0 + (vtxtem1->point()[0] + vtxtem2->point()[0]) / 8.0,
				3.0 * (v0->point()[1] + v1->point()[1]) / 8.0 + (vtxtem1->point()[1] + vtxtem2->point()[1]) / 8.0,
				3.0 * (v0->point()[2] + v1->point()[2]) / 8.0 + (vtxtem1->point()[2] + vtxtem2->point()[2]) / 8.0
			);
			nv->set_point(*np);
			//nv2->set_point(*np);
			newvtx.push_back(nv);
			//newnewvtx.push_back(nv);


			std::pair<CVertex*, CVertex*> pair(v0, v1);
			vtxpair.push_back(pair);
			he0->set_updated(true);
			he0->he_sym()->set_updated(true);

		}

		if (mesh.isBoundary(he1) && !he1->updated()) {
			CVertex* nv = mesh.createVertex(mesh.max_vertex_id() + 1);
			//CVertex* nv2 = mesh2.createVertex(mesh2.max_vertex_id() + 1);

			CPoint* np = new CPoint(
				(v1->point()[0] + v2->point()[0]) / 2.0,
				(v1->point()[1] + v2->point()[1]) / 2.0,
				(v1->point()[2] + v2->point()[2]) / 2.0
			);
			nv->set_point(*np);
			//nv2->set_point(*np);
			newvtx.push_back(nv);

			std::pair<CVertex*, CVertex*> pair(v1, v2);
			vtxpair.push_back(pair);

			he0->set_updated(true);
		}
		else if (!he1->updated()) {
			CHalfEdge* he1_prev = he1->he_prev();
			CHalfEdge* he1_prev_opp = he1->he_sym()->he_prev();

			CVertex* vtxtem1 = he1_prev->source();
			CVertex* vtxtem2 = he1_prev_opp->source();

			CVertex* nv = mesh.createVertex(mesh.max_vertex_id() + 1);
			//CVertex* nv2 = mesh2.createVertex(mesh2.max_vertex_id() + 1);

			CPoint* np = new CPoint(
				3.0 * (v1->point()[0] + v2->point()[0]) / 8.0 + (vtxtem1->point()[0] + vtxtem2->point()[0]) / 8.0,
				3.0 * (v1->point()[1] + v2->point()[1]) / 8.0 + (vtxtem1->point()[1] + vtxtem2->point()[1]) / 8.0,
				3.0 * (v1->point()[2] + v2->point()[2]) / 8.0 + (vtxtem1->point()[2] + vtxtem2->point()[2]) / 8.0
			);
			nv->set_point(*np);
			//nv2->set_point(*np);
			newvtx.push_back(nv);

			std::pair<CVertex*, CVertex*> pair(v1, v2);
			vtxpair.push_back(pair);
			he1->set_updated(true);
			he1->he_sym()->set_updated(true);

		}

		if (mesh.isBoundary(he2) && !he2->updated()) {
			CVertex* nv = mesh.createVertex(mesh.max_vertex_id() + 1);
			//CVertex* nv2 = mesh2.createVertex(mesh2.max_vertex_id() + 1);

			CPoint* np = new CPoint(
				(v2->point()[0] + v0->point()[0]) / 2.0,
				(v2->point()[1] + v0->point()[1]) / 2.0,
				(v2->point()[2] + v0->point()[2]) / 2.0
			);
			nv->set_point(*np);
			//nv2->set_point(*np);
			newvtx.push_back(nv);


			std::pair<CVertex*, CVertex*> pair(v2, v0);
			vtxpair.push_back(pair);

			he0->set_updated(true);
		}
		else if (!he2->updated()) {
			CHalfEdge* he2_prev = he2->he_prev();
			CHalfEdge* he2_prev_opp = he2->he_sym()->he_prev();

			CVertex* vtxtem1 = he2_prev->source();
			CVertex* vtxtem2 = he2_prev_opp->source();

			CVertex* nv = mesh.createVertex(mesh.max_vertex_id() + 1);
			//CVertex* nv2 = mesh2.createVertex(mesh2.max_vertex_id() + 1);

			CPoint* np = new CPoint(
				3.0 * (v2->point()[0] + v0->point()[0]) / 8.0 + (vtxtem1->point()[0] + vtxtem2->point()[0]) / 8.0,
				3.0 * (v2->point()[1] + v0->point()[1]) / 8.0 + (vtxtem1->point()[1] + vtxtem2->point()[1]) / 8.0,
				3.0 * (v2->point()[2] + v0->point()[2]) / 8.0 + (vtxtem1->point()[2] + vtxtem2->point()[2]) / 8.0
			);
			nv->set_point(*np);
			//nv2->set_point(*np);
			newvtx.push_back(nv);

			std::pair<CVertex*, CVertex*> pair(v2, v0);
			vtxpair.push_back(pair);
			he2->set_updated(true);
			he2->he_sym()->set_updated(true);
		}

		//if ( he0->updated() || he1->updated() || he2->updated() ) 
	}


	for (std::vector<CFace*>::iterator it = faces.begin(); it != faces.end(); ++it) {
		CFace* item = *it;
		auto he0 = item->halfedge();
		auto v0 = he0->source();

		auto he1 = he0->he_next();
		auto v1 = he1->source();

		auto he2 = he1->he_next();
		auto v2 = he2->source();

		int index0 = getNewVertexIndexByOldVertexes(vtxpair, v0, v1, vertexNum);
		int index1 = getNewVertexIndexByOldVertexes(vtxpair, v1, v2, vertexNum);
		int index2 = getNewVertexIndexByOldVertexes(vtxpair, v2, v0, vertexNum);

		if (index0 != -1 && index1 != -1 && index2 != -1) {

			CVertex* vi0 = mesh.idVertex(index0);
			CVertex* vi1 = mesh.idVertex(index1);
			CVertex* vi2 = mesh.idVertex(index2);
			CVertex* face0[3] = { v0, vi0, vi2 };
			CVertex* face1[3] = { vi0, v1, vi1 };
			CVertex* face2[3] = { vi0, vi1, vi2 };
			
			//The last face seems to be causing issues.
			CVertex* face3[3] = { vi1, v2, vi2 };
			//mesh.createFace(face0, mesh.max_face_id() + 1);
			//mesh.createFace(face1, mesh.max_face_id() + 1);
			//mesh.createFace(face2, mesh.max_face_id() + 1);
			//mesh.createFace(face3, mesh.max_face_id() + 1);

			mesh2.createFace(face0, mesh2.max_face_id() + 1);
			mesh2.createFace(face1, mesh2.max_face_id() + 1);
			mesh2.createFace(face2, mesh2.max_face_id() + 1);
			mesh2.createFace(face3, mesh2.max_face_id() + 1);
			//mesh2.write_m(outmeshfile.c_str());
		}
	}

	//mesh2.write_m("testing.m");

	//Update vertices
	for (std::vector<CVertex*>::iterator it = newvtx.begin(); it != newvtx.end(); ++it) {
		CVertex* item = *it;
		std::vector<CVertex*> adj_verts;
		getAdjacentVertices(*it, &adj_verts);
		if (item->boundary()) {
			CPoint* np = new CPoint(
				(adj_verts[0]->point()[0] + adj_verts[1]->point()[0] + 6.0 * item->point()[0] ) / 8.0,
				(adj_verts[0]->point()[1] + adj_verts[1]->point()[1] + 6.0 * item->point()[1]) / 8.0,
				(adj_verts[0]->point()[2] + adj_verts[1]->point()[2] + 6.0 * item->point()[2]) / 8.0
			);
			item->set_point(*np);
		}
		else {
			int n = adj_verts.size();
			if (n == 3) {
				CPoint* np = new CPoint(
					3.0 * (adj_verts[0]->point()[0]
							+ adj_verts[1]->point()[0]
							+ adj_verts[2]->point()[0]
							+ 7.0 * item->point()[0]) / 16.0,
					3.0 * (adj_verts[0]->point()[1]
						+ adj_verts[1]->point()[1]
						+ adj_verts[2]->point()[1]
						+ 7.0 * item->point()[1]) / 16.0,
					3.0 * (adj_verts[0]->point()[2]
						+ adj_verts[1]->point()[2]
						+ adj_verts[2]->point()[2]
						+ 7.0 * item->point()[2]) / 16.0
				);
				item->set_point(*np);

			}
			else if (n > 3) {
				float beta = 3.0 / (8.0 *n);
				float xt, yt, zt = yt = xt = 0.0;
				for (int j = 0; j < n; ++j) {
					xt += beta * adj_verts[j]->point()[0];
					yt += beta * adj_verts[j]->point()[1];
					zt += beta * adj_verts[j]->point()[2];
				}
				double npx = (xt + 5.0 * item->point()[0]) / 8.0;
				double npy = (yt + 5.0 * item->point()[1]) / 8.0;
				double npz = (zt + 5.0 * item->point()[2]) / 8.0;

				npx = roundDouble(npx);
				npy = roundDouble(npy);
				npz = roundDouble(npz);
				

				CPoint* np = new CPoint(
					(xt + 5.0 * item->point()[0]) / 8.0,
					(yt + 5.0 * item->point()[1]) / 8.0,
					(zt + 5.0 * item->point()[2]) / 8.0
				);
				
				item->set_point(*np);
			}
		}
	}

	// remove old faces, and then generate a completely new mesh.


	for (std::vector<CVertex*>::iterator it = newvtx.begin(); it != newvtx.end(); ++it) {
		CVertex* item = *it;
		CVertex* nv = mesh2.createVertex(item->id());
		CPoint* np = &(item->point());
		nv->set_point(*np);
	}

	mesh2.write_m(outmeshfile.c_str());


	//std::list<CVertex*> mvs2 = mesh.vertices();
	//std::vector<CVertex*> verts2{ std::begin(mvs2), std::end(mvs2) };
	////bool readMesh = true;

	//for (std::vector<CVertex*>::iterator it = verts2.begin(); it != verts2.end(); ++it) {
	//	//if (!readMesh) mesh.read_m(meshfile.c_str());
	//	CVertex* item = *it;
	//	std::vector<CVertex*> adj_verts;
	//	getAdjacentVertices(*it, &adj_verts);

	//	if (item->boundary())
	//	{

	//		//new_vertexes.at(i).x = 
	//		//	(new_vertexes.at(adjVtx_ids[0]).x + new_vertexes.at(adjVtx_ids[1]).x + 6 * (new_vertexes.at(i).x)) / 8.0;
	//		//new_vertexes.at(i).y = 
	//		//	(new_vertexes.at(adjVtx_ids[0]).y + new_vertexes.at(adjVtx_ids[1]).y + 6 * (new_vertexes.at(i).y)) / 8.0;
	//		//new_vertexes.at(i).z = 
	//		//	(new_vertexes.at(adjVtx_ids[0]).z + new_vertexes.at(adjVtx_ids[1]).z + 6 * (new_vertexes.at(i).z)) / 8.0;
	//	}
	//	else {

	//	}

	//	//int n = adj_verts.size();

	//	//if (n == 3){

	//	//	CPoint *point;

	//	//}
	//	//else if (n > 3) {

	//	//	//float beta = 3 / (8.0*n);
	//	//	//float x_tem = 0.0;
	//	//	//float y_tem = 0.0;
	//	//	//float z_tem = 0.0;

	//	//	//for (int j = 0; j<n; j++) {
	//	//	//	x_tem += beta*new_vertexes.at(adjVtx_ids[j]).x;
	//	//	//	y_tem += beta*new_vertexes.at(adjVtx_ids[j]).y;
	//	//	//	z_tem += beta*new_vertexes.at(adjVtx_ids[j]).z;

	//	//	//}
	//	//	//new_vertexes.at(i).x = x_tem + 5 * new_vertexes.at(i).x / 8.0;
	//	//	//new_vertexes.at(i).y = y_tem + 5 * new_vertexes.at(i).y / 8.0;
	//	//	//new_vertexes.at(i).z = z_tem + 5 * new_vertexes.at(i).z / 8.0;
	//	//}
	//	//readmesh = false;
	//}
	return;
}
#endif