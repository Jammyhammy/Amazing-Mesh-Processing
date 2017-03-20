#pragma once
#ifndef FILLHOLE_HPP
#define FILLHOLE_HPP

#define HH_CAT(a, b) a ## b
#define HH_CAT2(a, b) HH_CAT(a, b)
#define HH_ID(x) HH_CAT(_hh_id_, x) // private identifier in a macro definition
#define HH_UNIQUE_ID(x) HH_CAT2(HH_CAT2(HH_CAT2(_hh_id_, __COUNTER__), _), x)
#define for_range_aux(T, i, lb, ub, u) for (T u = ub, i = lb; i<u; i++)
#define for_range(T, i, lb, ub) for_range_aux(T, i, lb, ub, HH_UNIQUE_ID(u))
#define for_int(i, ub)      for_range(int, i, 0, ub)

#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>

#include "./lib/Mesh_Library/Mesh/mesh.h"
#include "./lib/Mesh_Library/Mesh/iterators.h"
#include "./lib/mesh_library/Geometry/Point.h"

using namespace MeshLib;
using namespace std;


void inline adjacentVertices(CVertex *v, std::vector<CVertex*> *adj_v)
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

double inline triangleArea(CVertex* v0, CVertex* v1, CVertex* v2) {
	CPoint p0 = v0->point();
	CPoint p1 = v1->point();
	CPoint p2 = v2->point();
	return (((p1 - p0) ^ (p2 - p1)).norm() * 0.5);
}

void inline getBoundaryVertices(vector<CVertex*> *v, CMesh *mesh)
{
	std::list<CVertex*> mvs = mesh->vertices();
	std::vector<CVertex*> verts{ std::begin(mvs), std::end(mvs) };

	for (std::vector<CVertex*>::iterator it = verts.begin(); it != verts.end(); ++it) {
		CVertex* item = *it;
		if (item->boundary()) {
			v->push_back(item);
		}
	}
}

void inline findHole(std::vector<CVertex*> *h, std::vector<CVertex*> *v) {
	while (v->size()) {
		CVertex* bv = v->at(0);
		int sid = bv->id();
		CHalfEdge *he = bv->halfedge();
		h->push_back(bv);
		bv = bv->most_ccw_in_halfedge()->source();
		while (bv->id() != sid) {
			h->push_back(bv);
			bv = bv->most_ccw_in_halfedge()->source();
		}
		for (std::vector<CVertex*>::iterator it = h->begin(); it != h->end(); ++it) {
			CVertex* item = *it;
			std::vector<CVertex*>::iterator res = std::find(v->begin(), v->end(), item);
			if (res != v->end())
				v->erase(res);
		}
	}
}

void inline calcMinWeight(int i, int k, double *weight, std::vector<CVertex*> *hole, CMesh* mesh, int *index)
{
	double minW = 100000.0;
	int res;
	int num = hole->size();
	CVertex* v[3];
	for (int m = i + 1; m < k; m++) {
		double tmp = 0;

		tmp += weight[i*num + m];

		tmp += weight[m*num + k];
		v[0] = hole->at(i);
		v[1] = hole->at(m);
		v[2] = hole->at(k);
		tmp += triangleArea(v[0], v[1], v[2]);
		if (tmp < minW) {
			minW = tmp;
			res = m;
		}
	}
	weight[i* num + k] = minW;
	index[i * num + k] = res;
}

CPoint inline vertexNormal(CVertex* v0, CVertex* v1, CVertex* v2) {
	auto zomg = (v1->point() - v0->point());
	auto zomg2 = (v2->point() - v1->point());
	auto normal = zomg ^ zomg2;
	auto normalized = normal / normal.norm();
	return normalized;
}

CPoint inline pointNormal(CPoint* v0, CPoint* v1, CPoint* v2) {
	auto normal = (*v1 - *v0) ^ (*v2 - *v1);
	auto normalized = normal / normal.norm();
	return normal;
}

double inline angleVectors(CPoint* p0, CPoint* p1) {
	auto in = (*p0 * *p1) / (p0->norm() * p1->norm());
	auto theta = acos(in);
	return theta;
}

CHalfEdge* queryhedgecheck(CVertex* v1, CVertex* v2) {
	std::vector<CHalfEdge*> arhes;
	std::vector<CVertex*> adj_v;
	bool foundIssue = false;

	auto he = v1->halfedge();
	arhes.push_back(he);

	he->set_updated(true);
	he = he->he_next();
	do {
		if (he->vertex() == v1) {
			arhes.push_back(he);
		}
		//arhes.push_back(he);

		he->set_updated(true);
		he = he->he_sym();
		if (!he) {
			foundIssue = true;
			break;
		}
		if (he->vertex() == v1) {
			arhes.push_back(he);
		}
		//arhes.push_back(he);
		he->set_updated(true);
		he = he->he_next();
		if (!he) {
			foundIssue = true;
			break;
		}
	} while (!(he->updated()));

	if (foundIssue) {
		he = v1->halfedge()->he_prev();
		if (he && !(he->updated())) {
			do {

				if (he->vertex() == v1) {
					arhes.push_back(he);
				}
				//arhes.push_back(he);

				he->set_updated(true);
				he = he->he_sym();
				if (!he) {
					break;
				}
				if (he->vertex() == v1) {
					arhes.push_back(he);
				}
				//arhes.push_back(he);

				he->set_updated(true);
				he = he->he_prev();
				if (!he) {
					break;
				}
			} while (!(he->updated()));
		}


	}
	//adjacentVertices(v1, &adj_v);

	for (std::vector<CHalfEdge*>::iterator it = arhes.begin(); it != arhes.end(); ++it) {
		CHalfEdge* item = *it;
		if (item->vertex() == v2) {
			return item;
		}
	}

	return nullptr;
}

bool inline legalcreateface(CFace* f) {
	CVertex* v[3];
	v[0] = f->halfedge()->he_prev()->vertex();
	v[1] = f->halfedge()->vertex();
	v[2] = f->halfedge()->he_next()->vertex();

	CVertex* v0 = v[2];
	for_int(i, 3) {
		if (queryhedgecheck(v0, v[i])) return false;
		v0 = v[i];
	}
	return true;
}

void inline triangulate(int i, int k, int *index, double *weight, std::vector<CVertex*> *hole, int bid, CMesh *mesh) {
	static int cnt = 0;
	CVertex* v[3];
	int num = hole->size();
	if (i + 2 == k) {
		v[0] = hole->at(i);
		v[1] = hole->at(i + 1);
		v[2] = hole->at(k);
		mesh->createFace(v, bid + cnt);
		cnt++;
	}
	else {
		int id = index[i*num + k];
		double w = weight[i* num + k];

		if (id != i + 1) triangulate(i, id, index, weight, hole, bid, mesh);
		v[0] = hole->at(i);
		v[1] = hole->at(id);
		v[2] = hole->at(k);
		auto vtest = hole->at(i);
		auto vyes1 = mesh->idVertex(vtest->id());
		auto vyes2 = vyes1->halfedge()->he_next()->target();
		auto vyes3 = vyes2->halfedge()->he_next()->target();
		//auto vyes1 = mesh->idVertex(vtest->id());
		//auto vyes2 = vyes1->halfedge()->he_next()->target();
		//auto vyes3 = vyes1->halfedge()->he_next()->he_next()->target();


		auto vc = vertexNormal(vyes1, vyes2, vyes3);

		//auto vc = vertexNormal(v1, v2, v3);
		//auto vcn = vc / vc.norm();
		auto pc = vertexNormal(v[0], v[1], v[2]);
		//auto pcn = pc / pc.norm();
		auto pcv = v[2]->point() + v[0]->point();
		auto ndv = pcv * pc;
		//auto ndvn = pcv * pcn;
		//auto vpc = vertexNormal(vtest, v[1], v[2]);
		////auto vpcn = vpc / vpc.norm();
		//auto vpcv = v[2]->point() + v[0]->point();
		//auto vndv = vpcv * vpc;
		//auto vndvn = vpcv * vpcn;
		auto testdot = vc * pc;
		//auto testv = pointNormal(&vpcv, &(v[1]->point()), &(v[2]->point()));
		//auto moretest = vtes
		auto testangles = angleVectors(&vc, &pc);
		if (testangles < 1.32 && ndv > 0.0) {
			if (bid + cnt == 69199)
				auto test2 = "wa";
			mesh->createFace(v, bid + cnt);
		}
		//if (testdot < 0.0) {
		//	auto test = "test";
		//}

		//if (ndv < 0.0)
		//{
		//	auto test = "test";
		//}
		//if (bid + cnt == 69198 || bid + cnt == 69199 || bid + cnt == 6001 || bid + cnt == 6002 || bid + cnt == 6003 || bid + cnt == 6004 
		//	|| bid + cnt >= 6005)
		//{
		//	CVertex* dv[3];
		//	auto test = "test";
		//	testangles = angleVectors(&vc, &pc);
		//	testdot = vc * pc;
		//	//if (bid + cnt == 6001) {
		//	dv[0] = v[0];
		//	dv[1] = v[2];
		//	dv[2] = v[1];
		//	//}
		//	pc = vertexNormal(dv[0], dv[1], dv[2]);
		//	testangles = angleVectors(&vc, &pc);
		//	testdot = vc * pc;
		//	if (testdot >= 0.0) {
		//		//if (bid + cnt == 6006) 
		//		//	triangulate(i, id, index, weight, hole, bid, mesh);
		//		//	//auto test4 = "test";
		//		//else {
		//		//	mesh->createFace(dv, bid + cnt);
		//		//}
		//	}
		//	else {
		//		CVertex* dvv[3];
		//		dv[0] = v[1];
		//		dv[1] = v[2];
		//		dv[2] = v[0];
		//		if (bid + cnt == 6007 || bid + cnt == 6006) 
		//			auto test4 = "test";
		//		pc = vertexNormal(dv[0], dv[1], dv[2]);
		//		testangles = angleVectors(&vc, &pc);
		//		testdot = vc * pc;
		//		auto test3 = "test";
		//		//if (testdot >= 0.0) {
		//		//	mesh->createFace(dv, bid + cnt);
		//		//}

		//	}


		//	auto test2 = "test";

		//	//auto testface = mesh->createFace(v, bid + cnt);
		//	//auto legal = legalcreateface(testface);
		//}

		cnt++;

		if (id != k - 1)
			triangulate(id, k, index, weight, hole, bid, mesh);
	}
}

void inline addFace(int *index, double *w, std::vector<CVertex*> *hole, CMesh *mesh) {
	int base_id = mesh->max_face_id() + 1;
	triangulate(0, hole->size() - 1, index, w, hole, base_id, mesh);
}

void inline faceFillHole(std::vector<CVertex*> *hole, CMesh *mesh) {
	int num = hole->size();
	CVertex* v[3];
	if (num == 3) {
		v[0] = hole->at(0);
		v[1] = hole->at(1);
		v[2] = hole->at(2);
		mesh->createFace(v, mesh->max_face_id() + 1);
		return;
	}
	int *index = (int *)malloc(num * num * sizeof(int));
	for (int i = 0; i < num * num; i++)
		index[i] = -1;
	double *w = (double *)malloc(num * num * sizeof(double));
	for (int i = 0; i < num; i++) {
		for (int j = 0; j < num; j++) {
			w[i * num + j] = -1;
		}
	}
	for (int i = 0; i <= num - 2; i++)
		w[i * num + i + 1] = 0;

	for (int i = 0; i <= num - 3; i++) {
		v[0] = hole->at(i);
		v[1] = hole->at(i + 1);
		v[2] = hole->at(i + 2);
		w[i * num + i + 2] = triangleArea(v[0], v[1], v[2]);
	}
	int cnt = 1;
	for (int j = 2; j <= num - 1; j++) {
		for (int i = 0; i <= num - j - 1; i++) {
			int k = i + j;
			calcMinWeight(i, k, w, hole, mesh, index);
		}
	}
	addFace(index, w, hole, mesh);
	free(w);
	free(index);
	index = NULL;
	w = NULL;
}


double inline vertexAngle(CHalfEdge *he, CMesh *mesh)
{
	double l1 = mesh->edgeLength(he->edge());
	double l2 = mesh->edgeLength(he->he_next()->edge());
	double l3 = mesh->edgeLength(he->he_next()->he_next()->edge());
	double cos = (l1*l1 + l2*l2 - l3*l3) / (2 * l1 * l2);
	//double cos = (l1*l1 + l2*l2 - l3*l3) / (2.0 * l1 * l2);
	return acos(cos);
}

double inline sumVertexAngle(CVertex *v, CMesh *mesh)
{
	CHalfEdge *he = v->halfedge();
	double va = vertexAngle(he, mesh);
	CHalfEdge *ne = he->clw_rotate_about_target();

	while (ne != NULL && !(ne->edge() == he->edge())) {
		va += vertexAngle(ne, mesh);

		ne = ne->clw_rotate_about_target();
	}

	return va;
}

void inline boundary_edge_update(CMesh *mesh)
{
	std::list<CEdge*> m_edges = mesh->edges();

	for (std::list<CEdge*>::iterator eiter = m_edges.begin(); eiter != m_edges.end(); ++eiter)
	{
		CEdge *edge = *eiter;
		CHalfEdge *he[2];

		he[0] = edge->halfedge(0);
		he[1] = edge->halfedge(1);

		assert(he[0] != NULL);
		
		if (he[1] != NULL)
		{
			//assert(he[0]->target() == he[1]->source() && he[0]->source() == he[1]->target());

			if (he[0]->target()->id() < he[0]->source()->id())
			{
				edge->halfedge(0) = he[1];
				edge->halfedge(1) = he[0];
			}

			assert(mesh->edgeVertex1(edge)->id() < mesh->edgeVertex2(edge)->id());

			he[0]->vertex()->boundary() = false;
			he[0]->he_prev()->vertex()->boundary() = false;
		}
		else
		{
			he[0]->vertex()->boundary() = true;
			he[0]->he_prev()->vertex()->boundary() = true;
		}

	}
}

void ccw_halfedge_update(CMesh *mesh, bool* loadAgain)
{
	//Arrange the boundary half_edge of boundary vertices, to make its halfedge
	//to be the most ccw in half_edge
	std::list<CVertex*> m_verts = mesh->vertices();


	for (std::list<CVertex*>::iterator viter = m_verts.begin(); viter != m_verts.end(); ++viter)
	{
		CVertex *v = *viter;
		if (!v->boundary())
			continue;

		CHalfEdge *he = v->halfedge();
		int max_tries = 25;

		while (he->he_sym() != NULL)
		{
			if (max_tries <= 0) {
				he->vertex()->boundary() = true;
				he->he_prev()->vertex()->boundary() = true;
				*loadAgain = false;
				break;
			}
			he = he->ccw_rotate_about_target();
			max_tries--;
		}
		v->halfedge() = he;
	}
}

void inline holeFill(CMesh* mesh, std::string outmeshfile, bool* loadAgain) {
	std::vector<CVertex*> boundaryVertex;
	getBoundaryVertices(&boundaryVertex, mesh);
	while (!boundaryVertex.empty()) {
		std::vector<CVertex*> hole;
		findHole(&hole, &boundaryVertex);
		faceFillHole(&hole, mesh);
	}
	boundary_edge_update(mesh);
	ccw_halfedge_update(mesh, loadAgain);
}


void fillHoles(std::string meshfile, std::string outmeshfile)
{
	CMesh mesh;
	bool loadAgain = true;
	mesh.read_m(meshfile.c_str());
	holeFill(&mesh, outmeshfile, &loadAgain);
	mesh.write_m(outmeshfile.c_str());
	for (int i = 0; i < 4; i++) {
		if (!loadAgain) break;
		CMesh mesh2;
		std::vector<CVertex*> boundaryVertex;
		mesh2.read_m(outmeshfile.c_str());
		getBoundaryVertices(&boundaryVertex, &mesh2);
		if (!boundaryVertex.empty())
		{
			holeFill(&mesh2, outmeshfile, &loadAgain);
			mesh2.write_m(outmeshfile.c_str());
		}
	}
}

#endif