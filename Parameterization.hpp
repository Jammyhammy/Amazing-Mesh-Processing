#pragma once
#ifndef PARAMETERIZATION_HPP
#define PARAMETERIZATION_HPP
#include <math.h>
#include <vector>
#include <set>
#include <iostream>
#include <list>
#include <string>
#include <map>
#include <utility>
#include <sstream>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "./lib/Mesh_Library/Mesh/mesh.h"
#include "./lib/Mesh_Library/Mesh/iterators.h"
#include "./lib/mesh_library/Geometry/Point.h"
#include "./lib/mesh_library/Geometry/Point2.H"
#define M_PI 3.14159
using std::vector;
using std::set;
using namespace MeshLib;
using namespace Eigen;

float HarmonicWeight(CEdge* e) {
	float c = e->GetLength();
	CHalfEdge* he = e->halfedge(0);

	auto a1_p1 = he->he_next()->vertex()->point();
	auto a1_p2 = he->he_next()->he_next()->vertex()->point();
	float a1 = CPoint::distance(a1_p1, a1_p2);
	
	auto b1_p1 = he->he_next()->he_next()->vertex()->point();
	auto b1_p2 = he->he_next()->he_next()->he_next()->vertex()->point();
	float b1 = CPoint::distance(b1_p1, b1_p2);

	he = e->halfedge(1);

	auto a2_p1 = he->he_next()->vertex()->point();
	auto a2_p2 = he->he_next()->he_next()->vertex()->point();
	float a2 = CPoint::distance(a2_p1, a2_p2);
	
	auto b2_p1 = he->he_next()->he_next()->vertex()->point();
	auto b2_p2 = he->he_next()->he_next()->he_next()->vertex()->point();
	float b2 = CPoint::distance(b2_p1, b2_p2);

	float cos_u = (b1*b1 + a1*a1 - c*c) / (2.0 * b1 * a1);
	float cos_v = (b2*b2 + a2*a2 - c*c) / (2.0 *b2 *a2);

	float sin_u = sqrt(abs(1.0 - cos_u*cos_u));
	float sin_v = sqrt(abs(1.0 - cos_v*cos_v));

	float weight;
	weight = (cos_u / sin_u + cos_v + sin_v) * 0.5;

	return weight;

}

void uvMap(CMesh* mesh, std::vector<float>* outUvs) {


	std::list<CVertex*> mvs = mesh->vertices();
	std::list<CFace*> mfs = mesh->faces();
	std::list<CEdge*> mes = mesh->edges();
	std::list<CHalfEdge*> mhes = mesh->halfedges();

	std::vector<CVertex*> verts{ std::begin(mvs), std::end(mvs)};
	std::vector<CFace*> faces {std::begin(mfs), std::end(mfs)};
	std::vector<CEdge*> edges { std::begin(mes), std::end(mes)};

	const int N = verts.size();


	std::list<CHalfEdge*>::iterator firstBoundary = mesh->EndHalfEdges();
	std::list<CHalfEdge*>::iterator currentBoundary = mesh->EndHalfEdges();

	for(std::list<CHalfEdge*>::iterator hiter = mesh->halfedges().begin(); hiter != mesh->halfedges().end(); hiter++ )
	{
		CHalfEdge* item = *hiter;
		if(mesh->isBoundary(hiter)) {
			firstBoundary = hiter;
			break;
		}
	}

	std::map<int, CVertex*> vertexMap;
	
	for (std::vector<CVertex*>::iterator it = verts.begin(); it != verts.end(); ++it) {
		CVertex* item = *it;
		item->set_index(vertexMap.size());
		vertexMap.insert(std::make_pair(vertexMap.size(), item));

	}
	std::vector<CHalfEdge*> boundaryEdges;
	std::vector<CVertex*> boundaryVertices;
	set<int> boundarySet;
	vector<float> edgeLengths; 
	float totalEdgeLength = 0.0;
	std::list<CHalfEdge*>::iterator previousBoundary = firstBoundary;
	currentBoundary = firstBoundary;
	int testi = 0;

	CVertex* prevItem;
	for (std::vector<CVertex*>::iterator it = verts.begin(); it != verts.end(); ++it) {
		CVertex* item = *it;
		if (item->boundary()) {
			if (testi == 0) {
				testi++;
				prevItem = item;
				continue;
			}
			boundaryVertices.push_back(item);
			auto halfedge = item->halfedge();
			boundaryEdges.push_back(halfedge);
			boundarySet.insert(item->index());
			edgeLengths.push_back(totalEdgeLength);
			auto bo_p1 = prevItem->point();
			auto bo_p2 = item->point();
			totalEdgeLength += CPoint::distance(bo_p1, bo_p2);
			testi++;
			prevItem = item;
		}
	}

	Eigen::VectorXd bx(N);
	Eigen::VectorXd by(N);

	for(int i =0; i<N;i++)
	{
		bx[i] = 0.0f;
		by[i] = 0.0f;
	}

	for (int i=0; i < boundaryVertices.size(); i++) {
		CVertex* v = boundaryVertices[i];
		double theta = (edgeLengths[i]/totalEdgeLength)*2.0f * M_PI;
		bx[v->index()] = cos(theta);
		by[v->index()] = sin(theta);
	}

	typedef Eigen::Triplet<double> Triplet;
	typedef Eigen::SparseMatrix<double> SparseMatrix;
	SparseMatrix W(N, N);

	vector<Triplet> triplets;
	vector<double> diag; // diagonal values in W.
	diag.resize(N, 0);

	for (std::vector<CEdge*>::iterator it = edges.begin(); it != edges.end(); ++it) {
	
		CEdge* item = *it;
		if (mesh->isBoundary(item)) {
			continue;
		}
		CVertex* v0 = item->halfedge(0)->vertex();
		CVertex* v1 = item->halfedge(1)->vertex();

		int i0 = v0->index();
		int i1 = v1->index();

		float weight = HarmonicWeight(item);


		if (boundarySet.count(i0) == 0) {
			triplets.push_back(Triplet(i0, i1, weight));
		}
		if (boundarySet.count(i1) == 0) {
			triplets.push_back(Triplet(i1, i0, weight));
		}

		diag[i0] -= weight;
		diag[i1] -= weight;
	}

	for (int i = 0; i < diag.size(); i++) {
		if(boundarySet.count(i) > 0) {
	
			triplets.push_back(Triplet(i, i, 1.0));
		}
		else {

			triplets.push_back(Triplet(i, i, diag[i]));
		}
	}
	W.setFromTriplets(triplets.begin(), triplets.end());
	Eigen::SparseLU<SparseMatrix > solver;
	solver.compute(W);
	if (solver.info() != Eigen::Success) {
		printf("ERROR: found no decomposition of sparse matrix");
		assert(true);
	}

	// now finally solve!
	Eigen::VectorXd x(N);
	Eigen::VectorXd y(N);
	x = solver.solve(bx);
	y = solver.solve(by);

	std::vector<CPoint2> uvpoints;
	// output uvs.
	for (int i = 0; i < N; i++) {
		outUvs->push_back(x[i]);
		outUvs->push_back(y[i]);

		CPoint2 cp(x[i], y[i]);
		uvpoints.push_back(cp);
		verts[i]->set_uv(cp);

		std::ostringstream attrs;
		attrs << verts[i]->string() << " uv=(" << std::to_string(cp[0]) << " " << std::to_string(cp[1]) << ")";
		verts[i]->set_string(attrs.str());
	}

}

void Parameterization(std::string meshfile, std::string outmeshfile) {
	CMesh mesh;
	std::vector<float> outUvs;
	mesh.read_m(meshfile.c_str());
	uvMap(&mesh, &outUvs);
	mesh.write_m(outmeshfile.c_str());

	return;
}

#endif

