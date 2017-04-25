#include "HalfEdge.h"
#include "Edge.h"
#include "Vertex.h"

using namespace MeshLib;


CHalfEdge* CHalfEdge::ccw_rotate_about_target()
{
    CHalfEdge* he_dual = he_sym();
    if(he_dual == NULL) return NULL;

    return he_dual->he_prev();
};

CHalfEdge* CHalfEdge::clw_rotate_about_target()
{
    CHalfEdge* he = he_next()->he_sym();
    return he;
};

CHalfEdge* CHalfEdge::ccw_rotate_about_source()
{

    CHalfEdge* he = he_prev()->he_sym();
    return he;
};

CHalfEdge* CHalfEdge::clw_rotate_about_source()
{
    CHalfEdge* he = he_sym();
    if(he == NULL) return NULL;
    return he->he_next();
}

double MeshLib::CHalfEdge::GetLength()
{
	auto p1 = m_vertex->point();
	auto p2 = m_edge->other(this)->vertex()->point();
	auto d = CPoint::distance(p1, p2);
	return d;
}
;
