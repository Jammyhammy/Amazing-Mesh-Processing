#include <assert.h>
#include "Edge.h"
#include "HalfEdge.h"
#include "Vertex.h"

using namespace MeshLib;

CEdgeKey::CEdgeKey(CVertex* v1, CVertex* v2)
{
    assert(v1->id() != v2->id());

    if(v1->id() < v2->id())
    {
        m_verts[0] = v1;
        m_verts[1] = v2;
    }
    else
    {
        m_verts[0] = v2;
        m_verts[1] = v1;
    }
}

bool CEdgeKey::operator<(const CEdgeKey& key) const
{
    if(m_verts[0]->id() < key.m_verts[0]->id()) return true;
    if(m_verts[0]->id() > key.m_verts[0]->id()) return false;

    if(m_verts[1]->id() < key.m_verts[1]->id()) return true;
    if(m_verts[1]->id() > key.m_verts[1]->id()) return false;

    return false;
}

bool CEdgeKey::operator==(const CEdgeKey& key) const
{
    if(m_verts[0]->id() != key.m_verts[0]->id()) return false;
    if(m_verts[1]->id() != key.m_verts[1]->id()) return false;

    return true;
}

double MeshLib::CEdge::GetLength()
{
	auto p1 = m_halfedge[0]->vertex()->point();
	auto p2 = m_halfedge[1]->vertex()->point();
	auto v = CPoint::distance(p1, p2);
	return v;
}
