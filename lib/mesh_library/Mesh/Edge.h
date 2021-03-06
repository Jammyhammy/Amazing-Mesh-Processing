#ifndef _MESHLIB_EDGE_H_
#define _MESHLIB_EDGE_H_

#include <assert.h>
#include <stdlib.h>
#include <cmath>
#include <string>
#include "../Geometry/Point.h"


namespace MeshLib {

class CHalfEdge;
class CVertex;

class CEdge
{
    public:
        CEdge() {
            m_halfedge[0] = NULL;
            m_halfedge[1] = NULL;
			m_alive = true;
        };
        ~CEdge() {};
		
		//get a half edge
        CHalfEdge*& halfedge(int i) {
            assert(0<=i && i < 2);
            return m_halfedge[i];
        };

		//check if an edge is a boundary edge or not
        bool         boundary() {
            return (m_halfedge[0] == NULL && m_halfedge[1] != NULL) || (m_halfedge[0] != NULL && m_halfedge[1] == NULL);
        };

		bool & status() {
			return m_alive;
		};

		//get the other half edge
        CHalfEdge*& other(CHalfEdge* he) {
            return (he != m_halfedge[0])?m_halfedge[0]:m_halfedge[1];
        };

        std::string& string() {
            return m_string;
        };

        double GetLength();

    private:

        CHalfEdge*   m_halfedge[2];
		bool m_alive;
        std::string      m_string;
        
};

class CEdgeKey
{
    public:
        CEdgeKey(CVertex* v1, CVertex* v2);
        ~CEdgeKey() {};
        bool operator< (const CEdgeKey& key) const;
        bool operator==(const CEdgeKey& key) const;

    protected:
        CVertex* m_verts[2];
};

}//name space MeshLib

#endif //_MESHLIB_EDGE_H_ defined