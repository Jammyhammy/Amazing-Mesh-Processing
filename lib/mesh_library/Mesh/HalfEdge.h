#ifndef _MESHLIB_HALFEDGE_H_
#define _MESHLIB_HALFEDGE_H_

#include  <assert.h>
#include <cmath>
#include <string>
#include "../Mesh/Edge.h"
#include "../Geometry/Point.h"



namespace MeshLib {

class CVertex;
class CEdge;
class CFace;

class CHalfEdge
{
    public:

        CHalfEdge() {
            m_edge = NULL;
            m_vertex = NULL;
            m_prev = NULL;
            m_next = NULL;
            m_face = NULL;
			m_alive = true;
            m_updated = false;
        };

        CHalfEdge(int id) {
            m_edge = NULL;
            m_vertex = NULL;
            m_prev = NULL;
            m_next = NULL;
            m_face = NULL;
			m_alive = true;
            m_updated = false;
            m_id = id;
        };        

        ~CHalfEdge() { };

        CEdge*&            edge()    {
            return m_edge;
        };
        CVertex*&          vertex()  {
            return m_vertex;
        };
        CVertex*&          target()  {
            return m_vertex;
        };
        CVertex*&          source()  {
            return m_prev->vertex();
        };
        CHalfEdge*&     he_prev() {
            return m_prev;
        };
        CHalfEdge*&     he_next() {
            return m_next;
        };

        CHalfEdge*& he_sym()  {
            return m_edge->other(this);
        };

        CFace*&      face()    {
            return m_face;
        };

		bool &	status() {
			return m_alive;
		};

		bool &	updated() {
			return m_updated;
		};

        void set_updated(bool updated) {
            m_updated = updated;
        }

        CHalfEdge*    ccw_rotate_about_target();
        CHalfEdge*    clw_rotate_about_target();

        CHalfEdge*    ccw_rotate_about_source();
        CHalfEdge*    clw_rotate_about_source();

        std::string& string() {
            return m_string;
        };

        int id() {
            return m_id;
        };
        
		double GetLength();
        
        
    private:

        CEdge*            m_edge;
        CFace*             m_face;
        CVertex*          m_vertex;     //target vertex
        CHalfEdge*     m_prev;
        CHalfEdge*      m_next;
		bool m_alive;
        bool m_updated;
        std::string          m_string;
        int m_id;
        
};


}//namespace MeshLib

#endif //_MESHLIB_HALFEDGE_H_ defined