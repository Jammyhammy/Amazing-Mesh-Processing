#ifndef  _MESHLIB_VERTEX_H_
#define _MESHLIB_VERTEX_H_

#include <stdlib.h>
#include <string>
#include "../Geometry/Point.h"
#include "../Geometry/Point2.h"


namespace MeshLib {

class CHalfEdge;

class CVertex
{
    public:
        CVertex() {
            m_halfedge = NULL;
            m_boundary = false;
			m_alive = true;
        };
        ~CVertex() {};

        CPoint& point()    {
            return m_point;
        };
        CPoint& normal()   {
            return m_normal;
        };
        CPoint2& uv()       {
            return m_uv;
        };
		CPoint& point_updated(){
			return m_point_update;
		};
        CHalfEdge* most_ccw_out_halfedge();
        CHalfEdge* most_clw_out_halfedge();
        CHalfEdge* most_ccw_in_halfedge();
        CHalfEdge* most_clw_in_halfedge();

        CHalfEdge*& halfedge() {
            return m_halfedge;
        };

        std::string& string() {
            return m_string;
        };

        int&   id() {
            return m_id;
        };
        bool& boundary() {
            return m_boundary;
        };
		bool& status() {
			return m_alive;
		};

		void set_id(int id){
			m_id = id;
		}

		void set_point(CPoint p)
		{
			m_point = p;
		}

		void set_point_update(CPoint update)
		{
			m_point_update = update;
		}

		//how many edges a vertex is shared with
		int num_neighbors()
		{

		}
		bool operator<(const CVertex& m_vertex)
		{
			return (m_vertex.m_id < m_id);
		}
		bool operator==(const CVertex& m_vertex)
		{
			return (m_vertex.m_id == m_id);
		}
		


    private:
        int  m_id ;

        CPoint m_point;
        CPoint m_normal;
		CPoint m_point_update;
        CPoint2 m_uv;

        CHalfEdge* m_halfedge;
        bool              m_boundary;
		bool		m_alive;
        std::string     m_string;
      

}; //class CVertex

}//name space MeshLib

#endif //_MESHLIB_VERTEX_H_defined