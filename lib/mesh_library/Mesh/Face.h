#ifndef _MESHLIB_FACE_H_
#define _MESHLIB_FACE_H_

#include <assert.h>
#include <string>


namespace MeshLib {


class CHalfEdge;

class CFace
{
    public:

        CFace() {
            m_halfedge = NULL;
			m_alive = true;
        };
        ~CFace() {};

        CHalfEdge*&     halfedge() {
            return m_halfedge;
        };
        int&                  id()          {
            return m_id;
        };
		bool &	status() {
			return m_alive;
		};
        const int           id() const {
            return m_id;
        };
        std::string&     string()     {
            return m_string;
        };
        
		void set_id(int id){
			m_id = id;
		}

    private:
        int                       m_id;
        CHalfEdge*     m_halfedge;
		bool m_alive;
        std::string        m_string;
        
};


}//name space MeshLib

#endif //_MESHLIB_FACE_H_ defined