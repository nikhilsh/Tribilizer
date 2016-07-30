#ifndef MESH_H
#define MESH_H

#include <map>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "CompFab.h"

namespace CompFab
{
	template <class T> inline T _max(T a, T b) { return a > b ? a : b; }
	template <class T> inline T _min(T a, T b) { return a < b ? a : b; }
	
	struct Mesh
	{
		std::vector<V3> v;
		std::vector<V3> n;
		std::vector<Vec2> tex;
		std::vector<V3i> texId;
		std::vector<V3i> t;
		
		inline Mesh() {};
		
		inline Mesh(const std::vector<V3> &_v,
					const std::vector<V3i> &_t):
		v(_v), t(_t) {}
		

		
		void rescale()
		{
			if (v.size() == 0) {
				std::cout << "empty mesh\n";
				return;
			}
			
			CompFab::V3 mn = v[0], mx = v[0];
			
			for(unsigned ii = 0; ii < v.size(); ++ii) {
				for (unsigned dim = 0; dim < 3; ++dim) {
					mn[dim] = _min(v[ii][dim], mn[dim]);
					mx[dim] = _max(v[ii][dim], mx[dim]);
				}
			}

			for(unsigned ii = 0; ii<v.size(); ++ii) {
				for (unsigned dim = 0; dim < 3; ++dim) {
					v[ii][dim] = v[ii][dim] - mn[dim];
				}
			}
		
			Real scale = 1 / (mx[0] - mn[0]);
			
			for(unsigned dim = 1; dim < 3; ++dim) {
				scale = _min(1 / (mx[dim] - mn[dim]), scale);
			}
			
			for(unsigned ii = 0; ii < v.size(); ++ii) {
				for (unsigned dim = 0; dim < 3; ++dim) {
					v[ii][dim] = v[ii][dim] * scale;
				}
			}
		}
		
		void computeNorms()
		{
			CompFab::V3 ZERO;
			n.resize(v.size(), ZERO);
			for(unsigned ii = 0; ii < t.size(); ++ii) {
				V3 a = v[t[ii][1]] - v[t[ii][0]];
				V3 b = v[t[ii][2]] - v[t[ii][0]];
                b = a.cross(b).normalize();
				for (unsigned jj = 0; jj < 3; ++jj) {
					n[t[ii][jj]] += b;
					if(t[ii][jj] >= (unsigned) n.size()) {
						std::cout << ii << " " << jj << " " << t[ii][jj] << " normal computation error\n";
					}
				}
			}
			for(unsigned int ii=0; ii<v.size(); ii++) {
				n[ii] = n[ii].normalize();
			}
		}
		
		Mesh(const CompFab::V3 *_v, const CompFab::V3i *_t)
		{
			v.assign(_v, _v + 8);
			t.assign(_t, _t + 12);
			
			computeNorms();
		}
		
		Mesh(const char * filename, bool normalize)
		{
			loadMesh(filename, normalize);
		}
		
		void loadMesh(const char * filename, bool normalize = true)
		{
			std::ifstream f;
			f.open(filename);
			if(!f.good()) {
				std::cout << "Error: cannot open mesh " << filename << "\n";
				return;
			}
			switch(filename[strlen(filename)-1]) {
				case 'y':
					readPly(f);
					break;
				case 'j':
					readObj(f);
					break;
				default:
					break;
			}
			if (normalize) {
				rescale();
			}
			
			computeNorms();
			
			f.close();
		}
		
		void save(const char * filename)
		{
			std::ofstream out;
			out.open(filename);
			save(out);
			out.close();
		}
		
		void save(std::ostream &out, std::vector<CompFab::V3> * vert=NULL)
		{
			std::string vTok("v");
			std::string fTok("f");
			std::string texTok("vt");
			char bslash='/';
			std::string tok;
			if(vert == NULL){
				vert = &v;
			}
			//Grid structure for Voxels
			for (size_t ii = 0; ii < vert->size(); ii++){
				out << vTok << " " << (*vert)[ii][0] << " " << (*vert)[ii][1] << " " << (*vert)[ii][2] << "\n";
			}
			if (tex.size() > 0){
				for (size_t ii = 0; ii < tex.size(); ++ii){
					out<<texTok << " " << tex[ii][0] << " " << tex[ii][1] << "\n";
				}
				for (size_t ii = 0; ii < t.size(); ++ii){
					out << fTok << " " << t[ii][0]+1 << bslash << texId[ii][0]+1 << " "
					<< t[ii][1]+1 << bslash << texId[ii][1]+1 << " "
					<< t[ii][2]+1 << bslash << texId[ii][2]+1 << "\n";
				}
			}else{
				for(size_t ii = 0; ii < t.size(); ++ii){
					out << fTok << " " << t[ii][0]+1 << " " <<
					t[ii][1]+1 << " " << t[ii][2]+1 << "\n";
				}
			}
			out << "#end\n";
		}

		void load(std::istream &in) { readObj(in); }
		
		void readPly(std::istream & f)
		{
			std::string line;
			std::string vertLine("element vertex");
			std::string faceLine("element face");
			std::string texLine("property float s");
			std::string endHeaderLine("end_header");

			while (true) {
				std::getline(f,line);
				if(std::string::npos!=line.find(vertLine)) {
					break;
				}
			}
			std::string token;
			std::stringstream ss(line);
			ss >> token >> token;
			unsigned nvert;
			ss >> nvert;
			bool hasTex = false;
			while (true) {
				std::getline(f,line);
				if(std::string::npos != line.find(faceLine)) {
					break;
				}
				if(std::string::npos != line.find(texLine)) {
					hasTex=true;
				}
			}
			std::stringstream ss1(line);
			ss1 >> token >> token;
			unsigned nface;
			ss1 >> nface;
			while (true) {
				std::getline(f, line);
				if(std::string::npos != line.find(endHeaderLine)) {
					break;
				}
			}
			
			v.resize(nvert);
			t.resize(nface);
			if (hasTex) {
				tex.resize(nvert);
			}
			for (unsigned ii = 0; ii < nvert; ++ii) {
				for (unsigned jj = 0; jj < 3; ++jj) {
					f >> v[ii][jj];
				}
				if (hasTex) {
					for (unsigned jj = 0; jj < 2; ++jj) {
						f >> tex[ii][jj];
					}
					tex[ii][1] = 1 - tex[ii][1];;
				}
			}
			for (unsigned ii = 0; ii < nface; ++ii) {
				unsigned nidx;
				f >> nidx;
				for (unsigned jj = 0; jj < 3; ++jj) {
					f >> t[ii][jj];
				}
			}
		}

		void readObj(std::istream &f)
		{
			std::string line;
			std::string vTok("v");
			std::string fTok("f");
			std::string texTok("vt");
			char bslash = '/' , space=' ';
			std::string tok;
			while (1) {
				std::getline(f,line);
				if (f.eof()) {
					break;
				}
				if (line == "#end"){
					break;
				}
				if (line.size() < 3) {
					continue;
				}
				if (line.at(0)=='#') {
					continue;
				}
				std::stringstream ss(line);
				ss >> tok;
				if(tok == vTok) {
					V3 vec;
					ss >> vec[0] >> vec[1] >> vec[2];
					v.push_back(vec);
				} else if (tok == fTok) {
					bool hasTexture = false;
					if (line.find(bslash) != std::string::npos) {
						std::replace(line.begin(), line.end(), bslash, space);
						hasTexture = true;
					}
					std::stringstream facess(line);
					facess >> tok;
					std::vector<int> vidx;
					std::vector<int> texIdx;
					unsigned x;
					while (facess >> x){
						vidx.push_back(x);
						if(hasTexture){
							facess >> x;
							texIdx.push_back(x);
						}
					}
					texIdx.resize(vidx.size());
					for (unsigned ii = 0;ii < vidx.size() - 2; ++ii){
						V3i trig, textureId;
						trig[0] = vidx[0]-1;
						textureId[0] = texIdx[0]-1;
						for (unsigned jj = 1; jj < 3; ++jj) {
							trig[jj] = vidx[ii+jj]-1;
							textureId[jj] = texIdx[ii+jj]-1;
						}
						t.push_back(trig);
						texId.push_back(textureId);
					}
				} else if (tok==texTok) {
					Vec2 texcoord;
					ss >> texcoord[0];
					ss >> texcoord[1];
					tex.push_back(texcoord);
				}
			}
			std::cout << "Num Triangles: " << t.size() << "\n";
		}

		
		void saveObj(const char * filename)
		{
			std::ofstream out(filename);
			if(!out.good()){
				std::cout<<"cannot open output file"<<filename<<"\n";
				return;
			}
			save(out);
			out.close();
		}
		
		inline void append(const Mesh &m)
		{
			unsigned offset = (unsigned) v.size();
			unsigned ot = (unsigned) t.size();
			
			v.insert(v.end(), m.v.begin(), m.v.end());
			t.insert(t.end(), m.t.begin(), m.t.end());
			
			for(unsigned ii = ot; ii < t.size(); ++ii){
				for(int jj = 0; jj < 3; ++jj){
					t[ii][jj] += offset;
				}
			}
		}
		
		inline Mesh & operator= (const Mesh& m)
		{
			v = m.v;
			t = m.t;
			n = m.n;
			return *this;
		}
		
	};
	
	inline bool ptInBox(const V3 &mn,
						const V3 mx,
						const V3 &x)
	{
		for (unsigned dim = 0; dim < 3; ++dim){
			if (x[dim] < mn[dim] || x[dim] > mx[dim]){
				return false;
			}
		}
		return true;
	}
	
	///@brief cube [0,1]^3
	const static V3 CUBE_VERT[8] = {
		V3 (0, 0, 0),
		V3 (1, 0, 0),
		V3 (1, 1, 0),
		V3 (0, 1, 0),
		V3 (0, 0, 1),
		V3 (1, 0, 1),
		V3 (1, 1, 1),
		V3 (0, 1, 1)
	};
	
	const static V3i CUBE_TRIG[12] = {
		V3i(0, 3, 1),
		V3i(1, 3, 2),
		V3i(5, 4, 0),
		V3i(5, 0, 1),
		V3i(6, 5, 1),
		V3i(1, 2, 6),
		V3i(3, 6, 2),
		V3i(3, 7, 6),
		V3i(4, 3, 0),
		V3i(4, 7, 3),
		V3i(7, 4, 5),
		V3i(7, 5, 6)};
	
	const static Mesh UNIT_CUBE(CUBE_VERT, CUBE_TRIG);
	
	inline void makeCube(Mesh &m,
						 const V3 &mn,
						 const V3 &mx)
	{
		V3 ss = mx - mn;
		m = UNIT_CUBE;
		for (unsigned ii = 0; ii < m.v.size(); ++ii) {
			m.v[ii][0] = mn[0] + ss[0] * m.v[ii][0];
			m.v[ii][1] = mn[1] + ss[1] * m.v[ii][1];
			m.v[ii][2] = mn[2] + ss[2] * m.v[ii][2];
		}
	}
	
	void BBox(const Mesh &m, V3 &mn, V3 &mx)
	{
		mn = m.v[0];
		mx = m.v[0];
		for (unsigned ii = 1; ii < m.v.size(); ++ii){
			for (unsigned dim = 0; dim < 3; ++dim)
				mn[dim] = _min(mn[dim], m.v[ii][dim]);
			for (unsigned dim = 0; dim < 3; ++dim)
				mx[dim] = _max(mx[dim], m.v[ii][dim]);
		}
	}
    
    V3 getBboxSize(const Mesh &m) {
        V3 min = m.v[0];
        V3 max = m.v[0];
        for (unsigned i = 1; i < m.v.size(); ++i)
        {
            if ( m.v[i].x < min.x ) min.x = m.v[i].x;
            if ( m.v[i].y < min.y ) min.y = m.v[i].y;
            if ( m.v[i].z < min.z ) min.z = m.v[i].z;
            if ( m.v[i].x > max.x ) max.x = m.v[i].x;
            if ( m.v[i].y > max.y ) max.y = m.v[i].y;
            if ( m.v[i].z > max.z ) max.z = m.v[i].z;
        }
        return V3(max.x-min.x, max.y-min.y, max.z-min.z);

    }
    V3 getMinSize(const Mesh &m) {
        V3 min = m.v[0];
        for (unsigned i = 1; i < m.v.size(); ++i)
        {
            if ( m.v[i].x < min.x ) min.x = m.v[i].x;
            if ( m.v[i].y < min.y ) min.y = m.v[i].y;
            if ( m.v[i].z < min.z ) min.z = m.v[i].z;
        }
        return min;
    }
}

#endif