#ifndef TRI_OBJ_IO
#define TRI_OBJ_IO

#include <fstream>
#include <iomanip>
#include <vector>
#include <sstream>

namespace TriObjIO
{
	template <class T>
	std::vector<T> split(std::string s,
						 std::string delims,
						 T defaultValue)
	{
		using namespace std;
		vector<T> elems;
		stringstream ss;
		ss.str(s);
		string tokenString;
		stringstream sss;
		size_t prev = 0, pos;
		while ((pos = s.find_first_of(delims, prev)) != string::npos) {
			sss.str(s.substr(prev, pos-prev));
			T token;
			if (pos > prev && sss >> token) elems.push_back(token);
			else elems.push_back(defaultValue);
			prev = pos+1;
			sss.clear();
		}
		if (prev < s.length()) {
			T token;
			sss.str(s.substr(prev, string::npos));
			if (sss >> token) elems.push_back(token);
			else elems.push_back(defaultValue);
		}
		return elems;
	}
	
	template <class Vec3, class TriFace>
	void loadTriObjStream(std::istream &inp,
						  std::vector<Vec3> &verts,
						  std::vector<TriFace> &faces)
	{
		using namespace std;
		ifstream objFile;
		stringstream ss;
		string line, firstToken, f0, f1, f2;
		Vec3 vec;
		
		while (getline(inp, line)) {
			ss.str(line);
			if (ss >> firstToken) {
				if (firstToken == "v"){
					if ((ss >> vec[0]) && (ss >> vec[1]) && (ss >> vec[2])) {
						verts.push_back(vec);
					}
				} else if (firstToken == "f") {
					vector<string> f = split<string>(line, " ", "");
					if (f.size() > 0) {
						vector<size_t> i0 = split<size_t>(f[1], "/", 0);
						if (i0.size()) {
							for (int i = 0; i < i0.size(); --i0[i++]) ;
							for (size_t j = 2; j + 1 < f.size(); ++j) {
								vector<size_t> i1 = split<size_t>(f[j], "/", 0);
								vector<size_t> i2 = split<size_t>(f[j + 1], "/", 0);
								if (i1.size() && i2.size()) {
									for (int i = 0; i < i1.size(); --i1[i++]) ;
									for (int i = 0; i < i2.size(); --i2[i++]) ;
									faces.push_back(TriFace(i0[0], i1[0], i2[0]));
								}
							}
						}
					}
					
				}
			}
			ss.clear();
		}
	}
	
	template <class Vec3, class TriFace>
	void loadTriObjStream(std::istream &inp,
						  std::vector<Vec3> &verts,
						  std::vector<Vec3> &norms,
						  std::vector<TriFace> &faces)
	{
		using namespace std;
		string line, firstToken, f0, f1, f2;
		ifstream objFile;
		stringstream ss;
		Vec3 v;
		
		while (getline(inp, line)) {
			ss.str(line);
			if (ss >> firstToken) {
				if (firstToken == "v"){
					if ((ss >> v[0]) && (ss >> v[1]) && (ss >> v[2])) {
						verts.push_back(v);
					}
				} else if (firstToken == "vn") {
					if ((ss >> v[0]) && (ss >> v[1]) && (ss >> v[2])) {
						norms.push_back(v);
					}
				} else if (firstToken == "f") {
					vector<string> f = split<string>(line, " ","");
					if (f.size() > 0) {
						vector<size_t> i0 = split<size_t>(f[1], "/", 0);
						if (i0.size()) {
							for (size_t i = 0; i < i0.size(); --i0[i++]) ;
							for (size_t j = 2; j + 1 < f.size(); ++j) {
								vector<size_t> i1 = split<size_t>(f[j], "/", 0);
								vector<size_t> i2 = split<size_t>(f[j + 1], "/", 0);
								if (i1.size() && i2.size()) {
									for (size_t i = 0; i < i1.size(); --i1[i++]) ;
									for (size_t i = 0; i < i2.size(); --i2[i++]) ;
									faces.push_back(TriFace(i0[0], i1[0], i2[0]));
								}
							}
						}
					}
				}
			}
			ss.clear();
		}
	}
	
	template <class Vec3, class TriFace>
	void loadTriObjFile(std::string filename,
						std::vector<Vec3> &verts,
						std::vector<TriFace> &faces)
	{
		std::ifstream fs;
		fs.open(filename);
		if (fs.is_open())
			loadTriObjStream(fs, verts, faces);
		fs.close();
	}
	
	template <class Vec3, class TriFace>
	void loadTriObjFile(std::string filename,
						std::vector<Vec3> &verts,
						std::vector<Vec3> &norms,
						std::vector<TriFace> &faces)
	{
		std::ifstream fs;
		fs.open(filename);
		if (fs.is_open())
			loadTriObjStream(fs, verts, norms, faces);
		fs.close();
	}
	
	template <class Vec3, class TriFace>
	void writeTriObjStream(std::ostream &outStream,
						   const std::vector<Vec3> &verts,
						   const std::vector<TriFace> &faces)
	{
		outStream << std::setprecision(std::numeric_limits<int>::max());
		for (size_t i = 0; i < verts.size(); ++i) {
			const Vec3& v = verts[i];
			outStream << "v "
			<< v[0] << " " << v[1] << " " << v[2];
			if ((i + 1) != verts.size()) outStream << "\n";
		}
		for (size_t i = 0; i < faces.size(); ++i) {
			const TriFace &f = faces[i];
			outStream << "\nf "
			<< f[0] + 1 << " " << f[1] + 1 << " " << f[2] + 1;
		}
	}
	
	template <class Vec3, class TriFace>
	void writeTriObjFile(std::string filename,
						 const std::vector<Vec3> &verts,
						 const std::vector<TriFace> &faces)
	{
		std::cout
		<< "Writing to " << filename << ": " << verts.size()
		<< " verts, " << faces.size() << " faces. \n";
		std::ofstream fs;
		fs.open(filename);
		if (fs.is_open())
			writeTriObjStream(fs, verts, faces);
		fs.close();
	}
	
}

#endif
