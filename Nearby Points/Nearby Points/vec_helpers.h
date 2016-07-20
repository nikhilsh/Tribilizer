#ifndef vec_helpers
#define vec_helpers

#include <vector>

#ifdef _WIN32
#include "GL/freeglut.h"
#elif __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

namespace DF
{
	struct TriFace
	{
		unsigned p[3];
		
		inline TriFace(unsigned p0, unsigned p1, unsigned p2)
		{ p[0] = p0; p[1] = p1; p[2] = p2; }
		
		inline TriFace() {}
		
		inline unsigned& operator [] (unsigned i) { return p[i]; };
		inline unsigned operator [] (unsigned i) const { return p[i]; };
	};
	
	template <class S, unsigned Dims>
struct Bounds
	{
		Vec<S, Dims> lower, upper;
	
		Bounds() {}
		
		Bounds(const Vec<S, Dims> &lower, const Vec<S, Dims> &upper):
		lower(lower), upper(upper) {};
	};
	
	template <class S, unsigned Dims>
	inline Bounds<S, Dims> getBounds(const std::vector<Vec<S, Dims> > &verts)
	{
		if (verts.size()) {
			Bounds<S, Dims> bounds(verts[0], verts[0]);
			for (unsigned i = 0; i < verts.size(); ++i) {
				bounds.lower = min(bounds.lower, verts[i]);
				bounds.upper = max(bounds.upper, verts[i]);
			}
			return bounds;
		}
		return Bounds<S, Dims>();
	}
	
	template <class S,
	class D0, unsigned R0, unsigned C0,
	class D1, unsigned R1, unsigned C1,
	class D2, unsigned R2, unsigned C2>
	inline Vec<S, 3> planeNorm(const tVec3Base<S, D0, R0, C0> &p0,
							   const tVec3Base<S, D1, R1, C1> &p1,
							   const tVec3Base<S, D2, R2, C2> &p2)
	{
		return (p1 - p0).cross(p2 - p0).normalized();
	}
	
	template <class S, class TriFace>
	void calculateNorms(const std::vector<Vec<S, 3> > &inVerts,
						const std::vector<TriFace> &inFaces,
						std::vector<Vec<S, 3> > &outNorms)
	{
		outNorms.clear();
		outNorms.resize(inVerts.size());
		for (unsigned i = 0; i < inFaces.size(); ++i) {
			TriFace f = inFaces[i];
			Vec<S, 3> n = planeNorm(inVerts[f[0]],
									inVerts[f[1]],
									inVerts[f[2]]);
			outNorms[f[0]] += n;
			outNorms[f[1]] += n;
			outNorms[f[2]] += n;
		}
		for (unsigned i = 0; i < outNorms.size(); ++i)
			outNorms[i] = outNorms[i].normalized();
	}
	
	template <class S, class TriFace>
	void drawTriFaces(const std::vector<Vec<S, 3> > &verts,
					  const std::vector<Vec<S, 3> > &norms,
					  const std::vector<TriFace> &faces)
	{
		glBegin(GL_TRIANGLES);
		for (unsigned i = 0; i < faces.size(); ++i) {
			const TriFace &f = faces[i];
			glNormal3f(norms[f[0]][0], norms[f[0]][1], norms[f[0]][2]);
			glVertex3f(verts[f[0]][0], verts[f[0]][1], verts[f[0]][2]);
			glNormal3f(norms[f[1]][0], norms[f[1]][1], norms[f[1]][2]);
			glVertex3f(verts[f[1]][0], verts[f[1]][1], verts[f[1]][2]);
			glNormal3f(norms[f[2]][0], norms[f[2]][1], norms[f[2]][2]);
			glVertex3f(verts[f[2]][0], verts[f[2]][1], verts[f[2]][2]);
		}
		glEnd();
	}
	
	template <class S, class TriFace>
	void drawTriEdges(const std::vector<Vec<S, 3> > &verts,
					  const std::vector<TriFace> &faces)
	{
		glBegin(GL_LINES);
		for (unsigned i = 0; i < faces.size(); ++i) {
			const TriFace &f = faces[i];
			glVertex3f(verts[f[0]][0], verts[f[0]][1], verts[f[0]][2]);
			glVertex3f(verts[f[1]][0], verts[f[1]][1], verts[f[1]][2]);
			glVertex3f(verts[f[1]][0], verts[f[1]][1], verts[f[1]][2]);
			glVertex3f(verts[f[2]][0], verts[f[2]][1], verts[f[2]][2]);
			glVertex3f(verts[f[2]][0], verts[f[2]][1], verts[f[2]][2]);
			glVertex3f(verts[f[0]][0], verts[f[0]][1], verts[f[0]][2]);
		}
		glEnd();
	}
	
	template <class S>
	Mat<S, 4> getCenterIntoUnitVolumeTransform(const std::vector<Vec<S, 3> > &verts)
	{
		struct { inline S operator()(S a, S b) { return a > b ? a : b; }; } max;
		
		Bounds<S, 3> bounds = getBounds(verts);
		Vec<S, 3> diff = bounds.upper - bounds.lower;
		Vec<S, 3> trans = (bounds.upper + bounds.lower) / -2.0;
		S s = 1 / max(max(diff[0], diff[1]), diff[2]);
		return scale(translate(Mat<S, 4>::Identity(), trans), Vec<S, 3>(s));
	}
	
	template <class S,
	class D0, unsigned R0, unsigned C0,
	class D1, unsigned R1, unsigned C1,
	class D2, unsigned R2, unsigned C2>
	inline S triArea(const tVec3Base<S, D0, R0, C0> &p0,
					 const tVec3Base<S, D1, R1, C1> &p1,
					 const tVec3Base<S, D2, R2, C2> &p2)
	{
		return 0.5 * (p1 - p0).cross(p2 - p0).norm();
	}
	
	template <class S>
	S totalTriArea(const std::vector<Vec<S, 3> > &verts,
				   const std::vector<TriFace> &faces)
	{
		S a = 0;
		for (unsigned i = 0; i < faces.size(); ++i) {
			const TriFace &f = faces[i];
			a += triArea(verts[f[0]], verts[f[1]], verts[f[2]]);
		}
		return a;
	}
	
}

#endif
