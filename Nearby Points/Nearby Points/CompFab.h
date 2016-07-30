//
//  CompFab.h
//  voxelizer
//
//  Created by David Levin on 2/3/14.
//
//

#ifndef voxelizer_CompFab_h
#define voxelizer_CompFab_h

#include "vec.h"

namespace CompFab
{
	typedef float Real;

	typedef DF::Vec<Real, 2> Vec2;
	typedef DF::Vec<Real, 3> V3;
	typedef DF::Vec<Real, 4> Vec4;

	typedef DF::Vec<unsigned, 2> Vec2i;
	typedef DF::Vec<unsigned, 3> V3i;
	typedef DF::Vec<unsigned, 4> Vec4i;
	
	struct Ray
    {
		V3 origin, direction;
		
		inline Ray() {}
		
		inline Ray(const V3 &origin, const V3 &direction): origin(origin), direction(direction) {}
    };
    
    class Plane {
    public:
        Plane() : mDistance(0) {}
        float distance() const { return mDistance; }
        float distanceToPoint(const V3 &vertex) const {
            return (vertex.dot(mNormal)) - mDistance;
        }
        void setNormal(V3 normal) { mNormal = normal; }
        void setDistance(float distance) { mDistance = distance; }
    protected:
        V3    mNormal;    // normalized Normal-Vector of the plane
        float   mDistance;  // shortest distance from plane to Origin
    };
    
    
    struct LineSegment {
        LineSegment(V3 p0=V3(), V3 p1=V3()) {
            v[0]=p0; v[1]=p1;
        }
        V3 v[2];
    };
    
	struct Triangle
    {
		V3 v0, v1, v2;
		
		inline Triangle() {}
		
		inline Triangle(const V3 &v0, const V3 &v1, const V3 &v2): v0(v0), v1(v1), v2(v2) {}
		
		inline V3 & operator [] (unsigned i)
		{
			if (i == 0) return v0;
			if (i == 1) return v1;
			if (i == 2) return v2;
			return v0;
		}
		
		inline V3 operator [] (unsigned i) const
		{
			if (i == 0) return v0;
			if (i == 1) return v1;
			if (i == 2) return v2;
			return v0;
		}
        int intersectPlane(const Plane &plane, LineSegment &ls) const{
            size_t countFront=0, countBack=0;
            float distance1 = plane.distanceToPoint(v0);
            float distance2 = plane.distanceToPoint(v1);
            float distance3 = plane.distanceToPoint(v2);
            if (distance1<0) {
                countBack++;
            }
            else {
                countFront++;
            }
            if (distance2<0) {
                countBack++;
            }
            else {
                countFront++;
            }
            if (distance3<0) {
                countBack++;
            }
            else {
                countFront++;
            }
            if (countBack == 3) {
                return -1;
            }
            else if (countFront == 3) {
                return 1;
            }
            size_t lines[] = {0,1,1,2,2,0};
            std::vector<V3> intersectPoints;
            
            Triangle tri;
            tri.v0 = v0;
            tri.v1 = v1;
            tri.v2 = v2;
            for (size_t i=0; i<3; ++i) {
                
                const V3 &a = tri[lines[i*2+0]];
                const V3 &b = tri[lines[i*2+1]];
                const float da = plane.distanceToPoint(a);
                const float db = plane.distanceToPoint(b);
                if (da*db<0) {
                    const float s = da/(da-db); // intersection factor (between 0 and 1)
                    V3 bMinusa = b-a;
                    intersectPoints.push_back(a+bMinusa*s);
                }
                else if (0==da) { // plane falls exactly on one of the three Triangle vertices
                    if (intersectPoints.size()<2) {
                        intersectPoints.push_back(a);
                    }
                }
                else if (0==db) { // plane falls exactly on one of the three Triangle vertices
                    if (intersectPoints.size()<2) {
                        intersectPoints.push_back(b);
                    }
                }
                
                if (2==intersectPoints.size()) {
                    // Output the intersecting line segment object
                    ls.v[0]=intersectPoints[0];
                    ls.v[1]=intersectPoints[1];
                    return 0;
                }
                return -2;
            }
            return 1;
        }

	};
	
	struct VoxelGrid
	{
		V3 lowerLeft;
		unsigned dimX, dimY, dimZ, size;
		float spacing;
		bool *insideArray;

		inline VoxelGrid(): insideArray(NULL) {}
		
		inline void copy(const VoxelGrid &g)
		{
			lowerLeft = g.lowerLeft;
			dimX = g.dimX;
			dimY = g.dimY;
			dimZ = g.dimZ;
			size = g.dimX * dimY * dimZ;
			spacing = g.spacing;
			insideArray = new bool[size];
		}
		
		inline VoxelGrid(const VoxelGrid &g)
		{
			copy(g);
		}
		
		inline VoxelGrid & operator = (const VoxelGrid &g)
		{
			if (insideArray) delete [] insideArray;
			copy(g);
			return *this;
		}
		
		inline VoxelGrid(const V3 &lowerLeft, unsigned dimX, unsigned dimY, unsigned dimZ, Real spacing):
		lowerLeft(lowerLeft),
		dimX(dimX),
		dimY(dimY),
		dimZ(dimZ),
		size(dimX * dimY * dimZ),
		spacing(spacing),
		insideArray(new bool[size])
		{
			for (unsigned ii = 0; ii < size; ++ii) insideArray[ii] = false;
		}
		
		inline ~VoxelGrid()
		{
			delete [] insideArray;
		}
		
		inline bool & isInside(unsigned i, unsigned j, unsigned k)
		{
			return insideArray[k * (dimX * dimY) + j * dimY + i];
		}
	};
}



#endif
