#ifndef CAMERA_H
#define CAMERA_H

#ifdef _WIN32
#include "GL/freeglut.h"
#include <GL/glu.h>
#elif __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#include "vec.h"
#include <vector>

namespace DF
{
	class Camera
	{
	private:
		static inline float PI()
		{ return 3.14159265358979323846264338327950288f; }
		
	public:
		
		typedef Mat<float, 4> M4;
		typedef Vec<float, 3> V3;
		typedef Vec<float, 2> V2;
		
		typedef enum { NONE, LEFT, MIDDLE, RIGHT } Button;
		
		Camera(): startRot(1), currentRot(1) {}
		
		// You must call all of the set*() functions before you use this!
		// I didn't put it into the constructor because it's inconvenient
		// to initialize stuff in my opengl application.
		
		Camera & setDimensions(int w, int h)
		{
			dimensions[0] = w;
			dimensions[1] = h;
			return *this;
		}
		
		Camera & setViewport(int x, int y, int w, int h)
		{
			viewport[0] = x;
			viewport[1] = y;
			viewport[2] = w;
			viewport[3] = h;
			persp[1] = float( w ) / h;
			return *this;
		}
		
		Camera & setPerspective(float fovy)
		{
			persp[0] = fovy;
			return *this;
		}
		
		// Call from whatever UI toolkit
		Camera & mouseClick(Button button, int x, int y)
		{
			startClick[0] = x;
			startClick[1] = y;
			
			buttonState = button;
			
			switch (button) {
				case LEFT: currentRot = startRot; break;
				case MIDDLE: currentCenter = startCenter; break;
				case RIGHT: currentDistance = startDistance; break;
				default: break;
			}
			return *this;
		}
		
		Camera & mouseDrag(int x, int y)
		{
			switch (buttonState) {
				case LEFT: arcBallRotation(x,y); break;
				case MIDDLE: planeTranslation(x,y); break;
				case RIGHT: distanceZoom(x,y); break;
				default: break;
			}
			return *this;
		}
		
		Camera & mouseRelease(int x, int y)
		{
			startRot = currentRot;
			startCenter = currentCenter;
			startDistance = currentDistance;
			buttonState = NONE;
			return *this;
		}
		
		// Apply viewport, perspective, and modeling
		// use these instead of
		Camera & applyViewport()
		{
			glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
			return *this;
		}
		
		const Camera & applyViewport() const
		{
			glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
			return *this;
		}
		
		M4 getProjectionMatrix() const
		{
			float fovYRadians = persp[0] * PI() / 180.0f;
			float aspect = persp[1];
			float near = 0.1f;
			float far = 1000.f;
			return perspective(fovYRadians, aspect, near, far);
		}
		
		M4 getModelViewMatrix() const
		{
			V3 eye(0, 0, currentDistance);
			V3 center(0);
			V3 up(0,1,0);
			return (lookAt(eye, center, up) * currentRot *
					translate(M4::Identity(), -currentCenter) *
					getTransform());
		}
		
		float getWindowWidth() { return dimensions[0]; }
		
		float getWindowHeight() { return dimensions[1]; }
		
		// set for relevant vars
		Camera & setCenter(const V3 &center)
		{
			startCenter = currentCenter = center;
			return *this;
		}
		
		Camera & setRotation(const M4& rotation)
		{
			startRot = currentRot = rotation;
			return *this;
		}
		
		Camera & setDistance(float distance)
		{
			startDistance = currentDistance = distance;
			return *this;
		}
		
		// get for relevant vars
		V3 getCenter() const
		{
			return currentCenter;
		}
		
		M4 getRotation() const
		{
			return currentRot;
		}
		
		float getDistance() const
		{
			return currentDistance;
		}
		
		M4 getTransform() const
		{
			M4 t = M4::Identity();
			for (unsigned i = 0; i < transforms.size(); ++i)
				t = transforms[i] * t;
			return t;
		}
		
		Camera & pushTransform(const M4 &transform)
		{
			transforms.push_back(transform);
			return *this;
		}
		
		Camera & popTransform()
		{
			transforms.pop_back();
			return *this;
		}
		
	private:
		
		// States
		int dimensions[2];
		int startClick[2];
		Button buttonState;
		
		// For rotation
		M4 startRot, currentRot;
		
		// For translation
		float persp[2];
		int   viewport[4];
		V3 startCenter, currentCenter;
		
		// For zoom
		float startDistance;
		float currentDistance;
		
		std::vector<M4> transforms;
		
		void arcBallRotation(int x, int y)
		{
			float sx, sy, sz, ex, ey, ez;
			float scale;
			float sl, el;
			float dotprod;
			
			// find vectors from center of window
			sx = startClick[0] - ( dimensions[0] / 2.f );
			sy = startClick[1] - ( dimensions[1] / 2.f );
			ex = x - ( dimensions[0] / 2.f );
			ey = y - ( dimensions[1] / 2.f );
			
			// invert y coordinates (raster versus device coordinates)
			sy = -sy;
			ey = -ey;
			
			// scale by inverse of size of window and magical sqrt2 factor
			if (dimensions[0] > dimensions[1]) {
				scale = (float) dimensions[1];
			} else {
				scale = (float) dimensions[0];
			}
			
			scale = 1.f / scale;
			
			sx *= scale;
			sy *= scale;
			ex *= scale;
			ey *= scale;
			
			// project points to unit circle
			sl = hypot(sx, sy);
			el = hypot(ex, ey);
			
			if (sl > 1.f) {
				sx /= sl;
				sy /= sl;
				sl = 1.0;
			}
			if (el > 1.f) {
				ex /= el;
				ey /= el;
				el = 1.f;
			}
			
			// project up to unit sphere - find Z coordinate
			sz = sqrt(1.0f - sl * sl);
			ez = sqrt(1.0f - el * el);
			
			// rotate (sx,sy,sz) into (ex,ey,ez)
			
			// compute angle from dot-product of unit vectors (and double it).
			// compute axis from cross product.
			dotprod = sx * ex + sy * ey + sz * ez;
			
			if( dotprod != 1 )
			{
				V3 axis(sy * ez - ey * sz,
						sz * ex - ez * sx,
						sx * ey - ex * sy);
				
				axis = axis.normalized();
				
				float angle = 2.0f * acos( dotprod );
				
				currentRot = rotate(M4::Identity(), angle, axis);
				currentRot = currentRot * startRot;
			}
			else
			{
				currentRot = startRot;
			}
		}
		
		void planeTranslation(int x, int y)
		{
			// map window x,y into viewport x,y
			
			// start
			int sx = startClick[0] - viewport[0];
			int sy = startClick[1] - viewport[1];
			
			// current
			int cx = x - viewport[0];
			int cy = y - viewport[1];
			
			
			// compute "distance" of image plane (wrt projection matrix)
			float d = float(viewport[3])/2.0f / tan(persp[0] * PI() / 180.0f / 2.0f);
			
			// compute up plane intersect of clickpoint (wrt fovy)
			float su = -sy + viewport[3]/2.0f;
			float cu = -cy + viewport[3]/2.0f;
			
			// compute right plane intersect of clickpoint (ASSUMED FOVY is 1)
			float sr = (sx - viewport[2]/2.0f);
			float cr = (cx - viewport[2]/2.0f);
			
			V2 move(cr-sr, cu-su);
			
			// this maps move
			move *= -currentDistance / d;
			
			currentCenter = startCenter +
			+ move[0] * currentRot.block<1, 3>(0, 0).transposed()
			+ move[1] * currentRot.block<1, 3>(1, 0).transposed();
		}
		
		void distanceZoom(int x, int y)
		{
			int sy = startClick[1] - viewport[1];
			int cy = y - viewport[1];
			
			float delta = float(cy - sy) / viewport[3];
			
			// exponential zoom factor
			currentDistance = startDistance * exp(delta);
		}
	};

}


#endif
