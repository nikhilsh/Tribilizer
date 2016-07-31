//
// For this Demo:
// Use the arrow keys to move the sphere up down left right.
// Use W and S keys to move the sphere forward and backward.
//
//

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

#ifdef _WIN32
#include "GL/freeglut.h"
#elif __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "vec.h"
#include "camera.h"
#include "bvh.h"
#include "vec_helpers.h"
#include "tri_obj_io.h"
#include "cork.h"
#include "CompFab.h"
#include "Mesh.h"

typedef float Real;

using namespace DF;

typedef Vec<Real, 2> V2;
typedef Vec<Real, 3> V3;
typedef Vec<Real, 4> V4;

typedef Mat<Real, 2> M2;
typedef Mat<Real, 3> M3;
typedef Mat<Real, 4> M4;

#define EPSILON 1e-9


// This is the camera
Camera camera;

// BVH for points.
ActiveRodMesh::BVH<Real, 3, true> bvh;

// BVH for triangles. Not used here.
ActiveRodMesh::BVH<Real, 3, false> triBVH;


V3 sphereCenter = 0;
Real sphereRadius = 0.5;

// The vertices, normals and faces
std::vector<V3> verts;
std::vector<V3> norms;
std::vector<TriFace> faces;

struct Center
{
	V3 position;
	bool hit;
	inline Center() {}
	inline Center(const V3 &position): position(position), hit(0) {}
};

std::vector<Center> centers;

struct MouseFunction { enum { NONE, ARC_BALL }; };
int mouseFunction = MouseFunction::NONE;

int mouseX = 0;
int mouseY = 0;

int cutPoint1 = 0;
int cutPoint2 = 0;

typedef std::vector<CompFab::Triangle> TriangleList;

TriangleList triangleList;
CompFab::VoxelGrid voxelGrid;

//-------------------------------------------------------------------

// Declarations of functions whose implementations occur later.
void keyboardFunc( unsigned char key, int x, int y);
void specialFunc( int key, int x, int y );
void mouseFunc(int button, int state, int x, int y);
void motionFunc(int x, int y);
void reshapeFunc(int w, int h);
void drawScene(void);
void initRendering();


// This function is called whenever a "Normal" key press is
// received.
void keyboardFunc(unsigned char key, int x, int y)
{
	if (key == 27) exit(0);
	
	key = tolower(key);
	if (key == 'r') {
	} else if (key == 't') {
	} else if (key == 'y') {
	} else if (key == 'u') {
	} else if (key == 'a') {
	} else if (key == 'z') {
	} else if (key == 'w') {
		sphereCenter.z -= 0.1;
	} else if (key == 's') {
		sphereCenter.z += 0.1;
	} else {
		std::cout << "Unhandled key press " << key << "." << std::endl;
	}
	
}

enum { KEY_LEFT = 100, KEY_UP = 101, KEY_RIGHT = 102, KEY_DOWN = 103 };

// This function is called whenever a "Special" key press is
// received.  Right now, it's handling the arrow keys.
void specialFunc(int key, int x, int y)
{
	if (key == KEY_LEFT) {
		sphereCenter.x -= 0.1;
	} else if (key == KEY_RIGHT) {
		sphereCenter.x += 0.1;
	} else if (key == KEY_DOWN) {
		sphereCenter.y -= 0.1;
	} else if (key == KEY_UP) {
		sphereCenter.y += 0.1;
	}

}

void specialUpFunc(int key, int x, int y)
{
	//cout << "Unhandled key press " << key << "." << endl;
}

//  Called when mouse button is pressed.
void mouseFunc(int button, int state, int x, int y)
{
	mouseX = x;
	mouseY = y;
	
	int modifiers = glutGetModifiers();
	bool shiftOn = modifiers & 1;
	//bool altOn = modifiers & 4;
	//bool ctrlOn = modifiers & 2;
	
	if (modifiers == 0) {
		if (state == GLUT_DOWN && mouseFunction == MouseFunction::NONE) {
			mouseFunction = MouseFunction::ARC_BALL;
			switch (button) {
				case GLUT_LEFT_BUTTON:
					camera.mouseClick(Camera::LEFT, x, y);
					break;
				case GLUT_RIGHT_BUTTON:
					camera.mouseClick(Camera::RIGHT, x,y);
				default:
					break;
			}
		} else {
			camera.mouseRelease(x,y);
			mouseFunction = MouseFunction::NONE;
		}
		glutPostRedisplay();
		
	}
}

// Called when mouse is moved while button pressed.
void motionFunc(int x, int y)
{
	mouseX = x;
	mouseY = y;
	
	if (mouseFunction == MouseFunction::ARC_BALL) {
		camera.mouseDrag(x,y);
	}
	glutPostRedisplay();
	
}

// Called when the window is resized
// w, h - width and height of the window in pixels.
void reshapeFunc(int w, int h)
{
	camera.setDimensions(w,h);
	
	camera.setViewport(0,0,w,h);
	camera.applyViewport();
	
	// Set up a perspective view, with square aspect ratio
	glMatrixMode(GL_PROJECTION);
	
	camera.setPerspective(50);
	glLoadMatrixf(camera.getProjectionMatrix());
	
}

// Initialize OpenGL's rendering modes
void initRendering()
{
	glEnable(GL_DEPTH_TEST);   // Depth testing must be turned on
	glEnable(GL_LIGHTING);     // Enable lighting calculations
	glEnable(GL_LIGHT0);       // Turn on light #0.
	
	glEnable(GL_NORMALIZE);
	
	// Setup polygon drawing
	glShadeModel(GL_SMOOTH);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	
	// Enable alpha
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	// Clear to black
	glClearColor(0,0,0,1);
}

// This function is responsible for displaying the object.
void drawScene(void)
{
	// Clear the rendering window
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	
	// Light color (RGBA)
	GLfloat Lt0diff[] = {1.0,1.0,1.0,1.0};
	GLfloat Lt0pos[] = {3.0,3.0,5.0,1.0};
	glLightfv(GL_LIGHT0, GL_DIFFUSE, Lt0diff);
	glLightfv(GL_LIGHT0, GL_POSITION, Lt0pos);
	
	glLoadMatrixf(camera.getModelViewMatrix());

	
	static bool sphereListInited = false;
	static GLuint smallSphereList = 0;
	if (!sphereListInited) {
		sphereListInited = true;
		smallSphereList = glGenLists(1);
		glNewList(smallSphereList, GL_COMPILE);
		glutSolidSphere(0.08, 8, 8);
		glEndList();
	}
	
	{
		GLfloat normalColor[] = {0.2f, 0.2f, 0.2f, 1.f};
		GLfloat hitColor[] = {0.9f, 0.2f, 0.2f, 1.f};
		GLfloat specColor[] = {1.f, 1.f, 1.f, 1.f};
		GLfloat shininess[] = {50.0f};
		
		glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, normalColor );
		glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, specColor );
		glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, shininess );

		std::vector<unsigned> pointsWithinSphere;
		bvh.findPointsWithinRadius(pointsWithinSphere, sphereCenter, sphereRadius);
		
		for (unsigned i = 0; i < centers.size(); ++i) {
			centers[i].hit = false;
		}
		for (unsigned i = 0; i < pointsWithinSphere.size(); ++i) {
			centers[pointsWithinSphere[i]].hit = true;
		}
		
		for (unsigned i = 0; i < centers.size(); ++i) {
			const V3 &p = centers[i].position;
			if (centers[i].hit) {
				glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, hitColor );
			} else {
				glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, normalColor );
			}
			glTranslatef(p.x, p.y, p.z);
			glCallList(smallSphereList);
			glTranslatef(-p.x, -p.y, -p.z);
		}

	}
	
	{
		GLfloat diffColor[] = {0.5f, 0.5f, 0.5f, 0.3f};
		GLfloat specColor[] = {0.8f, 0.8f, 0.8f, 1.f};
		GLfloat shininess[] = {50.0f};
		
		glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, diffColor );
		glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, specColor );
		glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, shininess );
		
		glTranslatef(sphereCenter.x, sphereCenter.y, sphereCenter.z);
		glutSolidSphere(sphereRadius, 10, 10);
		glTranslatef(-sphereCenter.x, -sphereCenter.y, -sphereCenter.z);

	}
	
	drawTriFaces(verts, norms, faces);
	
	// This draws the coordinate axes when you're rotating, to
	// keep yourself oriented.
	if( mouseFunction == MouseFunction::ARC_BALL )
	{
		glPushMatrix();
		V3 eye(M4(camera.getTransform().inverse()) * V4(camera.getCenter(), 1));
		glTranslatef( eye[0], eye[1], eye[2] );
		
		// Save current state of OpenGL
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		
		// This is to draw the axes when the mouse button is down
		glDisable(GL_LIGHTING);
		glLineWidth(2);
		glPushMatrix();
		glScaled(5.0,5.0,5.0);
		glBegin(GL_LINES);
		glColor4f(1,0.5,0.5,1); glVertex3f(0,0,0); glVertex3f(1,0,0);
		glColor4f(0.5,1,0.5,1); glVertex3f(0,0,0); glVertex3f(0,1,0);
		glColor4f(0.5,0.5,1,1); glVertex3f(0,0,0); glVertex3f(0,0,1);
		
		glColor4f(0.5,0.5,0.5,1);
		glVertex3f(0,0,0); glVertex3f(-1,0,0);
		glVertex3f(0,0,0); glVertex3f(0,-1,0);
		glVertex3f(0,0,0); glVertex3f(0,0,-1);
		
		glEnd();
		glPopMatrix();
		
		glPopAttrib();
		glPopMatrix();
	}
	
	//bvh.draw();
	
	// Dump the image to the screen.
	glutSwapBuffers();
}

void timerFunc(int t)
{
	glutPostRedisplay();
	glutTimerFunc(t, &timerFunc, t);
}

inline bool rayTriangleIntersection(const CompFab::Ray &ray, const CompFab::Triangle &triangle)
{
    /********* ASSIGNMENT *********/
    /* Ray-Triangle intersection test: Return 1 if ray intersects triangle,
     * 0 otherwise */
    V3 e1, e2;
    V3 P, Q, T;
    float t;
    
    e1 = triangle.v1 - triangle.v0;
    e2 = triangle.v2 - triangle.v0;
    
    P = ray.direction.cross(e2);
    float determinant = e1.dot(P);
    
    if (determinant > -EPSILON && determinant < EPSILON) {
        return 0;
    }
    
    float invDet = 1.0f / determinant;
    
    T = ray.origin - triangle.v0;
    float u = T.dot(P) * invDet;
    if (u < 0.0f || u > 1.0f) {
        return 0;
    }
    
    Q = T.cross(e1);
    
    float v = ray.direction.dot(Q) * invDet;
    if (v < 0.0f || u + v > 1.0f) {
        return 0;
    }
    
    t = e2.dot(Q) * invDet;
    
    if (t > EPSILON) {
        return 1;
    }
    return 0;
}

inline unsigned numSurfaceIntersections(const V3 &voxelPos, const V3 &dir)
{
    unsigned numHits = 0;
    
    /********* ASSIGNMENT *********/
    /* Check and return the number of times a ray cast in direction dir,
     * from voxel center voxelPos intersects the surface */
    CompFab::Ray r;
    for (int i = 0; i < triangleList.size(); i++) {
        r = CompFab::Ray(voxelPos, dir);
        if (rayTriangleIntersection(r, triangleList[i])) {
            numHits += 1;
        }
    }
    return numHits;
}

void loadModel(char *filename)
{
	TriObjIO::loadTriObjFile("bunny.obj", verts, faces);
	
	M4 transform = (scale(M4::Identity(), V3(1)) *
					getCenterIntoUnitVolumeTransform(verts));
	
	for (unsigned i = 0; i < verts.size(); ++i)
		verts[i] = V3(transform * V4(verts[i], 1));
	
	calculateNorms(verts, faces, norms);

	for (unsigned i = 0; i < faces.size(); ++i) {
		TriFace &f = faces[i];
		V3 p = (verts[f[0]] + verts[f[1]] + verts[f[2]]) / 3;
		Center c(p);
		centers.push_back(c);
		bvh.add(i, p); // Add the point
	}
    bvh.build(); // Build
	
	for (unsigned i = 0; i < faces.size(); ++i) {
		TriFace &f = faces[i];
		triBVH.add(i,
				   min(verts[f[0]], min(verts[f[1]], verts[f[2]])),
				   max(verts[f[0]], max(verts[f[1]], verts[f[2]])));
	}
	bvh.build(); // Build
	triBVH.build();
    
}

void inputMeshToCorkMesh ();
void calculateSplittingOfMesh() {
    
    //need to change dividing by also
    const unsigned dim = 32; //dimension of voxel grid (e.g. 32x32x32)
    int mass[dim] = {};
    int voxcount[dim] = {};
    float totalcount = 0;
    float massheight = 0;
    float totalmass = 0;
    int com[dim] = {};
    int average = 0;
    int cut1 = 0;
    int cut2 = 0;
    int pointer = 0;
    
    CompFab::Mesh tempMesh("bunny.obj", true);
    
    triangleList.resize(tempMesh.t.size());
    
    //copy triangles to global list
    for (unsigned i = 0; i < tempMesh.t.size(); ++i) {
        triangleList[i].v0 = tempMesh.v[tempMesh.t[i][0]];
        triangleList[i].v1 = tempMesh.v[tempMesh.t[i][1]];
        triangleList[i].v2 = tempMesh.v[tempMesh.t[i][2]];
    }
    
    //Create Voxel Grid
    V3 bbMax, bbMin;
    BBox(tempMesh, bbMin, bbMax);
    
    //Build Voxel Grid
    Real bbX = bbMax[0] - bbMin[0];
    Real bbY = bbMax[1] - bbMin[1];
    Real bbZ = bbMax[2] - bbMin[2];
    Real spacing;
    
    if (bbX > bbY && bbX > bbZ) {
        spacing = bbX / (Real)(dim - 2);
    } else if (bbY > bbX && bbY > bbZ) {
        spacing = bbY / (Real)(dim - 2);
    } else {
        spacing = bbZ / (Real)(dim - 2);
    }
    
    V3 hspacing(0.5 * spacing);

    voxelGrid = CompFab::VoxelGrid(bbMin - hspacing, dim, dim, dim, spacing);
    
    
    
    V3 voxelPos;
    V3 direction(1, 0, 0);
    /********* ASSIGNMENT *********/
    /* Iterate over all voxels in voxelGrid and test whether they are inside our outside of the
     * surface defined by the triangles in triangleList */
    for (int i = 0; i < dim; i++) {
        for (int k = 0; k < dim; k++) {
            for (int j = 0; j < dim; j++) {
                //find i-voxel
                V3 voxeli = V3(i * voxelGrid.spacing, 0, 0);
                //find j-voxel
                V3 voxelj = V3(0, j * voxelGrid.spacing, 0);
                //find k-voxel
                V3 voxelk = V3(0, 0, k * voxelGrid.spacing);
                //find the new position from iteration, start from origin then add each dimension
                voxelPos = voxelGrid.lowerLeft + voxeli + voxelj + voxelk;
                voxelGrid.isInside(i, j, k) = numSurfaceIntersections(voxelPos, direction) % 2;
                if (voxelGrid.isInside(i, j, k)) { mass[j]++; voxcount[j]++;}
            }
        }
    }
    //compilation of all COM after each layer is change to 0.2 mass
    
    for (int p = dim-1; p > 0-1; p--) {
        mass[p] = 0.2*mass[p];
        for (int q = 0; q < dim; q++) {
            totalmass = totalmass + mass[q];
            massheight = massheight + (mass[q] * q);
        }
        com[p]= ceil(massheight / totalmass);
        //printf("COM%d is %d\n", p, com[p]);
        totalmass = 0;
        massheight = 0;
    }
    int cut = dim - 1;
    for (int p = dim - 1; p > 0; p--) {
        if (com[cut] >= com[p]) {
            cut = p;
        }
    }
    
    for (int i = 0; i < dim; i++) {
        totalcount = totalcount + voxcount[i];
        printf("Layer%d voxcount=%d\n", i, voxcount[i]);
    }
    
    // sorts in ascending order
    std::sort(std::begin(voxcount), std::end(voxcount));
    
    //percentile calculation
    pointer = 0.25 * dim;
    if (floor(pointer) == pointer)
        average = (voxcount[pointer] + voxcount[pointer + 1]) * 0.5;
    else
        average = voxcount[(int) ceil(pointer)];
    
    //printf("cut@ %d, threshold=%d\n", cut,average);
    
    
    for (int q = cut; q < dim; q++) {
        cut1 = q;
        if (voxcount[q] > average) {
            break;
        }
    }
    for (int r = cut; r > 0; r--) {
        cut2 = r;
        if (voxcount[r] > average) {
            break;
        }
    }
    printf("lowestCOM cut@Layer%d with count%d. Threshold=%d therefore cut1@Layer%d with count%d cut2@Layer%d with count%d\n", cut,voxcount[cut],average,cut1,voxcount[cut1],cut2,voxcount[cut2]);
    
    cutPoint1 = cut1;
    cutPoint2 = cut2;
    
    inputMeshToCorkMesh();
}

void inputMeshToCorkMesh () {
    
    CorkTriMesh overallMesh;
    
    overallMesh.n_vertices = verts.size();
    overallMesh.n_triangles = faces.size();
    
    overallMesh.triangles = new uint[(overallMesh.n_triangles) * 3];
    overallMesh.vertices  = new float[(overallMesh.n_vertices) * 3];
    
    for(uint i=0; i < overallMesh.n_triangles; i++) {
        TriFace &f = faces[i];
        overallMesh.triangles[3*i+0] = f.p[0];
        overallMesh.triangles[3*i+1] = f.p[1];
        overallMesh.triangles[3*i+2] = f.p[2];
    }
    
    for(uint i=0; i<overallMesh.n_vertices; i++) {
        overallMesh.vertices[3*i+0] = verts[i][0];
        overallMesh.vertices[3*i+1] = verts[i][1];
        overallMesh.vertices[3*i+2] = verts[i][2];
    }
    
    Bounds<Real, 3> bounds = getBounds(verts);

    CorkTriMesh cubeMesh;
    cubeMesh.n_vertices = 8;
    cubeMesh.n_triangles = 12;
    cubeMesh.vertices = new float[3*8];
    cubeMesh.triangles = new uint[3*12];
    
    std::vector<TriFace> cubeFaces(12);
    std::vector<V3> cubeVerts(8);

    if (cutPoint1 != cutPoint2) {
#define DEF_CUBE_VERT(I, X,Y,Z) \
cubeMesh.vertices[I*3+0] = X; \
cubeMesh.vertices[I*3+1] = Y; \
cubeMesh.vertices[I*3+2] = Z;
        
        const Real &xMin = bounds.lower.x;
        const Real &yMin = bounds.lower.y;
        const Real &zMin = bounds.lower.z;
        
        const Real &xMax = bounds.upper.x;
        const Real &yMax = bounds.upper.y;
        const Real &zMax = bounds.upper.z;
        
        //lowest point
        //    const Real yCut = yMin + (yMax - yMin) * 1.0 / 3.0;
        
        const Real yCut = yMin + (yMax - yMin) * (cutPoint1/32.0);
        printf("cut point 1: %f", yCut);
        
        DEF_CUBE_VERT(0, xMin,yCut,zMin);
        DEF_CUBE_VERT(1, xMax,yCut,zMin);
        DEF_CUBE_VERT(2, xMax,yCut,zMax);
        DEF_CUBE_VERT(3, xMin,yCut,zMax);
        DEF_CUBE_VERT(4, xMin,yMin-0.0001,zMin);
        DEF_CUBE_VERT(5, xMax,yMin-0.0001,zMin);
        DEF_CUBE_VERT(6, xMax,yMin-0.0001,zMax);
        DEF_CUBE_VERT(7, xMin,yMin-0.0001,zMax);
        
#undef DEF_CUBE_VERT
        
        // clockwise
#define DEF_CUBE_FACE(I, A,B,C,D) \
cubeMesh.triangles[I*6+0] = A; \
cubeMesh.triangles[I*6+1] = D; \
cubeMesh.triangles[I*6+2] = B; \
cubeMesh.triangles[I*6+3] = B; \
cubeMesh.triangles[I*6+4] = D; \
cubeMesh.triangles[I*6+5] = C;
        
        DEF_CUBE_FACE(0, 0,1,2,3);
        DEF_CUBE_FACE(1, 1,5,6,2);
        DEF_CUBE_FACE(2, 7,6,5,4);
        DEF_CUBE_FACE(3, 3,7,4,0);
        DEF_CUBE_FACE(4, 3,2,6,7);
        DEF_CUBE_FACE(5, 0,4,5,1);
        
#undef DEF_CUBE_FACE
        
        CorkTriMesh cubeMesh2;
        cubeMesh2.n_vertices = 8;
        cubeMesh2.n_triangles = 12;
        cubeMesh2.vertices = new float[3*8];
        cubeMesh2.triangles = new uint[3*12];
        
        std::vector<TriFace> cubeFaces2(12);
        std::vector<V3> cubeVerts2(8);
        
#define DEF_CUBE_VERT(I, X,Y,Z) \
cubeMesh2.vertices[I*3+0] = X; \
cubeMesh2.vertices[I*3+1] = Y; \
cubeMesh2.vertices[I*3+2] = Z;
        
        //higher point
        //    const Real yCut2 = yMax - (yMax - yMin) * 1.0 / 3.0;
        const Real yCut2 = yMax - (yMax - yMin) * cutPoint2/32.0;
        printf("cut point 2: %f", yCut2);
        
        
        DEF_CUBE_VERT(0, xMin,yMax+0.0001,zMin);
        DEF_CUBE_VERT(1, xMax,yMax+0.0001,zMin);
        DEF_CUBE_VERT(2, xMax,yMax+0.0001,zMax);
        DEF_CUBE_VERT(3, xMin,yMax+0.0001,zMax);
        DEF_CUBE_VERT(4, xMin,yCut2,zMin);
        DEF_CUBE_VERT(5, xMax,yCut2,zMin);
        DEF_CUBE_VERT(6, xMax,yCut2,zMax);
        DEF_CUBE_VERT(7, xMin,yCut2,zMax);
        
#undef DEF_CUBE_VERT
        
        // clockwise
#define DEF_CUBE_FACE(I, A,B,C,D) \
cubeMesh2.triangles[I*6+0] = A; \
cubeMesh2.triangles[I*6+1] = D; \
cubeMesh2.triangles[I*6+2] = B; \
cubeMesh2.triangles[I*6+3] = B; \
cubeMesh2.triangles[I*6+4] = D; \
cubeMesh2.triangles[I*6+5] = C;
        
        DEF_CUBE_FACE(0, 0,1,2,3);
        DEF_CUBE_FACE(1, 1,5,6,2);
        DEF_CUBE_FACE(2, 7,6,5,4);
        DEF_CUBE_FACE(3, 3,7,4,0);
        DEF_CUBE_FACE(4, 3,2,6,7);
        DEF_CUBE_FACE(5, 0,4,5,1);
        
#undef DEF_CUBE_FACE
        
        //    TriObjIO::writeTriObjFile("outC.obj", cubeVerts, cubeFaces);
        //    TriObjIO::writeTriObjFile("outD.obj", cubeVerts2, cubeFaces2);
        
        CorkTriMesh outputMesh;
        CorkTriMesh outputMesh2;
        
        CorkTriMesh outputMesh3;
        CorkTriMesh outputMesh4;
        
        computeDifference(overallMesh, cubeMesh, &outputMesh);
        computeDifference(outputMesh, cubeMesh2, &outputMesh2);
        
        computeIntersection(cubeMesh, overallMesh, &outputMesh3);
        computeIntersection(cubeMesh2, overallMesh, &outputMesh4);
        
        delete [] overallMesh.triangles;
        delete [] overallMesh.vertices;
        
        delete [] cubeMesh.triangles;
        delete [] cubeMesh.vertices;
        
        delete [] cubeMesh2.triangles;
        delete [] cubeMesh2.vertices;
        
        std::vector<TriFace> outFaces(outputMesh2.n_triangles);
        std::vector<V3> outVerts(outputMesh2.n_vertices);
        
        for(uint i=0; i < outputMesh2.n_triangles; i++) {
            outFaces[i][0] = outputMesh2.triangles[3*i+0];
            outFaces[i][1] = outputMesh2.triangles[3*i+1];
            outFaces[i][2] = outputMesh2.triangles[3*i+2];
        }
        
        for(uint i=0; i<outputMesh2.n_vertices; i++) {
            outVerts[i][0] = outputMesh2.vertices[3*i+0];
            outVerts[i][1] = outputMesh2.vertices[3*i+1];
            outVerts[i][2] = outputMesh2.vertices[3*i+2];
        }
        
        freeCorkTriMesh(&outputMesh2);
        
        TriObjIO::writeTriObjFile("out1.obj", outVerts, outFaces);
        
        std::vector<TriFace> outFaces2(outputMesh3.n_triangles);
        std::vector<V3> outVerts2(outputMesh3.n_vertices);
        
        for(uint i=0; i < outputMesh3.n_triangles; i++) {
            outFaces2[i][0] = outputMesh3.triangles[3*i+0];
            outFaces2[i][1] = outputMesh3.triangles[3*i+1];
            outFaces2[i][2] = outputMesh3.triangles[3*i+2];
        }
        
        for(uint i=0; i<outputMesh3.n_vertices; i++) {
            outVerts2[i][0] = outputMesh3.vertices[3*i+0];
            outVerts2[i][1] = outputMesh3.vertices[3*i+1];
            outVerts2[i][2] = outputMesh3.vertices[3*i+2];
        }
        
        freeCorkTriMesh(&outputMesh3);
        
        TriObjIO::writeTriObjFile("out2.obj", outVerts2, outFaces2);
        
        std::vector<TriFace> outFaces3(outputMesh4.n_triangles);
        std::vector<V3> outVerts3(outputMesh4.n_vertices);
        
        for(uint i=0; i < outputMesh4.n_triangles; i++) {
            outFaces3[i][0] = outputMesh4.triangles[3*i+0];
            outFaces3[i][1] = outputMesh4.triangles[3*i+1];
            outFaces3[i][2] = outputMesh4.triangles[3*i+2];
        }
        
        for(uint i=0; i<outputMesh4.n_vertices; i++) {
            outVerts3[i][0] = outputMesh4.vertices[3*i+0];
            outVerts3[i][1] = outputMesh4.vertices[3*i+1];
            outVerts3[i][2] = outputMesh4.vertices[3*i+2];
        }
        
        freeCorkTriMesh(&outputMesh4);
        
        TriObjIO::writeTriObjFile("out3.obj", outVerts3, outFaces3);
        
    }
    else {
#define DEF_CUBE_VERT(I, X,Y,Z) \
cubeMesh.vertices[I*3+0] = X; \
cubeMesh.vertices[I*3+1] = Y; \
cubeMesh.vertices[I*3+2] = Z;
        
        const Real &xMin = bounds.lower.x;
        const Real &yMin = bounds.lower.y;
        const Real &zMin = bounds.lower.z;
        
        const Real &xMax = bounds.upper.x;
        const Real &yMax = bounds.upper.y;
        const Real &zMax = bounds.upper.z;
        
        //lowest point
        //    const Real yCut = yMin + (yMax - yMin) * 1.0 / 3.0;
        
        const Real yCut = yMin + (yMax - yMin) * (cutPoint1/32.0);
        printf("cut point 1: %f", yCut);
        
        DEF_CUBE_VERT(0, xMin,yCut,zMin);
        DEF_CUBE_VERT(1, xMax,yCut,zMin);
        DEF_CUBE_VERT(2, xMax,yCut,zMax);
        DEF_CUBE_VERT(3, xMin,yCut,zMax);
        DEF_CUBE_VERT(4, xMin,yMin-0.0001,zMin);
        DEF_CUBE_VERT(5, xMax,yMin-0.0001,zMin);
        DEF_CUBE_VERT(6, xMax,yMin-0.0001,zMax);
        DEF_CUBE_VERT(7, xMin,yMin-0.0001,zMax);
        
#undef DEF_CUBE_VERT
        
        // clockwise
#define DEF_CUBE_FACE(I, A,B,C,D) \
cubeMesh.triangles[I*6+0] = A; \
cubeMesh.triangles[I*6+1] = D; \
cubeMesh.triangles[I*6+2] = B; \
cubeMesh.triangles[I*6+3] = B; \
cubeMesh.triangles[I*6+4] = D; \
cubeMesh.triangles[I*6+5] = C;
        
        DEF_CUBE_FACE(0, 0,1,2,3);
        DEF_CUBE_FACE(1, 1,5,6,2);
        DEF_CUBE_FACE(2, 7,6,5,4);
        DEF_CUBE_FACE(3, 3,7,4,0);
        DEF_CUBE_FACE(4, 3,2,6,7);
        DEF_CUBE_FACE(5, 0,4,5,1);
        
#undef DEF_CUBE_FACE
        
        //    TriObjIO::writeTriObjFile("outC.obj", cubeVerts, cubeFaces);
        //    TriObjIO::writeTriObjFile("outD.obj", cubeVerts2, cubeFaces2);
        
        CorkTriMesh outputMesh;
        CorkTriMesh outputMesh2;
        
        
        computeIntersection(overallMesh, cubeMesh, &outputMesh);
        computeDifference(overallMesh, cubeMesh, &outputMesh2);
        
//        delete [] overallMesh.triangles;
//        delete [] overallMesh.vertices;
        
        delete [] cubeMesh.triangles;
        delete [] cubeMesh.vertices;
        
        
        std::vector<TriFace> outFaces(outputMesh.n_triangles);
        std::vector<V3> outVerts(outputMesh.n_vertices);
        
        for(uint i=0; i < outputMesh.n_triangles; i++) {
            outFaces[i][0] = outputMesh.triangles[3*i+0];
            outFaces[i][1] = outputMesh.triangles[3*i+1];
            outFaces[i][2] = outputMesh.triangles[3*i+2];
        }
        
        for(uint i=0; i<outputMesh.n_vertices; i++) {
            outVerts[i][0] = outputMesh.vertices[3*i+0];
            outVerts[i][1] = outputMesh.vertices[3*i+1];
            outVerts[i][2] = outputMesh.vertices[3*i+2];
        }
        
        
        TriObjIO::writeTriObjFile("out1.obj", outVerts, outFaces);
        
        std::vector<TriFace> outFaces2(outputMesh2.n_triangles);
        std::vector<V3> outVerts2(outputMesh2.n_vertices);
        
        for(uint i=0; i < outputMesh2.n_triangles; i++) {
            outFaces2[i][0] = outputMesh2.triangles[3*i+0];
            outFaces2[i][1] = outputMesh2.triangles[3*i+1];
            outFaces2[i][2] = outputMesh2.triangles[3*i+2];
        }
        
        for(uint i=0; i<outputMesh2.n_vertices; i++) {
            outVerts2[i][0] = outputMesh2.vertices[3*i+0];
            outVerts2[i][1] = outputMesh2.vertices[3*i+1];
            outVerts2[i][2] = outputMesh2.vertices[3*i+2];
        }
        
        freeCorkTriMesh(&outputMesh2);
        
        TriObjIO::writeTriObjFile("out2.obj", outVerts2, outFaces2);
        
    }
    
    //    std::vector<V3> vertsLink;
    //    std::vector<V3> normsLink;
    //    std::vector<TriFace> facesLink;
    //    
    //    TriObjIO::loadTriObjFile("linkset.obj", vertsLink, facesLink);
//    
//    M4 transform = (scale(M4::Identity(), V3(1)) *
//                    getCenterIntoUnitVolumeTransform(vertsLink));
//    
//    for (unsigned i = 0; i < vertsLink.size(); ++i)
//        vertsLink[i] = V3(transform * V4(vertsLink[i], 1));
//    
//    calculateNorms(vertsLink, facesLink, normsLink);
//
//    CorkTriMesh linkMesh;
//    
//    linkMesh.n_vertices = vertsLink.size();
//    linkMesh.n_triangles = vertsLink.size();
//    
//    linkMesh.triangles = new uint[(linkMesh.n_triangles) * 3];
//    linkMesh.vertices  = new float[(linkMesh.n_vertices) * 3];
//    
//    for(uint i=0; i < linkMesh.n_triangles; i++) {
//        TriFace &f = facesLink[i];
//        linkMesh.triangles[3*i+0] = f.p[0];
//        linkMesh.triangles[3*i+1] = f.p[1];
//        linkMesh.triangles[3*i+2] = f.p[2];
//    }
//    
//    for(uint i=0; i<linkMesh.n_vertices; i++) {
//        linkMesh.vertices[3*i+0] = verts[i][0];
//        linkMesh.vertices[3*i+1] = verts[i][1];
//        linkMesh.vertices[3*i+2] = verts[i][2];
//    }
//    
////    M4 transform2 = (translate(M4::Identity(), V3(1, ((yMax - yMin) * 1.0 / 3.0), 1)));
////    
////    for (unsigned i = 0; i < vertsLink.size(); ++i)
////        vertsLink[i] = V3(transform2 * V4(vertsLink[i], 1));
//    
//    CorkTriMesh combineMesh;
//
////    computeUnion(linkMesh, linkMesh, &combineMesh);
//
//    std::vector<TriFace> combineFaces(linkMesh.n_triangles);
//    std::vector<V3> combineVerts(linkMesh.n_vertices);
//    
//    for(uint i=0; i < linkMesh.n_triangles; i++) {
//        combineFaces[i][0] = linkMesh.triangles[3*i+0];
//        combineFaces[i][1] = linkMesh.triangles[3*i+1];
//        combineFaces[i][2] = linkMesh.triangles[3*i+2];
//    }
//    
//    for(uint i=0; i<linkMesh.n_vertices; i++) {
//        combineVerts[i][0] = linkMesh.vertices[3*i+0];
//        combineVerts[i][1] = linkMesh.vertices[3*i+1];
//        combineVerts[i][2] = linkMesh.vertices[3*i+2];
//    }
//    
//    freeCorkTriMesh(&linkMesh);
//    
//    TriObjIO::writeTriObjFile("outcombine.obj", combineVerts, combineFaces);
}

// Main routine.
// Set up OpenGL, define the callbacks and start the main loop
int main( int argc, char* argv[] )
{
	glutInit( &argc, argv );
	
	// We're going to animate it, so double buffer
	glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );
	
	// Initial parameters for window position and size
	glutInitWindowPosition( 60, 60 );
	glutInitWindowSize( 600, 600 );
	
	camera.setDimensions( 600, 600 );
	
	camera.setDistance( 10 );
	camera.setCenter(V3(0));
	
//	glutCreateWindow("BVH Test");
	
	// Load the model
	loadModel(argv[1]);
    calculateSplittingOfMesh();

	// Initialize OpenGL parameters.
//	initRendering();
//	
//	// Set up callback functions for key presses.
//	glutKeyboardFunc(keyboardFunc);   // Handles "normal" ascii symbols
//	glutSpecialFunc(specialFunc);     // Handles "special" keyboard keys on pressed down
//	glutSpecialUpFunc(specialUpFunc); // Handles "special" keyboard keys on release
//	
//	// Set up callback functions for mouse
//	glutMouseFunc(mouseFunc);
//	glutMotionFunc(motionFunc);
//	
//	// Set up the callback function for resizing windows
//	glutReshapeFunc( reshapeFunc );
//	
//	// Call this whenever window needs redrawing
//	glutDisplayFunc( drawScene );
//	
//	// Trigger timerFunc every 20 msec
//	glutTimerFunc(20, timerFunc, 20);
	
    
	// Start the main loop.  glutMainLoop never returns.
//	glutMainLoop();
	
	return 0;	// This line is never reached.
}
