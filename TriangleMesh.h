///////////////////////////////////////////////////////////////////////////////
//
// Representation and rendering facilities for a triangularized surface mesh
// 
// 2009 Karsten Noe
//
// Read the blog at cg.alexandra.dk for more information
//
///////////////////////////////////////////////////////////////////////////////


#pragma once

#include "Vector3.h"
#include <vector>

#include "glew/glew.h"

class BoundingSphere
{
public:
	float radius;
	Vector3 center;
	BoundingSphere(float radius, Vector3 center) : radius(radius),center(center) {};
	BoundingSphere() {};
};

class TriangleFace;

class MeshVertex
{

private:
	Vector3 position; 
	std::vector<TriangleFace *> faces;
	Vector3 normal;
	bool normalUpdateNeeded;

public:

	unsigned int index;

	MeshVertex(void);
	MeshVertex(Vector3 position);
	~MeshVertex(void);

	Vector3 &getPos() {return position;};
	void addFace(TriangleFace *face);
	const std::vector<TriangleFace*>& getFaces() {return faces;  };
	Vector3 getNormal();

	void setPos(Vector3 & pos) {position = pos; }
	bool isSurfaceParticle() {return faces.size()>0;}
	void updateNormal();
};

class TriangleFace
{
private:
	Vector3 normal;
	bool normalUpdateNeeded;

public:
	MeshVertex* particles[3];
	TriangleFace(void);
	TriangleFace(MeshVertex *p1, MeshVertex *p2, MeshVertex *p3);
	MeshVertex& operator[](int i) {  return *(particles[i]); }
	Vector3 getNormal();
	~TriangleFace(void);
	void updateNormal();
};

class TriangleMesh
{
protected:
	std::vector<MeshVertex> particles;
	std::vector<TriangleFace> faces;

public:
	TriangleMesh(string filename);
	~TriangleMesh(void);

	void reserveNumberOfFaces(unsigned int n) { faces.reserve(n); };
	void addFace(TriangleFace &f) {faces.push_back(f);};
	void addFace(MeshVertex *p1, MeshVertex *p2, MeshVertex *p3) {faces.push_back(TriangleFace(p1,p2,p3));};

	std::vector<TriangleFace>& getFaces() {return faces;};
	std::vector<MeshVertex>& getParticles() { return particles; };

	void updateNormals();
	void updateToBoundingSphere(BoundingSphere boundingSphere);

	BoundingSphere getBoundingSphere();
};


class RenderTriangleMesh
{
private: 
	TriangleMesh &m;

	GLuint vboid[2];

	GLfloat *vertices;
	GLfloat *normals;

public:
	RenderTriangleMesh(TriangleMesh &m);
	void draw();

private:
	void generateVBOs();
	void updateVBOs();

};