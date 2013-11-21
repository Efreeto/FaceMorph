#include <vector>
#include "Vector3.h"
#include "Vector2.h"

#ifndef _ObjLoader
#define _ObjLoader

using namespace std;


//f v/vt/vn v/vt/vn v/vt/vn 
// assume triangular faces
class Face
{
public:
	int v[3];
	int vt[3];
	int vn[3];

	unsigned char groupId;

	int tangent[3];//index to tangent and binormal - is per default same as vt, but may be changed by e.g. calculateTangentSpace

	Face();
	Face(int v1, int vt1, int vn1, int v2, int vt2, int vn2, int v3, int vt3, int vn3);
};


class ObjLoader
{
	public:
		
		vector<Vector3> v; // vertex
		vector<Vector2> vt; // vertex texture coordinates
		vector<Vector3> vn; // vertex normals
		vector<Face> f;  // faces

		vector<Vector2> vt1; // decal texture coordinates (MELODY SPECIFIC)
		vector<Vector3> tangent; // tangent (MELODY SPECIFIC)
		vector<Vector3> binormal; // binormals (MELODY SPECIFIC)

		ObjLoader() {};
		ObjLoader( const char *filename, unsigned short numGroups = 0, char** groupNames = 0x0 );


		void loadObj(const char *filename, unsigned short numGroups = 0, char** groupNames = 0x0);
		void loadStl(const char *filename);
		void saveObj(const char *filename);
		void saveStl(const char *filename);
		void transform(float m4x4[][4]);
		void translateVertices(Vector3 t);
		void scaleVertices(float s);
		void save();
		void grow(float factor);
		void SingleLaplacianSmooth(float factor=0.5);
		void ThomasSmooth(float factor=0.5, int zLimit=-1);
		void TaubinSmooth(float grow=0.6307, float iterations = 50, float smooth = -0.6732 );
//		void weld();
//		void weld(float EPS);
		void boundingBox(Vector3 &minV, Vector3&maxV);

		void interpolateNormals();

		void computeStaticTangentspace();
		void computeTangentSpace();

		void convertBumpMapToNormalMap(const char* bumpFilename, const char* normalMapFilename, float bumpFactor);
		void ObjLoader::convertBumpMapToObjectSpaceNormalMap(const char* bumpFilename, const char* normalMapFilename, float bumpFactor);
		void printInfo();
private:
	bool sameSide(Vector2 p1, Vector2 p2, Vector2 a, Vector2 b);
	bool pointIntTriangle(Vector2 p, Vector2 a, Vector2 b, Vector2 c);

};

#endif 