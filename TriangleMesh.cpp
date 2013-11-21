///////////////////////////////////////////////////////////////////////////////
//
// Representation and rendering facilities for a triangularized surface mesh
// 
// 2009 Karsten Noe
//
// Read the blog at cg.alexandra.dk for more information
//
///////////////////////////////////////////////////////////////////////////////

#include "TriangleMesh.h"
#include "objloader.h"
#include <limits>
#include <iostream>


RenderTriangleMesh::RenderTriangleMesh(TriangleMesh &m) : m(m)
{
	this->generateVBOs();
}

void RenderTriangleMesh::updateVBOs()
{
// cout<< "go"<<endl;
	for(unsigned int f=0;f<m.getFaces().size();f++)
	{
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				vertices[f*3*3+i*3+ j]= m.getFaces()[f][i].getPos().e[j];
				normals[f*3*3+i*3+ j]= m.getFaces()[f][i].getNormal().e[j];
			}
		}
	}

	glEnableClientState(GL_VERTEX_ARRAY);
	glBindBufferARB(GL_ARRAY_BUFFER_ARB, vboid[0]);
	glBufferDataARB(GL_ARRAY_BUFFER_ARB, ( m.getFaces().size() * 3 * 3 )*sizeof(GLfloat), vertices, GL_DYNAMIC_DRAW_ARB);

	glEnableClientState(GL_NORMAL_ARRAY);
	glBindBufferARB(GL_ARRAY_BUFFER_ARB, vboid[1]);
	glBufferDataARB(GL_ARRAY_BUFFER_ARB, ( m.getFaces().size() * 3 * 3 )*sizeof(GLfloat), normals, GL_DYNAMIC_DRAW_ARB);


}

void RenderTriangleMesh::generateVBOs()
{

	vertices = (GLfloat*) malloc(sizeof(GLfloat)*m.getFaces().size()*3*3);
	normals = (GLfloat*) malloc(sizeof(GLfloat)*m.getFaces().size()*3*3);

	for(unsigned int f=0;f<m.getFaces().size();f++)
	{
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<3;j++)
			{
				vertices[f*3*3+i*3+ j]= m.getFaces()[f][i].getPos().e[j];
				normals[f*3*3+i*3+ j]= m.getFaces()[f][i].getNormal().e[j];
			}
		}

	}

	glGenBuffersARB(2, vboid);

	glEnableClientState(GL_VERTEX_ARRAY);
	glBindBufferARB(GL_ARRAY_BUFFER_ARB, vboid[0]);
	glBufferDataARB(GL_ARRAY_BUFFER_ARB, ( m.getFaces().size() * 3 * 3 )*sizeof(GLfloat), vertices, GL_STATIC_DRAW_ARB);

	glEnableClientState(GL_NORMAL_ARRAY);
	glBindBufferARB(GL_ARRAY_BUFFER_ARB, vboid[1]);
	glBufferDataARB(GL_ARRAY_BUFFER_ARB, ( m.getFaces().size() * 3 * 3 )*sizeof(GLfloat), normals, GL_STATIC_DRAW_ARB);

}

#define BUFFER_OFFSET(i) ((char *)NULL + (i))


void RenderTriangleMesh::draw()
{
	this->updateVBOs();

	glBindBufferARB(GL_ARRAY_BUFFER_ARB, vboid[0]);
	glVertexPointer(3, GL_FLOAT, 0, BUFFER_OFFSET(0));
	glBindBufferARB(GL_ARRAY_BUFFER_ARB, vboid[1]);	
	glNormalPointer(GL_FLOAT,0,BUFFER_OFFSET(0));


	// Enable arrays
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);

	// Draw arrays
	glDrawArrays(GL_TRIANGLES, 0, (GLsizei) m.getFaces().size()*3);
	//glDrawElements(GL_TRIANGLES, m.getFaces().size(), GL_UNSIGNED_INT, NULL);

	// Disable arrays
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	//DisableClientState(COLOR_ARRAY);
}


TriangleFace::TriangleFace(void)
{
}

TriangleFace::~TriangleFace(void)
{
}

TriangleFace::TriangleFace(MeshVertex *p1, MeshVertex *p2, MeshVertex *p3)
{
	particles[0] = p1;
	particles[1] = p2;
	particles[2] = p3;

	normalUpdateNeeded = true;
}

void TriangleFace::updateNormal()
{
	normalUpdateNeeded = true;
}

Vector3 TriangleFace::getNormal()
{
	if (normalUpdateNeeded)
	{
		Vector3 &v1 = particles[0]->getPos();
		Vector3 &v2 = particles[1]->getPos();
		Vector3 &v3 = particles[2]->getPos();
		Vector3 tmp = (cross(v2-v1, v3-v1));
		tmp.makeUnitVector();
		normal = tmp;
		normalUpdateNeeded = false;
	}
	return normal;
}

MeshVertex::MeshVertex(void)
: position(Vector3(0,0,0))
{
}


MeshVertex::MeshVertex(Vector3 p)
{
	position = p; 
	normalUpdateNeeded = true;
}


MeshVertex::~MeshVertex(void)
{
}

void MeshVertex::addFace(TriangleFace *face)
{
	faces.push_back(face);
}

//!! force normal update
void MeshVertex::updateNormal()
{
	normalUpdateNeeded = true;
}

Vector3 MeshVertex::getNormal()
{
	if (normalUpdateNeeded )
	{
		int count = 0;
		Vector3 res(0,0,0);
		for (unsigned int i = 0; i<faces.size(); i++)
		{
			res = res + faces[i]->getNormal();

			count ++;	
		}
		res.makeUnitVector();
		normal = res;
		normalUpdateNeeded = false;
	}
	//	cout<<normal<< endl;
	return normal;

}


TriangleMesh::TriangleMesh(string filename)
{
	ObjLoader obj = ObjLoader(filename.c_str());

	particles.reserve(particles.size()+obj.v.size()*2);

	int startIndex = (int) particles.size();

	for(unsigned int i=0; i<obj.v.size(); i++)
	{
		particles.push_back(MeshVertex(obj.v[i]));
	}

	faces.reserve((int) obj.f.size()*2);

	for(unsigned int i=0; i<obj.f.size(); i++)
	{
		addFace(&particles[obj.f[i].v[0]-1+startIndex],&particles[obj.f[i].v[1]-1+startIndex],&particles[obj.f[i].v[2]-1+startIndex]);

		particles[obj.f[i].v[0]-1+startIndex].addFace( &getFaces()[getFaces().size()-1]);
		particles[obj.f[i].v[1]-1+startIndex].addFace( &getFaces()[getFaces().size()-1]);
		particles[obj.f[i].v[2]-1+startIndex].addFace( &getFaces()[getFaces().size()-1]);
	}

}

TriangleMesh::~TriangleMesh(void)
{
}

BoundingSphere TriangleMesh::getBoundingSphere()
{
	Vector3 center; 
	for( unsigned int i=0; i<this->getParticles().size(); i++)
	{
		center += this->getParticles()[i].getPos();;
	}
	center /= this->getParticles().size();

	float radius=0;

	for( unsigned int i=0; i<this->getParticles().size(); i++)
	{
		radius = std::max(radius,(center-this->getParticles()[i].getPos()).length());
	}

	return BoundingSphere(radius,center);

}

void TriangleMesh::updateToBoundingSphere(BoundingSphere boundingSphere)
{
	/* Same as getBoundingSphere here to get (maybe I should store boundingSphere on initialization) */
	Vector3 center;
	
	for( unsigned int i=0; i<this->getParticles().size(); i++)
	{
		center += this->getParticles()[i].getPos();;
	}
	center /= this->getParticles().size();

	float radius=0;

	for( unsigned int i=0; i<this->getParticles().size(); i++)
	{
		radius = std::max(radius,(center-this->getParticles()[i].getPos()).length());
	}

	/* Now update (This might take a long time if the mesh's too big... multiple multiplications) */
	float scale = boundingSphere.radius/radius;

	for( unsigned int i=0; i<this->getParticles().size(); i++)
	{
		Vector3 oldpos = this->getParticles()[i].getPos();

		Vector3 newpos;
		newpos[0] = (oldpos[0] - center[0])*scale+boundingSphere.center[0];
		newpos[1] = (oldpos[1] - center[1])*scale+boundingSphere.center[1];
		newpos[2] = (oldpos[2] - center[2])*scale+boundingSphere.center[2];

		this->getParticles()[i].setPos(newpos);;
	}
	//Vector3 centerDisplacement = boundingSphere.center - center;

	//for( unsigned int i=0; i<this->getParticles().size(); i++)
	//{
	//	Vector3 oldpos = this->getParticles()[i].getPos();

	//	Vector3 newpos;
	//	newpos[0] = oldpos[0] + centerDisplacement[0];
	//	newpos[1] = oldpos[1] + centerDisplacement[1];
	//	newpos[2] = oldpos[2] + centerDisplacement[2];

	//	this->getParticles()[i].setPos(newpos);;
	//}
}


void TriangleMesh::updateNormals(void)
{
	vector<MeshVertex> &particles = getParticles();
	vector<TriangleFace> &faces = getFaces();

	for( unsigned int i=0; i<faces.size(); i++)
	{
		TriangleFace & f = faces[i];
		f.updateNormal();
	}

	for( unsigned int i=0; i<particles.size(); i++)
	{
		MeshVertex & p = particles[i];
		p.updateNormal();
	}

}



