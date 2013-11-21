#include "ObjLoader.h"
#include <string>
#include <iostream>
#include <fstream>

#include "assert.h" 
using namespace std;

Face::Face()
{
	v[0]=0;v[1]=0;v[2]=0;
	vt[0]=0;vt[1]=0;vt[2]=0;
	vn[0]=0;vn[1]=0;vn[2]=0;
}

Face::Face(int v1, int vt1, int vn1, int v2, int vt2, int vn2, int v3, int vt3, int vn3)
{
	v[0]=v1;v[1]=v2;v[2]=v3;

	vt[0]=vt1;vt[1]=vt2;vt[2]=vt3;

	vn[0]=vn1;vn[1]=vn2;vn[2]=vn3;

	tangent[0]=vt1;
	tangent[1]=vt2;
	tangent[2]=vt3;
}

void ObjLoader::printInfo()
{
	cout << "v        " << v.size() << endl;
	cout << "vn       " << vn.size() << endl;
	cout << "vt       " << vt.size() << endl;
	cout << "f        " << f.size() << endl;
	cout << "tangent  " << tangent.size() << endl;
	cout << "binormal " << binormal.size() << endl;
}

void ObjLoader::grow(float factor)
{

	vector<Face>::iterator fi;  // faces
	for(fi=f.begin();fi!=f.end();fi++)
	{
		Face f = *fi;
		for(int i=0;i<3;i++)
		{
			v[ f.v[i]-1 ] = v[ f.v[i]-1 ]+ vn[ f.vn[i]-1 ]*factor;


		}

	}



}

void ObjLoader::saveStl(const char *filename)
{

	std::ofstream in(filename, ios::out | ios::binary);

	char header[80] = "Stl File, Generated through ObjLoader.cpp by Jesper Mosegaard";
	in.write(header,80);

	unsigned long numberOfFacets = f.size();
	unsigned int attribute = 0;
	in.write((char*)&numberOfFacets,4);
	Vector3 vec3;

	for(unsigned long i=0; i<numberOfFacets; i++)
	{
		
        vn.push_back(vec3);
		in.write((char*)&vn[f[i].vn[0]-1],4*3);

		for(int j=0; j<3; j++)
		{
			in.write((char*)&v[f[i].v[j]-1],4*3);
		}
		in.write((char*)&attribute,2);
	}
}

void ObjLoader::save()
{
	saveObj("ObjLoader_save.obj");
}

char* notZero(int i)
{
  char buffer[255];
  if (i!=0)
  {
	_itoa_s(i,buffer,10);
	return buffer;
  }
  else
  {
	return "";
  }

}

void ObjLoader::saveObj(const char *filename)
{
	std::ofstream out(filename);
	vector<Vector3>::iterator vi; // vertex
	vector<Vector2>::iterator vti; // vertex texture coordinates
	vector<Vector3>::iterator vni; // vertex normals
	vector<Face>::iterator fi;  // faces

	for(vi=v.begin();vi!=v.end();vi++)
	{
		Vector3 p = *vi;
		out << "v " << p.e[0] << " " << p.e[1] << " " <<p.e[2] << endl;
	
	}

	for(vti=vt.begin();vti!=vt.end();vti++)
	{
		Vector2 p = *vti;
		out << "vt " << p.e[0] << " " <<  p.e[1] << endl;
	}
	out << "#" << vt.size() << "texture coordinates" << endl;

	for(vni=vn.begin();vni!=vn.end();vni++)
	{
		Vector3 p = *vni;
		out << "vn " << p.e[0] << " " << p.e[1] << " " <<p.e[2] << endl;
	}

	for(fi=f.begin();fi!=f.end();fi++)
	{
		Face f = *fi;
		//&f.v[0],&f.vt[0],&(f.vn[0])

		out << "f " <<
			f.v[0] << "/" << notZero(f.vt[0]) << "/" << notZero(f.vn[0]) << " " <<
			f.v[1] << "/" << notZero(f.vt[1]) << "/" << notZero(f.vn[1]) << " " <<
			f.v[2] << "/" << notZero(f.vt[2]) << "/" << notZero(f.vn[2]) << endl;	
	}



	out.close();
}



ObjLoader::ObjLoader( const char *filename, unsigned short numGroups, char** groupNames )
{
	std::string f(filename);
	std::string kind = f.substr(f.length()-3,3);
	if ( kind=="obj") this->loadObj(filename, numGroups, groupNames );
	if ( kind=="stl") this->loadStl(filename);	

}

void ObjLoader::interpolateNormals()
{
	std::vector<Vector3>newNormal;
	std::vector<int> count;
	newNormal.resize(v.size());
	count.resize(v.size(),0);


	for(unsigned int i=0; i< f.size(); i++)
	{
		for(unsigned int j=0;j<3;j++)
		{	
			int nIndex = (f[i].v[j]-1);

			newNormal[nIndex]+=vn[ f[i].vn[j]-1 ];
			count[nIndex]++;
			f[i].vn[j]=f[i].v[j];
		}

	}

	vn.resize(count.size());
	for(int i=0;i<newNormal.size();i++)
	{
		vn[i] = newNormal[i]/float(count[i]);
	}



}

void ObjLoader::loadStl(const char *filename)
{

/* Binary STL, http://www.ennex.com/~fabbers/StL.asp */

	std::ifstream in(filename, ios::in | ios::binary);

	if(!in.good())
	{
		cerr << "ERROR: ObjLoader::loadStl(" << filename << ") file is not good" << endl;
		exit(-1);
	}
	char header[80];
	in.read(header,80);

	unsigned long numberOfFacets;
	unsigned int attribute;
	in.read((char*)&numberOfFacets,4);
	Vector3 vec3;

	for(unsigned long i=0; i<numberOfFacets; i++)
	{
		in.read((char*)&vec3,4*3);
        vn.push_back(vec3);

		for(int j=0; j<3; j++)
		{
			in.read((char*)&vec3,4*3);
			v.push_back(vec3);
		}
		f.push_back(Face(v.size()-1,0,vn.size(),v.size(),0,vn.size(),v.size()-2,0,vn.size()));
		//f.push_back(Face(v.size(),0,0,v.size()-1,0,0,v.size()-2,0,0));
		in.read((char*)&attribute,2);
	}
	//this->weld();
	//this->interpolateNormals();


/* ASCII STL */

		/* solid
 ...
 facet normal 0.00 0.00 1.00
    outer loop
      vertex  2.00  2.00  0.00
      vertex -1.00  1.00  0.00
      vertex  0.00 -1.00  0.00
    endloop
  endfacet
 ...
 endsolid*/

/*
	std::ifstream in(filename);

	if(!in.good())
	{
		cerr << "ERROR: ObjLoader::loadStl(" << filename << ") file is not good" << endl;
		exit(-1);
	}
	
	char buffer[255]; float f1,f2,f3;

	in.getline(buffer,255); // "solid"

	in.getline(buffer,255);
	while(!strcmp(buffer,"endsolid"))
		{
			if ( sscanf(buffer,"facet normal %f %f %f",f1,f2,f3)==3)
			{
				vn.push_back(Vector3(f1,f2,f3));
			}
			else
			{
				cerr << "error in format" << endl;
			}

			in.getline(buffer,255);

			for(int i=0; i<3; i++)
			{
				if ( sscanf(buffer,"vertex %f %f %f",f1,f2,f3)==3)
				{
					v.push_back(Vector3(f1,f2,f3));
				}
				else
				{
					cerr << "error in format" << endl;
				}
			}

			f.push_back(Face(v.size(),0,vn.size(),v.size()-1,0,vn.size(),v.size()-2,0,vn.size()));

			in.getline(buffer,255); // "endloop"
			in.getline(buffer,255); // "endfacet"

			in.getline(buffer,255); // get next buffer
		}*/
}




void ObjLoader::loadObj( const char *filename, unsigned short numGroups, char** groupNames )

{
	// Loading the OBJ file, matching it to the binary file

	std::ifstream in(filename);

	if(!in.good())
	{
		cerr << "ERROR: ObjLoader::loadObj(" << filename << ") file is not good" << endl;
		exit(-1);
	}

	char buffer[255], str[255];
	float f1,f2,f3;

	bool ignoreGroups = (numGroups == 0 );
	bool acceptGroup = true;
	unsigned char currentGroup = 0;

	while(!in.eof())
	{
		in.getline(buffer,255);

		// reading a vertex
		if (buffer[0]=='v' && (buffer[1]==' '  || buffer[1]==32) )
		{
			if ( sscanf(buffer,"v %f %f %f",&f1,&f2,&f3)==3)
			{
				v.push_back(Vector3(f1,f2,f3));
			}
			else
			{
				cerr << "ERROR: vertex not in wanted format in OBJLoader" << endl;
				cerr << buffer << endl;
				exit(-1);
			}
		}

		// reading vertex texture coordinat
		else if(buffer[0]=='v' && buffer[1]=='t')
		{
			if( sscanf(buffer,"vt %f %f ",&f1, &f2)==2)
			{
				//if (f1>0 && f2>0) // Why oh why was this neseccary it broke a sphere from 3dsmax and I spent a hot sommerday finding it !
				{
					vt.push_back(Vector2(f1,f2));
				}
			}
			else
			{
				cerr << "ERROR: vertex texture coordinate not in wanted format in OBJLoader" << endl;
				cerr << buffer << endl;
				exit(-1);
			}
		}

		// reading a vertex normal
		else if (buffer[0]=='v' && buffer[1]=='n')
		{
			if( sscanf(buffer,"vn %f %f %f",&f1, &f2, &f3)==3)
			{
				vn.push_back(Vector3(f1,f2,f3));
			}
			else
			{
				cerr << "ERROR: vertex normal not in wanted format in OBJLoader" << endl;
				cerr << buffer << endl;
				exit(-1);
			}
		}

		// reading faces 
		else if (buffer[0]=='f' && (buffer[1]==' ' || buffer[1]==32) )
		{
			if( ignoreGroups || acceptGroup ){
				//f v/vt/vn v/vt/vn v/vt/vn 
				Face f;
				int nt = sscanf(buffer,"f %d/%d/%d %d/%d/%d %d/%d/%d ",&f.v[0],&f.vt[0],&(f.vn[0]),&f.v[1],&f.vt[1],&(f.vn[1]),&f.v[2],&f.vt[2],&(f.vn[2]));
				if (nt!=9)
				{
					nt = sscanf(buffer,"f %d//%d %d//%d %d//%d ",&f.v[0],&f.vn[0],&f.v[1],&f.vn[1],&f.v[2],&f.vn[2]);
					if (nt!=6)
					{
						nt = sscanf(buffer,"f %d %d %d",&f.v[0],&f.v[1],&f.v[2]);
						if( nt!=3 )
						{
							cerr << "ERROR: I don't know the format of that face" << endl;
						}
					}
				}

				if( ignoreGroups )
					f.groupId = 0;
				else
					f.groupId = currentGroup;
				
				this->f.push_back(f);
			}
		}

		// reading a comment (or MELODY specific information)
		else if (buffer[0]=='#')
		{
			// a decal texture coordinat
			if ( sscanf(buffer,"#_#vt1 %f %f",&f1,&f2)==2  )
			{
				if(f1>0 && f2>0)
					vt1.push_back(Vector2(f1,f2));
			}
			// a tangent
			else if ( sscanf(buffer,"#_#tangent %f %f %f",&f1,&f2,&f3)==3  )
			{
				tangent.push_back(Vector3(f1,f2,f3));
			}
			// a binormal
			else if ( sscanf(buffer,"#_#binormal %f %f %f",&f1,&f2,&f3)==3  )
			{
				binormal.push_back(Vector3(f1,f2,f3));
			}
			else
			{
				// just a comment
			}
		}

		// new group
		else if (buffer[0]=='g')
		{
			str[0] = 0;
			sscanf( buffer,"g %s", str );
			acceptGroup = false;
			for( int i=0; i<numGroups; i++ )
				if( strcmp( str, groupNames[i] ) == 0 ){
					acceptGroup = true;
					currentGroup = i;
				}
		}

		else if (strlen(buffer)==0)
		{
		}

		else
		{
			/*cerr << "ObjLoader::ObjLoader("<<filename<<") encountered an unknown line:" << endl;
			cerr << "buffer: " << buffer << endl;
			cerr << "length: " << strlen(buffer) << endl;
			cerr << "first int: " << (int) buffer[0] << endl;*/
			//exit(-1);
		}
	}

	std::cout << "# f.size()        " << f.size() << endl;
	std::cout << "# v.size()        " << v.size() << endl;
	std::cout << "# vn.size()       " << vn.size() << endl;
	std::cout << "# vt.size()       " << vt.size() << endl;
	std::cout << "# vt1.size()      " << vt1.size() << endl;
	std::cout << "# tangent.size()  " <<  tangent.size() << endl;
	std::cout << "# binormal.size() " << binormal.size() << endl;
}

void ObjLoader::ThomasSmooth(float factor, int zLimit)
{
	if( zLimit < 0 ){
		SingleLaplacianSmooth( factor );
		return;
	}

	vector<Vector3> smoothingFactor;
	vector<float> numberOfSmoothing;

	for(unsigned int i=0; i<v.size(); i++)
	{
		smoothingFactor.push_back(Vector3(0,0,0));
		numberOfSmoothing.push_back(0);
	}

	vector<Face>::iterator fi;  // faces
	for(fi=f.begin();fi!=f.end();fi++)
	{
		Face f = *fi;
		for(int i=0;i<3;i++)
		{
			if( -1.0f*(v[f.v[i]-1].e[2]) >= zLimit ){

				smoothingFactor[ f.v[i]-1 ] = smoothingFactor[ f.v[i]-1 ] +
					v[ f.v[ ((i+1)%3 ) ]-1 ] - v[ f.v[i]-1 ];
				smoothingFactor[ f.v[i]-1 ] = smoothingFactor[ f.v[i]-1 ]+
					v[ f.v[ ((i+2)%3 ) ]-1 ] - v[ f.v[i]-1 ];
				numberOfSmoothing[ f.v[i]-1 ] = numberOfSmoothing[ f.v[i]-1 ]+
					2;
				Vector3 tmp(v[f.v[i]-1]);
				Vector3 tmp2(smoothingFactor[ f.v[i]-1 ]);
				float tmpf = -1.0f*(v[f.v[i]-1].e[2]);
				float ttt = tmpf;
			}	
		}
	}

	for(unsigned int i=0; i < v.size(); i++)
	{
		if( numberOfSmoothing[i]>0 )
			v[i] = v[i] + (smoothingFactor[i]/numberOfSmoothing[i])*factor;
	}
}

void ObjLoader::SingleLaplacianSmooth(float factor)
{
	vector<Vector3> smoothingFactor;
	vector<float> numberOfSmoothing;

	for(unsigned int i=0; i<v.size(); i++)
	{
		smoothingFactor.push_back(Vector3(0,0,0));
		numberOfSmoothing.push_back(0);
	}

	vector<Face>::iterator fi;  // faces
	for(fi=f.begin();fi!=f.end();fi++)
	{
		Face f = *fi;
		for(int i=0;i<3;i++)
		{
			smoothingFactor[ f.v[i]-1 ] = smoothingFactor[ f.v[i]-1 ] +
				v[ f.v[ ((i+1)%3 ) ]-1 ] - v[ f.v[i]-1 ];
			smoothingFactor[ f.v[i]-1 ] = smoothingFactor[ f.v[i]-1 ]+
				v[ f.v[ ((i+2)%3 ) ]-1 ] - v[ f.v[i]-1 ];
			numberOfSmoothing[ f.v[i]-1 ] = numberOfSmoothing[ f.v[i]-1 ]+
				2;
		}
	}

	for(unsigned int i=0; i < v.size(); i++)
	{
		v[i] = v[i] + (smoothingFactor[i]/numberOfSmoothing[i])*factor;
	}
}

void ObjLoader::TaubinSmooth(float growFactor, float iterations, float smooth  )
{
	for(int i=0;i<iterations;i++)
	{
		SingleLaplacianSmooth(smooth);
		SingleLaplacianSmooth(growFactor);
	}
	//grow(growFactor);
}

/*void ObjLoader::weld()
{
	weld( 0.1f );
}*/
/*
void ObjLoader::weld( float EPS )
{
	vector<Vector3> weldedV;	
	unsigned int nextBalance = 0;
	list<Vector3*>* found;

	KdTree_pfVec tree = KdTree_pfVec(v.size());

	for(unsigned int i=0; i < v.size(); i++)
	{
		found = tree.find(&v[i],EPS,1);
		if (found->begin() != found->end() ) // found 
		{
			list<Vector3*>::iterator it = found->begin();
			Vector3 *p = *it;

			v[i].e[0] = p->e[0];
			v[i].e[1] = -1;
			v[i].e[2] = -1;
		}
		else								// not found
		{
			tree.store(&v[i]);
			weldedV.push_back(Vector3(v[i]));
			v[i].e[0] = weldedV.size();
			v[i].e[1] = -1;
			v[i].e[2] = -1;

			if (i > nextBalance)
			{
				nextBalance = i+100;
				tree.balance();
			}
		}
	}

	for(unsigned int i=0; i < f.size(); i++)
	{
		for(unsigned int j=0; j<3; j++)
		{
			f[i].v[j] = (int) v[ f[i].v[j] -1 ].e[0];
		}
	}

	v.clear();
	v.resize(weldedV.size());

	for(unsigned int i=0; i < weldedV.size(); i++)
	{
		Vector3 t = weldedV[i];
		v[i] = t;
	}

	cout << "welded to " << v.size() << " nodes" << endl; 
}
*/
void ObjLoader::transform(float m4x4[][4])
{
	std::cout << "tranform with" << std::endl;
	std::cout << m4x4[0][0] << " " << m4x4[1][0] << " " << m4x4[2][0] << " " << m4x4[3][0] << std::endl;
	std::cout << m4x4[0][1] << " " << m4x4[1][1] << " " << m4x4[2][1] << " " << m4x4[3][1] << std::endl;
	std::cout << m4x4[0][2] << " " << m4x4[1][2] << " " << m4x4[2][2] << " " << m4x4[3][2] << std::endl;

	for(unsigned int i=0; i < v.size(); i++)
	{
		Vector3 t(0,0,0);
		for(int y=0;y<3;y++)
		{
			for(int x=0; x<4; x++)
			{
				//cout << x << " " << y << endl;
				if (x!=3)
				{
					t.e[y] += m4x4[x][y]*v[i].e[x];
					//	cout << m4x4[x][y] << " * " << v[i].e[x] << endl;
				}
				else
				{
					t.e[y] += m4x4[x][y]*1;
					//	cout << m4x4[x][y] << " * " << 1 << endl;
				}

			}
		}

		v[i]=t;
	}


	for(unsigned int i=0; i < vn.size(); i++)
	{
		Vector3 t(0,0,0);
		for(int y=0;y<3;y++)
		{
			for(int x=0; x<3; x++)
			{
				t.e[y] += m4x4[x][y]*vn[i].e[x];
				//	cout << m4x4[x][y] << " * " << v[i].e[x] << endl;


			}
		}

		vn[i]=t;
	}


}

void ObjLoader::computeTangentSpace()
{
	vector<int> triCount; 

	std::vector<Vector2> *tex;
	if(vt.size()>0)
	{
		tex = &vt;
	}
	else if (vt1.size()>0)
	{
		tex = &vt1;
	}
	else
	{
		std::cerr << "ERROR: ObjLoader::computeTangentSpace() no texture-coordinates" << std::endl;
		return;
	}


	/*//re-normalize winding of uv-triangles
	// calculate normals from binormal and tangent
	for(unsigned int i = 0; i < f.size(); i++)
	{
		Vector3 a((*tex)[f[i].vt[1]-1]-(*tex)[f[i].vt[0]-1]);
		Vector3 b((*tex)[f[i].vt[2]-1]-(*tex)[f[i].vt[0]-1]);
		a.makeUnitVector();
		b.makeUnitVector();
		if (a.dot(b)>0)
		{
			unsigned int t = f[i].vt[1];
			f[i].vt[1]=f[i].vt[0];
			f[i].vt[0]=t;
		}
	}*/


	triCount.resize(tex->size(),0);

	tangent.resize(tex->size());
	binormal.resize(tex->size());

	//http://www.netsoc.tcd.ie/~nash/tangent_note/tangent_note.html
	for(unsigned int i = 0; i < f.size(); i++)
	{
		int vt0 = f[i].vt[0]-1;
		int vt1 = f[i].vt[1]-1;
		int vt2 = f[i].vt[2]-1;

		if (vt0==-1 || vt1==-1 || vt2==-1)
		{
			cerr << "ERROR: a tex-coord is undefined" << endl;
		}

		int v0 = f[i].v[0]-1;
		int v1 = f[i].v[1]-1;
		int v2 = f[i].v[2]-1;

		float u13u = (*tex)[vt2][0] - (*tex)[vt0][0];
		float u13v = (*tex)[vt2][1] - (*tex)[vt0][1];
		float u12u = (*tex)[vt1][0] - (*tex)[vt0][0];
		float u12v = (*tex)[vt1][1] - (*tex)[vt0][1];

		Vector3 v1v2 = v[v1] - v[v0];
		Vector3 v1v3 = v[v2] - v[v0];

		float det_inv = 1.0f / (u12u * u13v - u13u * u12v);

		Vector3 q =  ( v1v2*u13v -  v1v3*u12v) * det_inv;
		Vector3 p =  ( v1v3*u12u -  v1v2*u13u ) * det_inv ;

		f[i].tangent[0]=vt0+1;
		f[i].tangent[1]=vt1+1;
		f[i].tangent[2]=vt2+1;

		tangent[vt0]+=p;
		binormal[vt0]+=q;

		tangent[vt1]+=p;
		binormal[vt1]+=q;

		tangent[vt2]+=p;
		binormal[vt2]+=q;

		int triCountSize = triCount.size();
		assert(vt0<triCountSize);
		assert(vt1<triCountSize);
		assert(vt2<triCountSize);

		triCount[vt0]++;
		triCount[vt1]++;
		triCount[vt2]++;



	}

	
	for(unsigned int i = 0; i < tangent.size(); i++)
	{
		if(triCount[i])
		{
			tangent[i] /= triCount[i];
			binormal[i] /= triCount[i];
		}
		tangent[i].makeUnitVector();
		binormal[i].makeUnitVector();
	}


	
	// calculate normals from binormal and tangent
/*	for(unsigned int i = 0; i < f.size(); i++)
	{
		for(unsigned int j=0; j<3; j++)
		{
			Vector3 newN = binormal[f[i].tangent[j]-1].cross(tangent[f[i].tangent[j]-1]);
			if(newN.dot(vn[f[i].vn[j]-1])<0)
			{
				binormal[f[i].tangent[j]-1] = binormal[f[i].tangent[j]-1]*-1;
			}
            vn[f[i].vn[j]-1] = newN;
		}

	}*/

	return;

}

void ObjLoader::computeStaticTangentspace()
{
//	this->weld();

	vn.resize(v.size(),Vector3(0,0,0));
	vector<int> triCount; triCount.resize(v.size(),0);

	tangent.resize(f.size());
	binormal.resize(f.size());


	//http://www.netsoc.tcd.ie/~nash/tangent_note/tangent_note.html
	for(unsigned int i = 0; i < f.size(); i++)
	{
		int vt0 = f[i].vt[0]-1;
		int vt1 = f[i].vt[1]-1;
		int vt2 = f[i].vt[2]-1;

		int v0 = f[i].v[0]-1;
		int v1 = f[i].v[1]-1;
		int v2 = f[i].v[2]-1;

		float u13u = vt[vt2][0] - vt[vt0][0];
		float u13v = vt[vt2][1] - vt[vt0][1];
		float u12u = vt[vt1][0] - vt[vt0][0];
		float u12v = vt[vt1][1] - vt[vt0][1];

		Vector3 v1v2 = v[v1] - v[v0];
		Vector3 v1v3 = v[v2] - v[v0];

		float det_inv = 1.0f / (u12u * u13v - u13u * u12v);

		Vector3 p =  ( v1v2*u13v -  v1v3*u12v) * det_inv;
		Vector3 q =  ( v1v3*u12u -  v1v2*u13u ) * det_inv ;

		tangent[i] = p;
		f[i].tangent[0]=i+1;
		f[i].tangent[1]=i+1;
		f[i].tangent[2]=i+1;
		binormal[i] = q; 

		// as many vertex normals as welded nodes
		f[i].vn[0] = f[i].v[0]; 
		f[i].vn[1] = f[i].v[1];
		f[i].vn[2] = f[i].v[2];

		for(int j=0;j<3;j++)
		{
			vn[f[i].vn[j]-1] = vn[f[i].vn[j]-1]+ cross(p,q);
			triCount[f[i].vn[j]-1]++;
		}


	}

	for(unsigned int i = 0; i < vn.size(); i++)
	{
		vn[i] = vn[i] / triCount[i];
		vn[i].makeUnitVector();
	}

	return;
}

void ObjLoader::translateVertices(Vector3 t)
{
	for(unsigned int i = 0; i < v.size(); i++)
	{
		v[i] = v[i] + t;
	}
}

void ObjLoader::scaleVertices(float s)
{
	for(unsigned int i = 0; i < v.size(); i++)
	{
		v[i] = v[i]*s;
	}

}

bool ObjLoader::sameSide(Vector2 p1, Vector2 p2, Vector2 a, Vector2 b)
{
/*	Vector3 cp1, cp2;

	cp1 = cross(b-a,p1-a);
	cp2 = cross(b-a,p2-a);

	if (dot(cp1,cp2) >= 0)
		return true;
	else*/
		return false;
}

bool ObjLoader::pointIntTriangle(Vector2 p, Vector2 a, Vector2 b, Vector2 c)
{
	if (sameSide(p,a, b,c) && sameSide(p,b, a,c) && sameSide(p,c, a,b)) 
		return true;
	else
		return false;
}


void ObjLoader::boundingBox(Vector3 &minV, Vector3&maxV)
{
	
	minV = v[0];
	maxV = v[0];

	for(unsigned int i = 1; i < v.size(); i++)
	{
		for(unsigned int j=0; j<3; j++)
		{
			minV[j] = std::min(v[i][j],minV[j]);
			maxV[j] = std::max(v[i][j],maxV[j]);
		}
	}


}





