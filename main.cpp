/*
    Title: Shoot 'Em All (V3)
    Author: Tashrif Ahmad Taif
*/

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <windows.h>
#include <mmsystem.h>
#include <bits/stdc++.h>
#include <chrono>

#ifndef BMPLOADER_H
#define BMPLOADER_H
#endif // BMPLOADER_H


#define ColorLiteWood {0.749, 0.565, 0.0}
#define ColorDarkWood {0.706, 0.373, 0.024}
#define ColorMidWood {0.902, 0.59, 0.22}
#define ColorExtraDarkWood {0.471, 0.25, 0.016}
#define ColorFanDark {0.0824, 0.3333, 0.3412}
#define ColorFanLite {0.1333, 0.612, 0.6313}
#define ColorWall {0.576471, 0.678431, 0.203922}
#define ColorDarkWall {0.450981, 0.541177, 0.129412}
#define Yellow 245.0 / 255, 255.0/255, 48.0/255
#define White 1.0, 1.0, 1.0
#define M_PI (2*acos(0.0))
#define double GLfloat
#define float GLfloat

using namespace std;

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
GLfloat ar = 16.0 / 9.0;
int windowHeight = 600;
int windowWidth = ar * windowHeight;

double txVal = 0, tyVal = 0;
bool flagRotate = true;
double rval = 0, fanSpeed = 0;
bool openDrawer = true;

bool eye = true, cen = false;
double sec, minute, hour;
bool testOpen;
bool lightLeft = false, lightUp = false, spotLight = false;
bool amb = true, diff = true, spec = true;
double ty = 0, tz = 0, tx = 0;
//double rx = 0, ry = 0, rz = 0;
double ds = 30;
GLfloat step = 1;
GLfloat stepAngle = 2.5;
const int nt = 30;				//number of slices along x-direction
const int ntheta = 10;
bool lamp = true;
bool shoot = false;
bool  Mainlight = true;

int curr = 0;

vector<bool> tarReach;
vector<int> cnt;
vector<unsigned int> ID;

float pS = 0.3;
float rS = 2;

auto start = std::chrono::high_resolution_clock::now();

class BmpLoader{
    public:
        unsigned char* textureData;
        int iWidth, iHeight;

        BmpLoader(const char* filename)
        {
            FILE *file=0;
            file=fopen(filename, "rb");
            if(!file)
                std::cout<<"File not found"<<std::endl;

            fread(&bfh, sizeof(BITMAPFILEHEADER),1,file);
            if(bfh.bfType != 0x4D42)
                std::cout<<"Not a valid bitmap"<<std::endl;

            fread(&bih, sizeof(BITMAPINFOHEADER),1,file);

            if(bih.biSizeImage==0)
                bih.biSizeImage=bih.biHeight*bih.biWidth*3;
            textureData = new unsigned char[bih.biSizeImage];
            fseek(file, bfh.bfOffBits, SEEK_SET);
            fread(textureData, 1, bih.biSizeImage, file);

            unsigned char temp;
            for(int i=0; i<bih.biSizeImage; i+=3)
            {
                temp = textureData[i];
                textureData[i] = textureData[i+2];
                textureData[i+2] = temp;

            }

            iWidth = bih.biWidth;
            iHeight = bih.biHeight;

            fclose(file);
        }

        ~BmpLoader()
        {
            delete [] textureData;
        }

    private:
        BITMAPFILEHEADER bfh;
        BITMAPINFOHEADER bih;
};

class point{
public:
	GLfloat x,y,z;
	point(){
        x = 0, y = 0, z = 0;
	}
    point(GLfloat X, GLfloat Y, GLfloat Z){
        x = X, y = Y, z = Z;
    }

    point moveIn(vector<float> v, float d){
        return point(x + v[0] * d, y + v[1] * d, z + v[2] * d);
    }

    float dist(point b){
        return sqrt((x - b.x) * (x - b.x) + (y - b.y) * (y - b.y) + (z - b.z) * (z - b.z));
    }
};

vector<point> tar;

class myVector{
public:
	GLfloat x, y, z;

	myVector(){
	}

	myVector(double px, double py, double pz){
		x = px, y = py, z = pz;
	}

	void copyIt(myVector a){
		x = a.x, y = a.y, z = a.z;
	}

    bool notEqual(myVector a){
        return (a.x != x || a.y != y || a.z != z);
    }

	myVector cross(myVector a){
		return myVector(y*a.z - z*a.y, z*a.x - x*a.z, x*a.y - y*a.x);
	}

	GLfloat dot(myVector a){
        return a.x * x + a.y * y + a.z * z;
	}

	myVector add(myVector a){
		return myVector(x + a.x, y + a.y, z + a.z);
	}

	myVector multiply(double val){
		return myVector(val*x, val*y, val*z);
	}

	void Move(myVector to, float steps){
		x += (to.x)*steps;
		y += (to.y)*steps;
		z += (to.z)*steps;
	}

	void Rotate(myVector per, float angle){
		myVector t = cross(per);
		myVector m = *this;
		m = m.multiply(cos(angle*M_PI/180.0));
		t = t.multiply(sin(angle*M_PI/180.0));
		m = t.add(m);
		copyIt(m);
	}

	GLfloat getVal(){
        return sqrt(x * x + y * y + z * z);
	}

	void makeUnitVec(){
        GLfloat val = getVal();
        x /= val;
        y /= val;
        z /= val;
	}

	GLfloat angleWith(myVector a){
        GLfloat A = a.getVal(), B = this->getVal();

        return acos(this->dot(a) / (A * B)) * 180 / M_PI;
	}

	GLfloat projectionRad(point A, point B){ //projection distance of point B with A in the direction of this vector
	    myVector AB = myVector(B.x - A.x, B.y - A.y, B.z - A.z);

	    GLfloat Bhumi = dot(AB);
	    GLfloat Otibhuj = AB.getVal();

	    return sqrt(Otibhuj * Otibhuj - Bhumi * Bhumi);
	}

	GLfloat projectionDist(point A, point B){ //distance from A to projection point of B in direction of this vector
	    myVector AB = myVector(B.x - A.x, B.y - A.y, B.z - A.z);

	    return dot(AB);
	}
};

myVector EYE(50, 1, 4.5), LOOK(0, 1, 0), UP(0,0,1), RIGHT(1, 0, 0);

vector<myVector> tarV;

void animate(){
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    sec = -((elapsed.count() * 6));
    minute = -((elapsed.count() / 60 * 6));
    hour = -((elapsed.count() / 3600 * 6));

    if(sec <= -360) sec += 360;
    if(minute <= -360) minute += 360;
    if(hour <= -360) hour += 360;

	glutPostRedisplay();
}

static void resize(int width, int height){
    GLfloat nheight = height;
    GLfloat nwidth = ar * height;
    glViewport((width - nwidth) / 2, (height - nheight) / 2, nwidth, nheight);
}

vector<GLfloat> getColor(GLfloat x, GLfloat y, GLfloat z){
    return {x / 255.0, y / 255.0, z / 255.0};
}

vector<vector<GLfloat>> generatePoints8p(float x, float y, float z){
    vector<vector<GLfloat>> points(8, vector<GLfloat> (3));
    points = {
        {0.0, 0.0, 0.0}, // 0
        {x, 0.0, 0.0}, // 1
        {x, y, 0.0}, // 2
        {0.0, y, 0.0}, // 3
        {0.0, 0.0, z}, // 4
        {x, 0.0, z}, // 5
        {x, y, z}, // 6
        {0.0, y, z} // 7
    };

    return points;
}

vector<vector<GLubyte>> getMyIndicesForCube(){
    vector<vector<GLubyte>> ind(6, vector<GLubyte> (4));

    ind = {
        {1, 0, 3, 2},  // X-Y-down
        {6, 7, 4, 5},  // X-Y-up
        {5, 4, 0, 1},  // X-Z-front
        {2, 3, 7, 6},  // X-Z-back
        {2, 6, 5, 1},  // Y-Z-right
        {7, 3, 0, 4},  // Y-Z-left
    };

    return ind;
}

static void getNormal3p(GLfloat x1, GLfloat y1,GLfloat z1, GLfloat x2, GLfloat y2,GLfloat z2, GLfloat x3, GLfloat y3,GLfloat z3){
    GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

    Ux = x2-x1;
    Uy = y2-y1;
    Uz = z2-z1;

    Vx = x3-x1;
    Vy = y3-y1;
    Vz = z3-z1;

    Nx = Uy*Vz - Uz*Vy;
    Ny = Uz*Vx - Ux*Vz;
    Nz = Ux*Vy - Uy*Vx;

    glNormal3f(Nx,Ny,Nz);
}

void colorMaterial(vector<GLfloat> color){
    GLfloat no_mat[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat mat_ambient[] = { color[0] * 0.5, color[1] * 0.5, color[2] * 0.5, 1.0 };
    GLfloat mat_diffuse[] = { color[0], color[1], color[2], 1.0 };
    GLfloat mat_specular[] = { color[0], color[1], color[2], 1.0 };
    GLfloat mat_shininess[] = {30};
    GLfloat mat_em[] = {0.0, 0.0, 0.0, 1.0};

    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess);
    glMaterialfv(GL_FRONT, GL_EMISSION, mat_em);
}

void colorMaterialEmit(vector<GLfloat> color){
    GLfloat no_mat[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat mat_ambient[] = { color[0] * 0.5, color[1] * 0.5, color[2] * 0.5, 1.0};
    GLfloat mat_diffuse[] = { color[0], color[1], color[2], 1.0};
    GLfloat mat_specular[] = { color[0], color[1], color[2], 1.0};
    GLfloat em[] = {color[0], color[1], color[2], 1.0};
    GLfloat mat_shininess[] = {30};

    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess);
    glMaterialfv( GL_FRONT, GL_EMISSION, em);
}

void drawCube(float x, float y, float z, vector<GLfloat> color){
    vector<vector<GLfloat>> points = generatePoints8p(x, y, z);
    vector<vector<GLubyte>> ind = getMyIndicesForCube();

    colorMaterial(color);

    glBegin(GL_QUADS);

        for(int i = 0; i < 6; i++){

            getNormal3p(
                        points[ind[i][0]][0], points[ind[i][0]][1], points[ind[i][0]][2],
                        points[ind[i][1]][0], points[ind[i][1]][1], points[ind[i][1]][2],
                        points[ind[i][2]][0], points[ind[i][2]][1], points[ind[i][2]][2]
            );

            for(int j = 0; j < 4; j++){
                glVertex3fv(&points[ind[i][j]][0]);
            }
        }

    glEnd();
}

void drawSphere(double radius,int slices,int stacks, vector<GLfloat> color){

    colorMaterial(color);

	point points[slices + 1][stacks + 1];

	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(M_PI/2));
		r=radius*cos(((double)i/(double)stacks)*(M_PI/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*M_PI);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*M_PI);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		//glColor3f(1, 1, 1);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{

			    getNormal3p(points[i][j].x,points[i][j].y,points[i][j].z,
                    points[i][j+1].x,points[i][j+1].y,points[i][j+1].z,
                    points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);

			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);

				getNormal3p(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z,
                    points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z,
                    points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);

                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}glEnd();
		}
	}
}

void drawSphere(double radius,int slices,int stacks, vector<GLfloat> color, bool glow){

    colorMaterialEmit(color);

	point points[slices + 1][stacks + 1];

	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(M_PI/2));
		r=radius*cos(((double)i/(double)stacks)*(M_PI/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*M_PI);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*M_PI);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		//glColor3f(1, 1, 1);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{

			    getNormal3p(points[i][j].x,points[i][j].y,points[i][j].z,
                    points[i][j+1].x,points[i][j+1].y,points[i][j+1].z,
                    points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);

			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);

				getNormal3p(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z,
                    points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z,
                    points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);

                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}glEnd();
		}
	}
}

void drawHalfSphere(double radius,int slices,int stacks, vector<GLfloat> color, bool glow = false){
    if(!glow) colorMaterial(color);
    else colorMaterialEmit(color);

	struct point points[slices + 1][stacks + 1];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(M_PI/2));
		r=radius*cos(((double)i/(double)stacks)*(M_PI/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*M_PI);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*M_PI);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		//glColor3f(1, 1, 1);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);

			    getNormal3p(points[i][j].x,points[i][j].y,points[i][j].z,
                    points[i][j+1].x,points[i][j+1].y,points[i][j+1].z,
                    points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);

			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
			glEnd();
		}
	}
}

void drawCylinder(float radius, float height, int hSlices, int vSlices, vector<GLfloat> color){
	colorMaterial(color);

	vector<vector<point>> points(hSlices+2, vector<point>(vSlices+1));

	for(int j = 0; j <= vSlices; j++){
		double angle = (2*M_PI) * j/(double) vSlices;
		points[0][j].x = radius * cos(angle);
		points[0][j].y = radius * sin(angle);
		points[0][j].z = 0;
	}

	for(int j = 0; j <= vSlices; j++){
        glBegin(GL_TRIANGLES);
            getNormal3p(
                        0, 0, 0,
                        points[0][j].x, points[0][j].y, points[0][j].z,
                        points[0][(j + 1) % (vSlices + 1)].x, points[0][(j + 1) % (vSlices + 1)].y, points[0][(j + 1) % (vSlices + 1)].z
                        );
            glVertex3f(0, 0, 0);
            glVertex3f(points[0][j].x, points[0][j].y, points[0][j].z);
            glVertex3f(points[0][(j + 1) % (vSlices + 1)].x, points[0][(j + 1) % (vSlices + 1)].y, points[0][(j + 1) % (vSlices + 1)].z);

            getNormal3p(
                        0, 0, height,
                        points[0][(j + 1) % (vSlices + 1)].x, points[0][(j + 1) % (vSlices + 1)].y, height,
                        points[0][j].x, points[0][j].y, height
                        );

            glVertex3f(0, 0, height);
            glVertex3f(points[0][(j + 1) % (vSlices + 1)].x, points[0][(j + 1) % (vSlices + 1)].y, height);
            glVertex3f(points[0][j].x, points[0][j].y, height);
        glEnd();
	}

	for(int i = 1; i <= hSlices; i++){
		for(int j = 0; j <= vSlices; j++){
			points[i][j].x = points[0][j].x;
			points[i][j].y = points[0][j].y;
			points[i][j].z = i * (height / hSlices);
		}
	}

	for(int i = 0; i < hSlices; i++){
        for(int j = 0; j <= vSlices; j++){
            glBegin(GL_QUADS);
                getNormal3p(
                            points[i][(j + 1) % (vSlices + 1)].x, points[i][(j + 1) % (vSlices + 1)].y, points[i][(j + 1) % (vSlices + 1)].z,
                            points[i + 1][(j + 1) % (vSlices + 1)].x, points[i + 1][(j + 1) % (vSlices + 1)].y, points[i + 1][(j + 1) % (vSlices + 1)].z,
                            points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z
                            );

                glVertex3f(points[i][(j + 1) % (vSlices + 1)].x, points[i][(j + 1) % (vSlices + 1)].y, points[i][(j + 1) % (vSlices + 1)].z);
                glVertex3f(points[i + 1][(j + 1) % (vSlices + 1)].x, points[i + 1][(j + 1) % (vSlices + 1)].y, points[i + 1][(j + 1) % (vSlices + 1)].z);
                glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z);
                glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
            glEnd();
        }
	}
}

void drawDrawer(float len, float width, float height, float thin, float Handlegap){
    //glBegin(GL_QUAD_STRIP);
        vector<GLfloat> color = ColorDarkWood;
        glPushMatrix();
            glTranslatef(0.0, 0.0, thin);
            drawCube(thin, width, height, color); // side Left
        glPopMatrix();

        glPushMatrix();
            glTranslatef(len - thin, 0.0, thin);
            drawCube(thin, width, height, color); // side Right
        glPopMatrix();

        color.clear();
        color = ColorLiteWood;

        glPushMatrix();
            glTranslatef(thin, 0.0, thin);
            drawCube(len - 2 * thin, thin, height - thin, color); // front
        glPopMatrix();

        glPushMatrix();
            glTranslatef(thin, width - thin, thin);
            drawCube(len - 2 * thin, thin, height - thin, color); // back
        glPopMatrix();

        color = ColorMidWood;

        glPushMatrix();
            drawCube(len, width, thin, color);  // floor
        glPopMatrix();

        color = ColorExtraDarkWood;

        glPushMatrix();
            glTranslatef(len, -2.0 * thin, 0.0);
            drawCube(2 * thin, width + 4 * thin, height + 2 * thin, color);  // Main Plate
        glPopMatrix();

        // Handle
        glPushMatrix();
            glTranslatef(len + 2 * thin, width / 3, 3 * height / 7);
            drawCube(Handlegap, thin, height / 7, color);  // Handle side front
        glPopMatrix();

        color = ColorDarkWood;

        glPushMatrix();
            glTranslatef(len + 2 * thin, 2 * width / 3 - thin, 3 * height / 7);
            drawCube(Handlegap, thin, height / 7, color);  // Handle side Back
        glPopMatrix();

        glPushMatrix();
            glTranslatef(len + 2 * thin + Handlegap, width / 3 - thin, 3 * height / 7 - thin);
            drawCube(1.5 * thin, width / 3 + 2 * thin, height / 7 + 2 * thin, color);  // Handle plate
        glPopMatrix();

    //glEnd();
}

void DrawHalfBedWithDrawer(float len, float width, float heigth, bool Open, float HandleGap = 0.1){
    float drawerLen = 2 * width / 3;
    float drawerWidth = 2 * len / 3;
    float drawerHeigth = 2 * heigth / 3;
    float drawerThin = drawerWidth / 80.0;
    vector<GLfloat> color = ColorDarkWood;

    glPushMatrix();
        glTranslatef( width / 3 + (0.5 * Open), len / 6, heigth / 6);
        drawDrawer(drawerLen, drawerWidth, drawerHeigth, drawerThin, HandleGap);
    glPopMatrix();

    glPushMatrix();
        drawCube(width, len, heigth / 6, color); // down - floor
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0.0, 0.0, 5 * heigth / 6);
        drawCube(width, len, heigth / 6, color); // top - floor
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0.0, 0.0, heigth / 6);
        drawCube(width, len / 6, 2 * heigth / 3, color); // front side
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0.0, 5 * len / 6, heigth / 6);
        drawCube(width, len / 6, 2 * heigth / 3, color); // front side
    glPopMatrix();
}

void drawFlatPyramid(float dLen, float dWidth, float extraWidth, vector<GLfloat> colorNotUsed, vector<GLfloat> color){
    vector<vector<GLfloat>> points = generatePoints8p(dLen, dWidth, extraWidth);
    vector<vector<GLubyte>> ind = getMyIndicesForCube();

    points[4][0] -= extraWidth; points[4][1] -= extraWidth;
    points[5][0] += extraWidth; points[5][1] -= extraWidth;
    points[6][0] += extraWidth; points[6][1] += extraWidth;
    points[7][0] -= extraWidth; points[7][1] += extraWidth;

    GLfloat no_mat[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat mat_ambient[] = { color[0] * 0.5, color[1] * 0.3, color[2] * 0.3, 1.0 };
    GLfloat mat_diffuse[] = { color[0], color[1], color[2], 1.0 };
    GLfloat mat_specular[] = { color[0], color[1], color[2], 1.0 };
    GLfloat mat_shininess[] = {30};

    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess);

    glBegin(GL_QUADS);

        for(int i = 0; i < 6; i++){

            getNormal3p(points[ind[i][0]][0], points[ind[i][0]][1], points[ind[i][0]][2],
                    points[ind[i][1]][0], points[ind[i][1]][1], points[ind[i][1]][2],
                    points[ind[i][2]][0], points[ind[i][2]][1], points[ind[i][2]][2]);

            for(int j = 0; j < 4; j++){
                glVertex3fv(&points[ind[i][j]][0]);
            }
        }

    glEnd();
}

void drawFlatPyramidTex(float dLen, float dWidth, float extraWidth, vector<GLfloat> colorNotUsed, vector<GLfloat> color){
    vector<vector<GLfloat>> points = generatePoints8p(dLen, dWidth, extraWidth);
    vector<vector<GLubyte>> ind = getMyIndicesForCube();

    points[4][0] -= extraWidth; points[4][1] -= extraWidth;
    points[5][0] += extraWidth; points[5][1] -= extraWidth;
    points[6][0] += extraWidth; points[6][1] += extraWidth;
    points[7][0] -= extraWidth; points[7][1] += extraWidth;

    GLfloat no_mat[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat mat_ambient[] = { color[0] * 0.5, color[1] * 0.3, color[2] * 0.3, 1.0 };
    GLfloat mat_diffuse[] = { color[0], color[1], color[2], 1.0 };
    GLfloat mat_specular[] = { color[0], color[1], color[2], 1.0 };
    GLfloat mat_shininess[] = {30};

    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess);

    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glBindTexture(GL_TEXTURE_2D, ID[1]);
    glBegin(GL_QUADS);

        for(int i = 0; i < 6; i++){

            getNormal3p(
                        points[ind[i][0]][0], points[ind[i][0]][1], points[ind[i][0]][2],
                        points[ind[i][1]][0], points[ind[i][1]][1], points[ind[i][1]][2],
                        points[ind[i][2]][0], points[ind[i][2]][1], points[ind[i][2]][2]
            );

            for(int j = 0; j < 4; j++){
                glVertex3fv(&points[ind[i][j]][0]);
                if(j == 0) glTexCoord2f(1,0);
                if(j == 1) glTexCoord2f(1,1);
                if(j == 2) glTexCoord2f(0,1);
                if(j == 3) glTexCoord2f(0,0);
            }
        }

    glEnd();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
}

void drawFlatPyramid(float dLen, float dWidth,float dheight, float extraWidth, vector<GLfloat> colorNotUsed, vector<GLfloat> color){
    vector<vector<GLfloat>> points = generatePoints8p(dLen, dWidth, dheight);
    vector<vector<GLubyte>> ind = getMyIndicesForCube();

    points[4][0] -= extraWidth; points[4][1] -= extraWidth;
    points[5][0] += extraWidth; points[5][1] -= extraWidth;
    points[6][0] += extraWidth; points[6][1] += extraWidth;
    points[7][0] -= extraWidth; points[7][1] += extraWidth;

    GLfloat no_mat[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat mat_ambient[] = { color[0] * 0.5, color[1] * 0.3, color[2] * 0.3, 1.0 };
    GLfloat mat_diffuse[] = { color[0], color[1], color[2], 1.0 };
    GLfloat mat_specular[] = { color[0], color[1], color[2], 1.0 };
    GLfloat mat_shininess[] = {30};

    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess);

    glBegin(GL_QUADS);

        for(int i = 0; i < 6; i++){

            getNormal3p(points[ind[i][0]][0], points[ind[i][0]][1], points[ind[i][0]][2],
                    points[ind[i][1]][0], points[ind[i][1]][1], points[ind[i][1]][2],
                    points[ind[i][2]][0], points[ind[i][2]][1], points[ind[i][2]][2]);

            for(int j = 0; j < 4; j++){
                glVertex3fv(&points[ind[i][j]][0]);
            }
        }

    glEnd();
}

void drawFlatPyramid(float dLen, float dWidth,float dheight, float extraWidth, vector<GLfloat> color){
    vector<vector<GLfloat>> points = generatePoints8p(dLen, dWidth, dheight);
    vector<vector<GLubyte>> ind = getMyIndicesForCube();

    points[4][0] -= extraWidth; points[4][1] -= extraWidth;
    points[5][0] += extraWidth; points[5][1] -= extraWidth;
    points[6][0] += extraWidth; points[6][1] += extraWidth;
    points[7][0] -= extraWidth; points[7][1] += extraWidth;

    colorMaterial(color);

    glBegin(GL_QUADS);

        for(int i = 0; i < 6; i++){

            getNormal3p(points[ind[i][0]][0], points[ind[i][0]][1], points[ind[i][0]][2],
                    points[ind[i][1]][0], points[ind[i][1]][1], points[ind[i][1]][2],
                    points[ind[i][2]][0], points[ind[i][2]][1], points[ind[i][2]][2]);

            for(int j = 0; j < 4; j++){
                glVertex3fv(&points[ind[i][j]][0]);
            }
        }

    glEnd();
}

void drawBalish(float len, float width, float height){
    vector<GLfloat> color;
    color = {0.6, 0, 0};

    glPushMatrix();
        glTranslatef(0.2, 0, 0);
        drawCube(len - 0.4, width, height, color);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0,0.1,0.1);
        glRotatef(90, 0,1,0);
        glRotatef(90, 0, 0, 1);
        drawFlatPyramid(width - 0.2, height - 0.2, 0.2, 0.1, color);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(len, width - 0.1,0.1);
        glRotatef(-90, 0,1,0);
        glRotatef(-90, 0, 0, 1);
        drawFlatPyramid(width - 0.2, height - 0.2, 0.2, 0.1, color);
    glPopMatrix();
}

void drawBed(float len, float width, float height){
    vector<GLfloat> color = ColorExtraDarkWood;

    glPushMatrix();
        DrawHalfBedWithDrawer(len / 2, width, height, 0);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0.0, len / 2, 0.0);
        DrawHalfBedWithDrawer(len / 2, width, height, 0);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0.0, len, 0.0);
        drawCube(width, 0.2, 3 * height, color);
    glPopMatrix();

    color = {230.0 / 255, 235.0 / 255, 234.0 / 255};

    glPushMatrix(); // toshok
        glTranslatef(0.04, 0.04, height);
        drawCube(width - 0.08, len - 0.08, 0.3, color);
    glPopMatrix();

    float blen = width / 2 - (width / 10);
    float bwidth = 2 * blen / 3;
    float bheight = 0.3;

    glPushMatrix();
        glTranslatef(width / 20, len - bwidth - width / 20, height + 0.3);
        drawBalish(blen, bwidth, bheight);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(width / 2 + width / 20, len - bwidth - width / 20, height + 0.3);
        drawBalish(blen, bwidth, bheight);
    glPopMatrix();
}

void drawWardrobe(float length, float width, float heigth, int NumberOfDrawer, float extraWidth){
    float dHeigth = (length - extraWidth) / NumberOfDrawer;
    float dLen = length;
    float dWidth = width;

    for(int i = 0; i < NumberOfDrawer; i++){
        glPushMatrix();
            glTranslatef(0, 0, i * dHeigth);
            DrawHalfBedWithDrawer(dLen, dWidth, dHeigth, 0);
        glPopMatrix();
    }

    glPushMatrix();
        glTranslatef(0,0, NumberOfDrawer * dHeigth);
        drawFlatPyramid(dWidth, dLen, extraWidth, ColorDarkWood, ColorExtraDarkWood);
    glPopMatrix();
}

void drawPizzaSlice(float len, float side, float heigth){
    vector<vector<GLfloat>> points(6, vector<GLfloat>(3));

    points = {
        {0, 0, 0},
        {side, len, 0},
        {-side, len, 0},
        {0, 0, heigth},
        {side, len, heigth},
        {-side, len, heigth}
    };

    vector<vector<GLubyte>> TriInd(2, vector<GLubyte> (3));
    TriInd = {
        {0, 2, 1},
        {3, 4, 5}
    };

    vector<vector<GLubyte>> RectInd(3, vector<GLubyte> (4));
    RectInd = {
        {0, 3, 5, 2},
        {0, 1, 4, 3},
        {1, 2, 5, 4}
    };

    vector<GLfloat> color = ColorFanDark;

    GLfloat no_mat[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat mat_ambient[] = { color[0] * 0.5, color[1] * 0.3, color[2] * 0.3, 1.0 };
    GLfloat mat_diffuse[] = { color[0], color[1], color[2], 1.0 };
    GLfloat mat_specular[] = { color[0], color[1], color[2], 1.0 };
    GLfloat mat_shininess[] = {30};

    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess);

    glBegin(GL_QUADS);
        for(int i = 0; i < 3; i++){

            getNormal3p(points[RectInd[i][0]][0], points[RectInd[i][0]][1], points[RectInd[i][0]][2],
                    points[RectInd[i][1]][0], points[RectInd[i][1]][1], points[RectInd[i][1]][2],
                    points[RectInd[i][2]][0], points[RectInd[i][2]][1], points[RectInd[i][2]][2]);

            for(int j = 0; j < 4; j++){
                glVertex3fv(&points[RectInd[i][j]][0]);
            }
        }

    glEnd();

    glBegin(GL_TRIANGLES);
        for(int i = 0; i < 2; i++){

            getNormal3p(points[TriInd[i][0]][0], points[TriInd[i][0]][1], points[TriInd[i][0]][2],
                    points[TriInd[i][1]][0], points[TriInd[i][1]][1], points[TriInd[i][1]][2],
                    points[TriInd[i][2]][0], points[TriInd[i][2]][1], points[TriInd[i][2]][2]);

            for(int j = 0; j < 3; j++){
                glVertex3fv(&points[TriInd[i][j]][0]);
            }
        }

    glEnd();
}

void drawRoundPlate(float radius, float heigth, int numberOfTriangle = 360){
    float theta = 2.0 * M_PI / numberOfTriangle;
    float len = radius * cos(theta / 2.0);
    float side = radius * sin(theta / 2.0);

    for(int i = 0; i < numberOfTriangle; i++){
        glPushMatrix();
            glRotatef(i * theta * 180.0 / M_PI, 0, 0, 1);
            drawPizzaSlice(len, side, heigth);
        glPopMatrix();
    }
}

void drawBlade(float length, float width, float heigth = 0.05){
    float outLength = 7 * length / 8;
    float inLength = length / 8  + 2 * length / 10;
    float inWidth = width / 5;

    vector<GLfloat> color = ColorFanLite;

    glPushMatrix();
        drawCube(width, outLength, heigth, color);
    glPopMatrix();

    color = ColorFanDark;

    glPushMatrix();
        glTranslatef(inWidth, -inLength + (length / 10), heigth);
        drawCube(inWidth, inLength, heigth, color);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(3 * inWidth, -inLength + (length / 10), heigth);
        drawCube(inWidth, inLength, heigth, color);
    glPopMatrix();
}

void drawFan(float longHeigth, float shortHeigth, float radius, float bladeWidth, float numberOfBlades = 3){
    glPushMatrix();
        glTranslatef(0, 0, shortHeigth);
        drawRoundPlate(radius / 32, longHeigth);
    glPopMatrix();

    glPushMatrix();
        glRotatef(fanSpeed, 0, 0, 1);
        //glTranslatef()
        glPushMatrix();
            drawRoundPlate(radius / 4, shortHeigth);
        glPopMatrix();

        glPushMatrix();
            glTranslatef(-(bladeWidth / 2), -(radius / 3 - 3 * radius / 40), shortHeigth);

            float theta = 360 / numberOfBlades;

            for(int i = 0; i < numberOfBlades; i++){
                glPushMatrix();
                    glTranslatef( bladeWidth / 2, +(radius / 3 - 3 * radius / 40), -shortHeigth);
                    glRotatef(i * theta, 0, 0, 1);
                    glTranslatef( -bladeWidth / 2, (radius / 3 - 3 * radius / 40), shortHeigth);
                    drawBlade(3 * radius / 4, bladeWidth);
                glPopMatrix();
            }

        glPopMatrix();
    glPopMatrix();
}

void drawClock(float len){
    float extraWidth = 0.1;
    len -= 2 * extraWidth;
    glPushMatrix(); //Whole Clock
        glTranslatef(len / 2, len / 2, 0);
        glRotatef(90, 0, 0, 1);

        glPushMatrix(); // Frame
            glTranslatef(len / 2, -len / 2, 0);

            glPushMatrix(); // Nicher Palla
                glRotatef(180, 0, 1, 0);
                drawFlatPyramid(len, len, extraWidth, ColorDarkWood, {1,1,1});
            glPopMatrix();

            glPushMatrix(); //Border
                //glTranslatef(-len, extraWidth, 0);
                glRotatef(90, 0, 1, 0);
                glRotatef(90, 0, 0, 1);
                drawFlatPyramid(len, 0.1, extraWidth, ColorLiteWood, ColorExtraDarkWood);
            glPopMatrix();

            glPushMatrix(); //Border
                glTranslatef(-len, 0, extraWidth);
                glRotatef(180, 0, 1, 0);
                glRotatef(90, 0, 1, 0);
                glRotatef(90, 0, 0, 1);
                drawFlatPyramid(len, 0.1, extraWidth, ColorLiteWood, ColorExtraDarkWood);
            glPopMatrix();

            glPushMatrix(); //Border
                glTranslatef(0, len, 0);
                glRotatef(90, 0, 0, 1);
                glPushMatrix();
                    //glTranslatef(-len, extraWidth, 0);
                    glRotatef(90, 0, 1, 0);
                    glRotatef(90, 0, 0, 1);
                    drawFlatPyramid(len, 0.1, extraWidth, ColorLiteWood, ColorExtraDarkWood);
                glPopMatrix();

                glPushMatrix();
                    glTranslatef(-len, 0, extraWidth);
                    glRotatef(180, 0, 1, 0);
                    glRotatef(90, 0, 1, 0);
                    glRotatef(90, 0, 0, 1);
                    drawFlatPyramid(len, 0.1, extraWidth, ColorLiteWood, ColorExtraDarkWood);
                glPopMatrix();
            glPopMatrix();
        glPopMatrix();

        float frameLen = len;
        len -= 4 * extraWidth;
        len /= 6;

        vector<GLfloat> color = {0,0,0}; //black

        glPushMatrix(); // second-stick
            glRotatef(sec, 0, 0, 1);
            glTranslatef(0, -len / 24, extraWidth);
            drawCube(3 * len, len / 12, extraWidth / 3, color);
        glPopMatrix();

        glPushMatrix(); // Minute-stick
            glRotatef(minute, 0, 0, 1);
            glTranslatef(0, -len / 16, extraWidth + extraWidth / 3);
            drawCube(2 * len, len / 8, extraWidth / 3, color);
        glPopMatrix();

        glPushMatrix(); //hour - stick
            glRotatef(hour, 0, 0, 1);
            glTranslatef(0, -len / 8, extraWidth + 2*extraWidth / 3);
            drawCube(1.5 * len, len / 4, extraWidth / 3, color);
        glPopMatrix();

        glPushMatrix(); //knob
            glTranslatef( 0, 0, extraWidth + 3*extraWidth / 3);
            drawRoundPlate(len / 2.5, 0.01);
        glPopMatrix();

        color = ColorExtraDarkWood;

        for(int i = 0; i < 12; i++){
            glPushMatrix();
                glRotatef(i * 30, 0, 0, 1);
                glTranslatef(3 * len - 0.1, 0, extraWidth);
                drawCube(0.1, 0.1, 0.02, color);
            glPopMatrix();
        }

    glPopMatrix();
}

void drawBoxThatCanOpen(float len, float width, float height, bool openDoor = false){
    vector<GLfloat> color;

    color = {64.0 / 255, 28.0 / 255, 1.0 / 255};

    float palla = len / 10;

    glPushMatrix();
        drawCube(len, width, palla, color);
    glPopMatrix();

    color = {38.0/255, 17.0/255, 2.0/255};

    glPushMatrix();
        glTranslatef(0,0,palla);
        drawCube(palla, width, height - 2 * palla, color);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(9 * palla, 0, palla);
        drawCube(palla, width, height - 2 * palla, color);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(palla, width - palla, palla);
        drawCube(len - 2 * palla, palla, height - 2 * palla, color);
    glPopMatrix();

    //color[1] = {38.0/255, 17.0/255, 2.0/255};
    glPushMatrix();
        glTranslatef(0,0, height - palla);
        drawCube(len, width, palla, color);
    glPopMatrix();

    color.clear();

    color = {168.0/255, 122.0/255, 22.0/255};

    glPushMatrix();
        if(openDoor){
            glTranslatef(palla, 0.05, palla);
            glRotatef(-30, 0, 0, 1);
            glTranslatef(-palla, -0.05, -palla);
        }

        glTranslatef(palla, 0.05, palla);
        drawCube(len - 2 * palla, palla / 4, height - 2 * palla, color);

        glPushMatrix();
            if(openDoor){
                glTranslatef(len - palla - palla - 0.05 - palla, -0.05, height / 2);
                glRotatef(-85, 0, 0, 1);
                glTranslatef(-(len - palla - palla - 0.05 - palla), 0.05, -height / 2);
            }

            glTranslatef(len - palla - palla - 0.1 - palla, -0.05, height / 2);
            drawCube(palla / 2, palla / 2, palla / 2, {0.1, 0.1, 0.1});
        glPopMatrix();
    glPopMatrix();

//    color = {(float)(rand() % 255)/255, (float)(rand() % 255)/255, (float)(rand() % 255)/255};
//
//    glPushMatrix();
//        glTranslatef(2 * palla, palla, palla);
//        drawCube(len / 2, width / 2, height / 2, color);
//    glPopMatrix();
}

void drawDressingTable(float len, float width, float height, int numberOfDrawer){
    float dlen = len / numberOfDrawer;
    float dheight = height / 3;

    for(int i = 0; i < numberOfDrawer; i++){
        glPushMatrix();
            glTranslatef(i * dlen, 0, 0);
            drawBoxThatCanOpen(dlen, width, dheight, 0);
        glPopMatrix();
    }

    float palla = width / 6;

    vector<GLfloat> color = {89.0 / 255, 38.0 / 255, 8.0 / 255};

    float mirrorWidth = palla / 6;

    glPushMatrix();
        glTranslatef(0, width - palla / 2, height / 3);
        drawCube(palla, palla / 2, 2 * height / 3, color);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0, width - palla - mirrorWidth, height / 3);
        drawCube(palla, palla / 2, 2 * height / 3, color);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(len - palla, width - palla / 2, height / 3);
        drawCube(palla, palla / 2, 2 * height / 3, color);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(len - palla, width - palla - mirrorWidth, height / 3);
        drawCube(palla, palla / 2, 2 * height / 3, color);
    glPopMatrix();

    color = {192.0 / 255, 235.0 / 255, 234.0 / 255};

    glPushMatrix();
        glTranslatef(0, width - palla / 2 - mirrorWidth, height / 3);
        drawCube(len, mirrorWidth, 2 * height / 3, color);
    glPopMatrix();
}

void fromThetatoTheta(float r, float from, float to){
    glPushMatrix();
        glBegin(GL_TRIANGLES);
            for(float theta = from; theta < to; theta++){
                getNormal3p(0, 0, 0, r * cos(theta * M_PI / 180), r * sin(theta * M_PI / 180), 0, r * cos((theta + 1) * M_PI / 180), r * sin((theta + 1) * M_PI / 180), 0);
                glVertex3f(0, 0, 0);
                glVertex3f(r * cos(theta * M_PI / 180), r * sin(theta * M_PI / 180), 0);
                glVertex3f(r * cos((theta + 1) * M_PI / 180), r * sin((theta + 1) * M_PI / 180), 0);
            }
        glEnd();
    glPopMatrix();
}

void myLight(float r, float h, vector<float> color, bool glow){
    glPushMatrix();

    glTranslatef(0, 0, -(r * cos(20 * M_PI / 180) + h));

    if(glow) colorMaterialEmit(color);
    else {
        colorMaterial(color);
    }

    glPushMatrix();
        drawSphere(r, 80, 80, color, 1);
    glPopMatrix();

    colorMaterial({0.2, 0.2, 0.2});
    glPushMatrix();
        glTranslatef(0, 0, r * cos(20 * M_PI / 180));
        drawCylinder(r * sin(20 * M_PI / 180), h, 60, 60, {0.2, 0.2, 0.2});
    glPopMatrix();
    glPopMatrix();
}

void drawLamp(float H, float R, float r, float h, float floorH, float floorS, float standR){
    vector<vector<float>> p1, p2;
    if(R < r) swap(R, r);

    for(float theta = 0; theta <= 360; theta += 3){
            p1.push_back({r * cos(theta * M_PI / 180), r * sin(theta * M_PI / 180), h});
            p2.push_back({R * cos(theta * M_PI / 180), R * sin(theta * M_PI / 180), 0});
    }

    float height = H - h + h / 2 - floorH;

    glPushMatrix();
        glTranslatef(0, 0, H - h);
        colorMaterial({0.2, 0.2, 0.2});
        glPushMatrix(); // Uporer jali
            glBegin(GL_LINES);
                for(int i = 0; i < p1.size() - 1; i++){
                    glVertex3f(0, 0, h - h / 20);
                    glVertex3fv(&p1[i][0]);
                }
            glEnd();
        glPopMatrix();

        colorMaterial({0.8, 0, 0});

        glPushMatrix();
            glBegin(GL_QUADS);
                for(int i = 0; i < p1.size() - 1; i++){
                    getNormal3p(p2[i][0], p2[i][1], p2[i][2],
                                p2[i + 1][0], p2[i + 1][1], p2[i + 1][2],
                                p1[i + 1][0], p1[i + 1][1], p1[i + 1][2]);

                    glVertex3fv(&p2[i][0]);
                    glVertex3fv(&p2[i + 1][0]);
                    glVertex3fv(&p1[i + 1][0]);
                    glVertex3fv(&p1[i][0]);
                }
            glEnd();
        glPopMatrix();
    glPopMatrix();


    glPushMatrix();
        glTranslatef(-floorS / 2, -floorS / 2, 0);
        drawCube(floorS,floorS, floorH, {0.2, 0.2, 0.2});
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0, 0, floorH);
        drawCylinder(standR, height, 10, 120, {0.2, 0.2, 0.2});
    glPopMatrix();

    float newR = r + (R - r) / 2;

    colorMaterial({0.2, 0.2, 0.2});
    glPushMatrix();
        glTranslatef(0, 0, H - h / 2);
        fromThetatoTheta(newR, 80, 100);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0, 0, H - h / 2);
        glRotatef(180, 0, 0, 1);
        fromThetatoTheta(newR, 80, 100);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0, 0, H - h / 2);
        glRotatef(90, 0, 0, 1);
        fromThetatoTheta(newR, 80, 100);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0, 0, H - h / 2);
        glRotatef(-90, 0, 0, 1);
        fromThetatoTheta(newR, 80, 100);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0, 0, H - h / 2);
        glRotatef(-180, 0, 1, 0);
        myLight(h / 8, h / 16, {Yellow}, lamp);
    glPopMatrix();
}

void light(){
    GLfloat no_light[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat light_ambient[]  = { 0.2, 0.2, 0.2, 1.0 };
    GLfloat light_diffuse[]  = { 1, 1, 1, 1.0 };
    GLfloat light_specular[] = { 1, 1, 1, 1.0 };
    GLfloat light_position[] = {3 + 4 + 15 + (27 - 15 / 2), 7 + 7, 8 + 3, 1.0 };

    glEnable( GL_LIGHT0);

    if(Mainlight){
        glLightfv( GL_LIGHT0, GL_AMBIENT, light_ambient);
        glLightfv( GL_LIGHT0, GL_DIFFUSE, light_diffuse);
        glLightfv( GL_LIGHT0, GL_SPECULAR, light_specular);
        glLightfv( GL_LIGHT0, GL_POSITION, light_position);
    }
    else{
        glLightfv( GL_LIGHT0, GL_AMBIENT, no_light);
        glLightfv( GL_LIGHT0, GL_DIFFUSE, no_light);
        glLightfv( GL_LIGHT0, GL_SPECULAR, no_light);
        glLightfv( GL_LIGHT0, GL_POSITION, light_position);
    }

    GLfloat light1_ambient[]  = { 0, 0, 0.2, 1.0 };
    GLfloat light1_diffuse[]  = { 0, 0, 0.6, 1.0 };
    GLfloat light1_specular[] = { 0, 0, 0, 1.0 };
    GLfloat light1_position[] = { 25, 25, 40, 1.0 };

    glEnable( GL_LIGHT1);

    glLightfv( GL_LIGHT1, GL_AMBIENT, light1_ambient);
    glLightfv( GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
    glLightfv( GL_LIGHT1, GL_SPECULAR, light1_specular);
    glLightfv( GL_LIGHT1, GL_POSITION, light1_position);

    GLfloat light2_ambient[]  = { Yellow, 1.0 };
    GLfloat light2_diffuse[]  = { Yellow, 1.0 };
    GLfloat light2_specular[] = { Yellow, 1.0 };
    GLfloat light2_position[] = {8.5, 15.5, 10, 1.0}; //Lamp

    glEnable( GL_LIGHT2);

    if(lamp){
        glLightfv( GL_LIGHT2, GL_AMBIENT, light2_ambient);
        glLightfv( GL_LIGHT2, GL_DIFFUSE, light2_diffuse);
        glLightfv( GL_LIGHT2, GL_SPECULAR, light2_specular);
        glLightfv( GL_LIGHT2, GL_POSITION, light2_position);
    }
    else{
        glLightfv( GL_LIGHT2, GL_AMBIENT, no_light);
        glLightfv( GL_LIGHT2, GL_DIFFUSE, no_light);
        glLightfv( GL_LIGHT2, GL_SPECULAR, no_light);
    }

    GLfloat spot_direction[] = { 0, 0.0, -1.0};
    glLightfv(GL_LIGHT2, GL_SPOT_DIRECTION, spot_direction);
    glLightf( GL_LIGHT2, GL_SPOT_CUTOFF, 30.0);
}

vector<vector<long long>> nCr;

void calnCr(int n, int r){
    for(int i = 1; i < n; i++){
        long long curr = i;
        nCr[i][0] = 1;
        for(int j = 1; j <= i; j++){
            nCr[i][j] = curr;
            curr *= i - j;
            curr /= (j + 1);
        }
    }
}

void BezierCurve ( double t,  float xy[2], vector<vector<GLfloat>> ctrlpoints){
    double y=0;
    double x=0;
    t = t > 1.0 ? 1.0 : t;

    long long L = ctrlpoints.size() - 1;

    for(int i=0; i<=L; i++)
    {
        int ncr = nCr[L][i];
        double oneMinusTpow=pow(1-t,double(L-i));
        double tPow=pow(t,double(i));
        double coef=oneMinusTpow*tPow*ncr;
        x+=coef*ctrlpoints[i][0];
        y+=coef*ctrlpoints[i][1];

    }
    xy[0] = float(x);
    xy[1] = float(y);

    //return y;
}

void setNormal(GLfloat x1, GLfloat y1,GLfloat z1, GLfloat x2, GLfloat y2,GLfloat z2, GLfloat x3, GLfloat y3,GLfloat z3){
    GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

    Ux = x2-x1;
    Uy = y2-y1;
    Uz = z2-z1;

    Vx = x3-x1;
    Vy = y3-y1;
    Vz = z3-z1;

    Nx = Uy*Vz - Uz*Vy;
    Ny = Uz*Vx - Ux*Vz;
    Nz = Ux*Vy - Uy*Vx;

    glNormal3f(-Nx,-Ny,-Nz);
}

void draw_curve(vector<vector<GLfloat>> ctrlpoints){
    int i, j;
    float x, y, z, r;				//current coordinates
    float x1, y1, z1, r1;			//next coordinates
    float theta;
    int L = ctrlpoints.size() - 1;

    const float startx = 0, endx = ctrlpoints[L][0];
    //number of angular slices
    const float dx = (endx - startx) / nt;	//x step size
    const float dtheta = 2*M_PI / ntheta;		//angular step size

    float t=0;
    float dt=1.0/nt;
    float xy[2];

    BezierCurve( t,  xy, ctrlpoints);
    x = xy[0];
    r = xy[1];
    //rotate about z-axis
    float p1x,p1y,p1z,p2x,p2y,p2z;
    for ( i = 0; i < nt; ++i )  			//step through x
    {
        theta = 0;
        t+=dt;
        BezierCurve( t,  xy, ctrlpoints);
        x1 = xy[0];
        r1 = xy[1];

        //draw the surface composed of quadrilaterals by sweeping theta
        glBegin( GL_QUAD_STRIP );
        for ( j = 0; j <= ntheta; ++j )
        {
            theta += dtheta;
            double cosa = cos( theta );
            double sina = sin ( theta );
            y = r * cosa;
            y1 = r1 * cosa;	//current and next y
            z = r * sina;
            z1 = r1 * sina;	//current and next z

            //edge from point at x to point at next x
            glVertex3f (x, y, z);

            if(j>0) setNormal(p1x,p1y,p1z,p2x,p2y,p2z,x, y, z);
            else p1x = x, p1y = y, p1z = z, p2x = x1, p2y = y1, p2z = z1;

            glVertex3f (x1, y1, z1);
            //forms quad with next pair of points with incremented theta value
        }
        glEnd();
        x = x1;
        r = r1;
    } //for i
}

void draw_baluster_up(){
    vector<vector<GLfloat>> ctrlpoints = {
        { 0.0, 0.0, 0.0}, { 0, 0.15, 0.0},
        { 0.3, 0.15, 0.0},{ 0.5, 0.15, 0.0},
        { 0.7, 0.1, 0}, { 1, 0.07, 0}
    };

    glPushMatrix();
        glRotatef(-90, 0, 1, 0);
        draw_curve(ctrlpoints);
    glPopMatrix();
}

void draw_baluster_down(){
    vector<vector<GLfloat>> ctrlpoints = {
        {0,0,0},{ 0, 0.1, 0.0},
        { 0, 0.2, 0}, { 0, 0.3, 0},
        {0.3, 0.3, 0}, {0.4, 0.3, 0},
        {0.8, 0.2, 0},
        {1.2, 0.07, 0}, {1.3, 0.07, 0},
        {1.8, 0.07, 0}, {1.9, 0.09, 0}, {2, 0.16, 0},
        {2.1, 0.22, 0},
        {2.1, 0.2, 0}, {2.1, 0.18, 0},
        {2.1, 0.1, 0}, {2.08, 0, 0}
    };

    glPushMatrix();
        glRotatef(-90, 0, 1, 0);
        draw_curve(ctrlpoints);
    glPopMatrix();
}

void draw_baluster(vector<GLfloat> color){
    colorMaterial(color);

    glPushMatrix();
        draw_baluster_down();
        glTranslatef(0, 0, 2.08);
        draw_baluster_up();
    glPopMatrix();
}

void draw_jointer(float r, vector<GLfloat> color){
    vector<point> p1;
    colorMaterial(color);

    for(float theta = 0; theta <= 180; theta+=1){
        p1.push_back(point(r * cos(theta * M_PI / 180), r * cos(theta * M_PI / 180), r * sin(theta * M_PI / 180)));
    }

    point tmp = p1[0];
    tmp.z -= 0.2; // extra height
    p1.insert(p1.begin(), tmp);
    tmp = p1.back();
    tmp.z -= 0.2;
    p1.push_back(tmp);

//    glPushMatrix();
//    glBegin(GL_LINES);
//    for(int i = 0; i < p1.size() - 1; i++){
//        glVertex3f(p1[i].x, p1[i].y, p1[i].z);
//        glVertex3f(p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);
//    }
//    glEnd();
//    glPopMatrix();

    vector<point> p2;
    for(float theta = 0; theta <= 180; theta+=1){
        p2.push_back(point(r * cos(theta * M_PI / 180), -r, r * sin(theta * M_PI / 180)));
    }

    tmp = p2[0];
    tmp.z -= 0.2; // extra height
    p2.insert(p2.begin(), tmp);
    tmp = p2.back();
    tmp.z -= 0.2;
    p2.push_back(tmp);

    vector<point> p3;
    for(float theta = 0; theta <= 180; theta+=1){
        p3.push_back(point(-r, r * cos(theta * M_PI / 180), r * sin(theta * M_PI / 180)));
    }

    tmp = p3[0];
    tmp.z -= 0.2; // extra height
    p3.insert(p3.begin(), tmp);
    tmp = p3.back();
    tmp.z -= 0.2;
    p3.push_back(tmp);

    glPushMatrix();
        glBegin(GL_QUADS);
            for(int i = 0; i < p1.size() - 1; i++){
                getNormal3p(p2[i].x, p2[i].y, p2[i].z,
                            p1[i].x, p1[i].y, p1[i].z,
                            p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);

                glVertex3f(p2[i].x, p2[i].y, p2[i].z);
                glVertex3f(p1[i].x, p1[i].y, p1[i].z);
                glVertex3f(p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);
                glVertex3f(p2[i + 1].x, p2[i + 1].y, p2[i + 1].z);
            }

            for(int i = 0; i < p1.size() - 1; i++){
                getNormal3p(p3[i].x, p3[i].y, p3[i].z,
                            p3[i + 1].x, p3[i + 1].y, p3[i + 1].z,
                            p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);

                glVertex3f(p3[i].x, p3[i].y, p3[i].z);
                glVertex3f(p3[i + 1].x, p3[i + 1].y, p3[i + 1].z);
                glVertex3f(p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);
                glVertex3f(p1[i].x, p1[i].y, p1[i].z);
            }
        glEnd();
    glPopMatrix();

    glPushMatrix();
        glBegin(GL_TRIANGLES);
            getNormal3p(p2[0].x, p2[0].y, p2[0].z, p2.back().x, p2.back().y, p2.back().z, p1[0].x, p1[0].y, p1[0].z);
            glVertex3f(p2[0].x, p2[0].y, p2[0].z);
            glVertex3f(p2.back().x, p2.back().y, p2.back().z);
            glVertex3f(p1[0].x, p1[0].y, p1[0].z);

            getNormal3p(p3.back().x, p3.back().y, p3.back().z, p3[0].x, p3[0].y, p3[0].z, p1[0].x, p1[0].y, p1[0].z);
            glVertex3f(p3.back().x, p3.back().y, p3.back().z);
            glVertex3f(p3[0].x, p3[0].y, p3[0].z);
            glVertex3f(p1[0].x, p1[0].y, p1[0].z);
        glEnd();
    glPopMatrix();
}

void draw_stair_head(float theta, float len, float r, float ex, vector<GLfloat> color){
    vector<point> p1, p2;
    colorMaterial(color);

    for(GLfloat ctheta = 0; ctheta <= 180; ctheta += 1){
        p1.push_back(point(0, r * cos(ctheta * M_PI / 180), r * sin(ctheta * M_PI / 180)));
    }

    for(GLfloat ctheta = 0; ctheta <= 180; ctheta += 1){
        p2.push_back(point(len * cos(theta * M_PI / 180), r * cos(ctheta * M_PI / 180), r * sin(ctheta * M_PI / 180) - len * sin(theta * M_PI / 180)));
    }

    point c1 = p1.back();
    point c2 = p2.back();
    c2.z -= 0.2;
    c1.z -= 0.2; // extra height
    p1.push_back(c1);
    p2.push_back(c2);

    c1.y += 2 *  r;
    c2.y += 2 * r;
    p1.insert(p1.begin(), c1);
    p2.insert(p2.begin(), c2);

    glBegin(GL_QUADS);
        for(int i = 0; i < p1.size() - 1; i++){
            getNormal3p(p2[i].x, p2[i].y, p2[i].z,
                        p1[i].x, p1[i].y, p1[i].z,
                        p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);

            glVertex3f(p2[i].x, p2[i].y, p2[i].z);
            glVertex3f(p1[i].x, p1[i].y, p1[i].z);
            glVertex3f(p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);
            glVertex3f(p2[i + 1].x, p2[i + 1].y, p2[i + 1].z);
        }

        for(int i = p2.size() - 1, j = 0; i - 1 >= j + 1; i--, j++){
            getNormal3p(p2[i].x, p2[i].y, p2[i].z, p2[j].x, p2[j].y, p2[j].z, p2[j + 1].x, p2[j + 1].y, p2[j + 1].z);

            glVertex3f(p2[i].x, p2[i].y, p2[i].z);
            glVertex3f(p2[j].x, p2[j].y, p2[j].z);
            glVertex3f(p2[j + 1].x, p2[j + 1].y, p2[j + 1].z);
            glVertex3f(p2[i - 1].x, p2[i - 1].y, p2[i - 1].z);
        }

        for(int i = p1.size() - 1, j = 0; i - 1 >= j + 1; i--, j++){
            getNormal3p(p1[j].x, p1[j].y, p1[j].z, p1[i].x, p1[i].y, p1[i].z, p1[i - 1].x, p1[i - 1].y, p1[i - 1].z);

            glVertex3f(p1[j].x, p1[j].y, p1[j].z);
            glVertex3f(p1[i].x, p1[i].y, p1[i].z);
            glVertex3f(p1[i - 1].x, p1[i - 1].y, p1[i - 1].z);
            glVertex3f(p1[j + 1].x, p1[j + 1].y, p1[j + 1].z);
        }

        getNormal3p(p1.back().x, p1.back().y, p1.back().z, p1[0].x, p1[0].y, p1[0].z, p2[0].x, p2[0].y, p2[0].z);
        glVertex3f(p1.back().x, p1.back().y, p1.back().z);
        glVertex3f(p1[0].x, p1[0].y, p1[0].z);
        glVertex3f(p2[0].x, p2[0].y, p2[0].z);
        glVertex3f(p2.back().x, p2.back().y, p2.back().z);

        for(int i = 0; i < p1.size() - 1; i++){
            getNormal3p(p1[i].x - (ex - r), p1[i].y, p1[i].z, p1[i + 1].x - (ex - r), p1[i + 1].y, p1[i + 1].z, p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);

            glVertex3f(p1[i].x - (ex - r), p1[i].y, p1[i].z);
            glVertex3f(p1[i + 1].x - (ex - r), p1[i + 1].y, p1[i + 1].z);
            glVertex3f(p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);
            glVertex3f(p1[i].x, p1[i].y, p1[i].z);
        }

        getNormal3p(p1.back().x - (ex - r), p1.back().y, p1.back().z, p1[0].x - (ex - r), p1[0].y, p1[0].z, p1[0].x, p1[0].y, p1[0].z);
        glVertex3f(p1.back().x - (ex - r), p1.back().y, p1.back().z);
        glVertex3f(p1[0].x - (ex - r), p1[0].y, p1[0].z);
        glVertex3f(p1[0].x, p1[0].y, p1[0].z);
        glVertex3f(p1.back().x, p1.back().y, p1.back().z);
    glEnd();
}

void draw_stair(GLfloat x, GLfloat y, GLfloat z, GLfloat ex, int n){
    GLfloat perX = (x - ex) / n;
    GLfloat perZ = z / n;

    colorMaterial(getColor(255, 255, 255));

    glPushMatrix();
        glEnable(GL_TEXTURE_2D);
        glPushMatrix();
        glBindTexture(GL_TEXTURE_2D, ID[0]);
        glBegin(GL_QUADS);
            for(int i = 0; i < n; i++){ // top Floor
                getNormal3p((i != 0) * ex + i * perX, 0, z - i * perZ,
                            ex + (i + 1) * perX + 0.05, 0, z - i * perZ,
                            ex + (i + 1) * perX + 0.05, y, z - i * perZ);

                glVertex3f((i != 0) * ex + i * perX, 0, z - i * perZ);                glTexCoord2f(0, 1);
                glVertex3f(ex + (i + 1) * perX + 0.05, 0, z - i * perZ);               glTexCoord2f(1, 1);
                glVertex3f(ex + (i + 1) * perX + 0.05, y, z - i * perZ);               glTexCoord2f(1, 0);
                glVertex3f((i != 0) * ex + i * perX, y, z - i * perZ);                glTexCoord2f(0, 0);
            }

            for(int i = 0; i < n; i++){ //Side front
                getNormal3p((i != 0) * ex + i * perX, 0, 0,
                            ex + (i + 1) * perX, 0, 0,
                            ex + (i + 1) * perX, 0, z - (i * perZ));

                glVertex3f((i != 0) * ex + i * perX, 0, 0);                 glTexCoord2f(0, 1);
                glVertex3f(ex + (i + 1) * perX, 0, 0);                      glTexCoord2f(1, 1);
                glVertex3f(ex + (i + 1) * perX, 0, z - (i * perZ));         glTexCoord2f(1, 0);
                glVertex3f((i != 0) * ex + i * perX, 0, z - (i * perZ));    glTexCoord2f(0, 0);
            }

            for(int i = 0; i < n; i++){ //Side back
                getNormal3p((i != 0) * ex + i * perX, y, z - (i * perZ),
                            ex + (i + 1) * perX, y, z - (i * perZ),
                            ex + (i + 1) * perX, y, 0);

                glVertex3f((i != 0) * ex + i * perX, y, z - (i * perZ));    glTexCoord2f(0, 1);
                glVertex3f(ex + (i + 1) * perX, y, z - (i * perZ));         glTexCoord2f(1, 1);
                glVertex3f(ex + (i + 1) * perX, y, 0);                      glTexCoord2f(1, 0);
                glVertex3f((i != 0) * ex + i * perX, y, 0);                 glTexCoord2f(0, 0);
            }
        glEnd();
        glPopMatrix();
        glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    colorMaterial(getColor(8, 4, 0));

    glPushMatrix();
        glBegin(GL_QUADS);
        for(int i = 0; i < n; i++){
            getNormal3p(ex + (i + 1) * perX, 0, z - (i + 1) * perZ,
                        ex + (i + 1) * perX, y, z - (i + 1) * perZ,
                        ex  + (i + 1) * perX, y, z - i * perZ);

            glVertex3f(ex + (i + 1) * perX, 0, z - (i + 1) * perZ);
            glVertex3f(ex + (i + 1) * perX, y, z - (i + 1) * perZ);
            glVertex3f(ex  + (i + 1) * perX, y, z - i * perZ);
            glVertex3f(ex  + (i + 1) * perX, 0, z - i * perZ);
        }
        glEnd();
    glPopMatrix();

    for(int i = 0; i < n; i++){
        glPushMatrix();
        glTranslatef(ex + i * perX + perX / 2, 0.5, z - i * perZ);
        glScalef(0.9, 0.9, 0.9);
        draw_baluster(getColor(18, 14, 3));
        glPopMatrix();
    }

    for(int i = 0; i < n; i++){
        glPushMatrix();
        glTranslatef(ex + i * perX + perX / 2, y - 0.5, z - i * perZ);
        glScalef(0.9, 0.9, 0.9);
        draw_baluster(getColor(18, 14, 3));
        glPopMatrix();
    }

    point p1, p2;
    p1.x = ex + perX / 2 - 0.07;
    p1.y = 0.5;
    p1.z = 6 + 0.085;

    p2.z = perZ + 3 + 0.085;
    p2.y = 0.5;
    p2.x = ex + (n - 1) * perX - 0.07 + perX / 2;

    float thetaRad = atan(abs((p2.z - p1.z) / (p2.x - p1.x)));

    float theta = (180.0 * thetaRad / M_PI);
    float len = sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.z - p2.z) * (p1.z - p2.z)) + perX;


    glPushMatrix();
        glTranslatef(ex, 0.5, 6 + 0.085 + 0.2);
        draw_stair_head(theta, len, 0.3, ex + 0.5, getColor(46, 35, 3));
    glPopMatrix();

    glPushMatrix();
        glTranslatef(ex, y - 0.5, 6 + 0.085 + 0.2);
        draw_stair_head(theta, len, 0.3, ex + 0.5, getColor(46, 35, 3));
    glPopMatrix();

    glPushMatrix();
        glTranslatef(-0.5, 0.5, 6 + 0.085 + 0.2);
        glRotatef(90, 0, 0, 1);
        draw_jointer(0.3, getColor(46, 35, 3));
    glPopMatrix();

    glPushMatrix();
        glTranslatef(-0.5, y - 0.5, 6 + 0.085 + 0.2);
        glRotatef(90, 0, 0, 1);
        glScalef(-1, 1, 1);
        draw_jointer(0.3, getColor(46, 35, 3));
    glPopMatrix();
}

void draw_ralling(float len, float r, vector<GLfloat> color){
    colorMaterial(color);
    vector<point> p1;
    p1.push_back(point(r, 0, -0.2));
    for(float theta = 0; theta <= 180; theta++){
        p1.push_back(point(r * cos(theta * M_PI / 180), 0, r * sin(theta * M_PI / 180)));
    }
    p1.push_back(point(p1.back().x, p1.back().y, p1.back().z - 0.2));

    vector<point> p2;
    p2.push_back(point(r, len, -0.2));
    for(float theta = 0; theta <= 180; theta++){
        p2.push_back(point(r * cos(theta * M_PI / 180), len, r * sin(theta * M_PI / 180)));
    }
    p2.push_back(point(p2.back().x, p2.back().y, p2.back().z - 0.2));

    glPushMatrix();
        glBegin(GL_QUADS);
            for(int i = 0; i < p1.size() - 1; i++){
                getNormal3p(p1[i].x, p1[i].y, p1[i].z, p2[i].x, p2[i].y, p2[i].z, p2[i + 1].x, p2[i + 1].y, p2[i + 1].z);

                glVertex3f(p1[i].x, p1[i].y, p1[i].z);
                glVertex3f(p2[i].x, p2[i].y, p2[i].z);
                glVertex3f(p2[i + 1].x, p2[i + 1].y, p2[i + 1].z);
                glVertex3f(p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);
            }

            for(int i = p1.size() - 1, j = 0; i - 1 >= j + 1; i--, j++){
                getNormal3p(p1[i].x, p1[i].y, p1[i].z, p1[j].x, p1[j].y, p1[j].z, p1[j + 1].x, p1[j + 1].y, p1[j + 1].z);

                glVertex3f(p1[i].x, p1[i].y, p1[i].z);
                glVertex3f(p1[j].x, p1[j].y, p1[j].z);
                glVertex3f(p1[j + 1].x, p1[j + 1].y, p1[j + 1].z);
                glVertex3f(p1[i - 1].x, p1[i - 1].y, p1[i - 1].z);
            }

            for(int i = p2.size() - 1, j = 0; i - 1 >= j + 1; i--, j++){
                getNormal3p(p2[i].x, p2[i].y, p2[i].z, p2[i - 1].x, p2[i - 1].y, p2[i - 1].z, p2[j + 1].x, p2[j + 1].y, p2[j + 1].z);

                glVertex3f(p2[i].x, p2[i].y, p2[i].z);
                glVertex3f(p2[i - 1].x, p2[i - 1].y, p2[i - 1].z);
                glVertex3f(p2[j + 1].x, p2[j + 1].y, p2[j + 1].z);
                glVertex3f(p2[j].x, p2[j].y, p2[j].z);
            }

            getNormal3p(p1.back().x, p1.back().y, p1.back().z, p2.back().x, p2.back().y, p2.back().z, p2[0].x, p2[0].y, p2[0].z);

            glVertex3f(p1.back().x, p1.back().y, p1.back().z);
            glVertex3f(p2.back().x, p2.back().y, p2.back().z);
            glVertex3f(p2[0].x, p2[0].y, p2[0].z);
            glVertex3f(p1[0].x, p1[0].y, p1[0].z);
        glEnd();
    glPopMatrix();
}

void draw_baluster_series(int n, float len, vector<GLfloat> color){ // draws in X direction
    float perX = len / (n - 1);

    glPushMatrix();
    for(int i = 0; i < n; i++){
        glTranslatef((i > 0) * perX, 0, 0);
        draw_baluster(color);
    }
    glPopMatrix();
}

void draw_cube_textured(float x, float y, float z, vector<float> color, int ti){
    colorMaterial({1, 1, 1});
    vector<vector<GLfloat>> points = generatePoints8p(x, y, z);
    vector<vector<GLubyte>> ind = getMyIndicesForCube();

    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glBindTexture(GL_TEXTURE_2D, ID[ti]);
    glBegin(GL_QUADS);

        for(int i = 0; i < 6; i++){

            getNormal3p(
                        points[ind[i][0]][0], points[ind[i][0]][1], points[ind[i][0]][2],
                        points[ind[i][1]][0], points[ind[i][1]][1], points[ind[i][1]][2],
                        points[ind[i][2]][0], points[ind[i][2]][1], points[ind[i][2]][2]
            );

            for(int j = 0; j < 4; j++){
                glVertex3fv(&points[ind[i][j]][0]);

                if(j == 0) glTexCoord2f(1,0);
                if(j == 1) glTexCoord2f(1,1);
                if(j == 2) glTexCoord2f(0,1);
                if(j == 3) glTexCoord2f(0,0);
            }
        }

    glEnd();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
}

void drawPhoto(float len){
    float extraWidth = 0.1;
    len -= 2 * extraWidth;
    glPushMatrix();
        glTranslatef(len / 2, len / 2, 0);
        glRotatef(90, 0, 0, 1);

        glPushMatrix(); // Frame
            glTranslatef(len / 2, -len / 2, 0);

            glPushMatrix(); // Nicher Palla
                glRotatef(180, 0, 1, 0);
                draw_cube_textured(len, len, extraWidth, {0, 0, 0}, 1);
            glPopMatrix();

            glPushMatrix(); //Border
                //glTranslatef(-len, extraWidth, 0);
                glRotatef(90, 0, 1, 0);
                glRotatef(90, 0, 0, 1);
                drawFlatPyramid(len, 0.1, extraWidth, ColorLiteWood, ColorExtraDarkWood);
            glPopMatrix();

            glPushMatrix(); //Border
                glTranslatef(-len, 0, extraWidth);
                glRotatef(180, 0, 1, 0);
                glRotatef(90, 0, 1, 0);
                glRotatef(90, 0, 0, 1);
                drawFlatPyramid(len, 0.1, extraWidth, ColorLiteWood, ColorExtraDarkWood);
            glPopMatrix();

            glPushMatrix(); //Border
                glTranslatef(0, len, 0);
                glRotatef(90, 0, 0, 1);
                glPushMatrix();
                    //glTranslatef(-len, extraWidth, 0);
                    glRotatef(90, 0, 1, 0);
                    glRotatef(90, 0, 0, 1);
                    drawFlatPyramid(len, 0.1, extraWidth, ColorLiteWood, ColorExtraDarkWood);
                glPopMatrix();

                glPushMatrix();
                    glTranslatef(-len, 0, extraWidth);
                    glRotatef(180, 0, 1, 0);
                    glRotatef(90, 0, 1, 0);
                    glRotatef(90, 0, 0, 1);
                    drawFlatPyramid(len, 0.1, extraWidth, ColorLiteWood, ColorExtraDarkWood);
                glPopMatrix();
            glPopMatrix();
        glPopMatrix();
    glPopMatrix();
}

void draw_basement(){
    /*
    length (x) = 35
         + 5 (stairs):
                x = 5
                y = 15
                z = 3
                posx = 35
                posy = 15
    width (y) = 45
    Height (z) = 3
    */

    draw_cube_textured(35, 45, 3, {1, 1, 1}, 0);

    //Stairs
    glPushMatrix();
        glTranslatef(35, 15, 0);
        draw_stair(5, 15, 3, 1, 4);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(35 - 0.5, 0.5 + 0.3, 6 + 0.085 + 0.2);
        draw_ralling(15 - 2 * 0.3, 0.3, getColor(46, 35, 3));
    glPopMatrix();

    glPushMatrix();
        glTranslatef(35 - 0.5, 0.5, 6 + 0.085 + 0.2);
        glRotatef(-90, 0, 0, 1);
        draw_jointer(0.3, getColor(46, 35, 3));
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0.5 + 0.3, 0.5, 6 + 0.085 + 0.2);
        glRotatef(-90, 0, 0, 1);
        draw_ralling(35 - 2 * 0.3 - 2 * 0.5, 0.3, getColor(46, 35, 3));
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0.5, 0.5, 6 + 0.085 + 0.2);
        glRotatef(90, 0, 0, 1);
        glScalef(-1, 1, 1);
        draw_jointer(0.3, getColor(46, 35, 3));
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0.5, 0.5 + 0.3, 6 + 0.085 + 0.2);
        draw_ralling(45 - 2 * 0.3 - 2 * 0.5, 0.3, getColor(46, 35, 3));
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0.5, 45 - 0.5, 6 + 0.085 + 0.2);
        glRotatef(90, 0, 0, 1);
        draw_jointer(0.3, getColor(46, 35, 3));
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0.5 + 0.3, 45 - 0.5, 6 + 0.085 + 0.2);
        glRotatef(-90, 0, 0, 1);
        draw_ralling(35 - 2 * 0.3 - 2 * 0.5, 0.3, getColor(46, 35, 3));
    glPopMatrix();

    glPushMatrix();
        glTranslatef(35 - 0.5, 45 - 0.5, 6 + 0.085 + 0.2);
        draw_jointer(0.3, getColor(46, 35, 3));
    glPopMatrix();

    glPushMatrix();
        glTranslatef(35 - 0.5, 30 - 0.5 + 0.3, 6 + 0.085 + 0.2);
        draw_ralling(15 - 2 * 0.3, 0.3, getColor(46, 35, 3));
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0.5, 0.5, 3);
        draw_baluster_series(18,34, getColor(18, 14, 3));
    glPopMatrix();

    glPushMatrix();
        glTranslatef(35 - 0.5, 0.5 + (15.0 / 7.0), 3);
        glRotatef(90, 0, 0, 1);
        draw_baluster_series(8,15.0 - (15.0 / 7.0), getColor(18, 14, 3));
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0.5, 45 - 0.5, 3);
        draw_baluster_series(18,34, getColor(18, 14, 3));
    glPopMatrix();

    glPushMatrix();
        glTranslatef(35 - 0.5, 45 - (0.5 + (15.0 / 7.0)) - (15.0 - (15.0 / 7.0)), 3);
        glRotatef(90, 0, 0, 1);
        draw_baluster_series(7,15.0 - (15.0 / 7.0), getColor(18, 14, 3));
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0.5, 0.5 + (44.0 / 19.0), 3);
        glRotatef(90, 0, 0, 1);
        draw_baluster_series(20,44.0 - 2 * (44.0 / 19.0), getColor(18, 14, 3));
    glPopMatrix();
}

void draw_axes(vector<GLfloat> color){
    colorMaterial(color);

    glPushMatrix();
        glBegin(GL_LINES);
            glColor3f(1, 0, 0);
            glVertex3f(-5, 0, 0);
            glVertex3f(5, 0, 0);

            glColor3f(0, 1, 0);
            glVertex3f(0, -5, 0);
            glVertex3f(0, 5, 0);

            glColor3f(0, 0, 1);
            glVertex3f(0, 0, -5);
            glVertex3f(0, 0, 5);
        glEnd();
    glPopMatrix();
}

void draw_area(GLfloat x, GLfloat y, vector<GLfloat> color){
    colorMaterial(color);

    glBegin(GL_QUADS);
        getNormal3p(0, 0, 0, x, 0, 0, x, y, 0);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(x, 0.0, 0.0);
        glVertex3f(x, y, 0.0);
        glVertex3f(0.0, y, 0.0);
    glEnd();
}

void teleJali(vector<point> p1){
    colorMaterial({0.1, 0.1, 0.1});
    glBegin(GL_LINES);
        for(int i = 30 + 1, j = 30 - 1; j >= 0; j--, i++){
            glVertex3f(p1[i].x, p1[i].y, p1[i].z);
            glVertex3f(p1[j].x, p1[j].y, p1[j].z);
        }
        for(int i = 61, j = 119; j > 90; j--, i++){
            glVertex3f(p1[i].x, p1[i].y, p1[i].z);
            glVertex3f(p1[j].x, p1[j].y, p1[j].z);
        }
        for(int i = 59, j = 61; j <= 90; j++, i--){
            glVertex3f(p1[i].x, p1[i].y, p1[i].z);
            glVertex3f(p1[j].x, p1[j].y, p1[j].z);
        }
        for(int i = 91, j = 29; j > 0; i++, j--){
            glVertex3f(p1[i].x, p1[i].y, p1[i].z);
            glVertex3f(p1[j].x, p1[j].y, p1[j].z);
        }
    glEnd();
}

void teleChaka(float r, float d, float h, vector<point> p1){
    vector<point> p2;
    colorMaterial({.2, .2, .2});
    for(float theta = 0; theta < 360; theta += 3){
        p2.push_back(point((r - d) * cos(theta * M_PI / 180), (r - d) * sin(theta * M_PI / 180), 0));
    }

    glPushMatrix();
        glBegin(GL_QUADS);
            for(int i = 0, n = p1.size(); i < p1.size(); i++){
                getNormal3p(p1[i].x, p1[i].y, p1[i].z,
                            p1[i].x, p1[i].y, p1[i].z - h,
                            p1[(i + 1) % n].x, p1[(i + 1) % n].y, p1[(i + 1) % n].z - h);

                glVertex3f(p1[i].x, p1[i].y, p1[i].z);
                glVertex3f(p1[i].x, p1[i].y, p1[i].z - h);
                glVertex3f(p1[(i + 1) % n].x, p1[(i + 1) % n].y, p1[(i + 1) % n].z - h);
                glVertex3f(p1[(i + 1) % n].x, p1[(i + 1) % n].y, p1[(i + 1) % n].z);
            }

            for(int i = 0, n = p1.size(); i < p2.size(); i++){
                getNormal3p(p2[(i + 1) % n].x, p2[(i + 1) % n].y, p2[(i + 1) % n].z,
                            p2[(i + 1) % n].x, p2[(i + 1) % n].y, p2[(i + 1) % n].z - h,
                            p2[i].x, p2[i].y, p2[i].z - h);

                glVertex3f(p2[(i + 1) % n].x, p2[(i + 1) % n].y, p2[(i + 1) % n].z);
                glVertex3f(p2[(i + 1) % n].x, p2[(i + 1) % n].y, p2[(i + 1) % n].z - h);
                glVertex3f(p2[i].x, p2[i].y, p2[i].z - h);
                glVertex3f(p2[i].x, p2[i].y, p2[i].z);
            }

            for(int i = 0, n = p1.size(); i < p1.size(); i++){
                getNormal3p(p2[i].x, p2[i].y, p2[i].z, p1[i].x, p1[i].y, p1[i].z, p1[(i + 1) % n].x, p1[(i + 1) % n].y, p1[(i + 1) % n].z);

                glVertex3f(p2[i].x, p2[i].y, p2[i].z);
                glVertex3f(p1[i].x, p1[i].y, p1[i].z);
                glVertex3f(p1[(i + 1) % n].x, p1[(i + 1) % n].y, p1[(i + 1) % n].z);
                glVertex3f(p2[(i + 1) % n].x, p2[(i + 1) % n].y, p2[(i + 1) % n].z);
            }

            for(int i = 0, n = p1.size(); i < p1.size(); i++){
                getNormal3p(p2[i].x, p2[i].y, p2[i].z - h, p1[i].x, p1[i].y, p1[i].z - h, p1[(i + 1) % n].x, p1[(i + 1) % n].y, p1[(i + 1) % n].z - h);

                glVertex3f(p2[i].x, p2[i].y, p2[i].z - h);
                glVertex3f(p1[i].x, p1[i].y, p1[i].z - h);
                glVertex3f(p1[(i + 1) % n].x, p1[(i + 1) % n].y, p1[(i + 1) % n].z - h);
                glVertex3f(p2[(i + 1) % n].x, p2[(i + 1) % n].y, p2[(i + 1) % n].z - h);
            }
        glEnd();
    glPopMatrix();
}

void teleWire(float r, float h, float thOwn, float thPos){
    float d = r * sin(thOwn * M_PI / 180);
    colorMaterial({0.2, 0.2, 0.2});

    vector<point> p1;
    for(float theta = 0; theta < 360; theta += 3){
        p1.push_back(point(d * cos(theta * M_PI / 180), d * sin(theta * M_PI / 180), 0));
    }

    glPushMatrix();
        glTranslatef(-sqrt(r * r - d * d) * cos((thPos + thOwn) * M_PI / 180), 0, sqrt(r * r - d * d) * sin((thPos + thOwn) * M_PI / 180));
        glRotatef(-(90 - (thPos + thOwn)), 0, 1, 0);

        glBegin(GL_QUADS);
            for(int i = 0, n = p1.size(); i < p1.size(); i++){
                getNormal3p(p1[i].x, p1[i].y, p1[i].z, p1[(i + 1) % n].x, p1[(i + 1) % n].y, p1[(i + 1) % n].z, p1[(i + 1) % n].x, p1[(i + 1) % n].y, p1[(i + 1) % n].z + h);

                glVertex3f(p1[i].x, p1[i].y, p1[i].z);
                glVertex3f(p1[(i + 1) % n].x, p1[(i + 1) % n].y, p1[(i + 1) % n].z);
                glVertex3f(p1[(i + 1) % n].x, p1[(i + 1) % n].y, p1[(i + 1) % n].z + h);
                glVertex3f(p1[i].x, p1[i].y, p1[i].z + h);
            }

            for(int i = 1, j = p1.size() - 1; i + 1 < 60; i++, j--){
                getNormal3p(p1[i].x, p1[i].y, p1[i].z + h, p1[i + 1].x, p1[i + 1].y, p1[i + 1].z + h, p1[j - 1].x, p1[j - 1].y, p1[j - 1].z + h);

                glVertex3f(p1[i].x, p1[i].y, p1[i].z + h);
                glVertex3f(p1[i + 1].x, p1[i + 1].y, p1[i + 1].z + h);
                glVertex3f(p1[j - 1].x, p1[j - 1].y, p1[j - 1].z + h);
                glVertex3f(p1[j].x, p1[j].y, p1[j].z + h);
            }
        glEnd();

        glBegin(GL_TRIANGLES);
            getNormal3p(p1[0].x, p1[0].y, p1[0].z + h, p1[1].x, p1[1].y, p1[1].z + h, p1.back().x, p1.back().y, p1.back().z + h);

            glVertex3f(p1[0].x, p1[0].y, p1[0].z + h);
            glVertex3f(p1[1].x, p1[1].y, p1[1].z + h);
            glVertex3f(p1.back().x, p1.back().y, p1.back().z + h);

            getNormal3p(p1[61].x, p1[61].y, p1[61].z + h, p1[59].x, p1[59].y, p1[59].z + h, p1[60].x, p1[60].y, p1[60].z + h);
            glVertex3f(p1[61].x, p1[61].y, p1[61].z + h);
            glVertex3f(p1[59].x, p1[59].y, p1[59].z + h);
            glVertex3f(p1[60].x, p1[60].y, p1[60].z + h);
        glEnd();
    glPopMatrix();
}

void draw_telephone(float len, float r){
    vector<point> p1;
    for(float theta = 60; theta <= 120; theta++){
        p1.push_back(point(len * cos(theta * M_PI / 180), 0, len * sin(theta * M_PI / 180)));
    }

    colorMaterial({.2, .2, .2});

    float totalHeight = r + r / 4 + len - len * sin(M_PI / 3);

    glPushMatrix();
    glRotatef(-(asin((r - r / 4) / totalHeight) * 180 / M_PI), 1, 0, 0);
    glTranslatef(0, -(len - len * sin(M_PI / 3)), 0);
    glRotatef(-90, 1, 0, 0);
    glTranslatef(0, -r / 4, 0); //radius of sphere + height of handle helper + height of chaka

    glPushMatrix();
        glTranslatef(0, -r / 4, -len * cos(M_PI / 6));
        glBegin(GL_QUADS);
            for(int i = p1.size() - 1, j = 0; i - 1 >= j + 1; i--, j++){
                getNormal3p(p1[i].x, p1[i].y, p1[i].z, p1[j].x, p1[j].y, p1[j].z, p1[j + 1].x, p1[j + 1].y, p1[j + 1].z);

                glVertex3f(p1[i].x, p1[i].y, p1[i].z);
                glVertex3f(p1[j].x, p1[j].y, p1[j].z);
                glVertex3f(p1[j + 1].x, p1[j + 1].y, p1[j + 1].z);
                glVertex3f(p1[i - 1].x, p1[i - 1].y, p1[i - 1].z);
            }

            for(int i = p1.size() - 1, j = 0; i - 1 >= j + 1; i--, j++){
                getNormal3p(p1[i - 1].x, p1[i - 1].y + r / 2, p1[i - 1].z, p1[j + 1].x, p1[j + 1].y + r / 2, p1[j + 1].z, p1[j].x, p1[j].y + r / 2, p1[j].z);

                glVertex3f(p1[i - 1].x, p1[i - 1].y + r / 2, p1[i - 1].z);
                glVertex3f(p1[j + 1].x, p1[j + 1].y + r / 2, p1[j + 1].z);
                glVertex3f(p1[j].x, p1[j].y + r / 2, p1[j].z);
                glVertex3f(p1[i].x, p1[i].y + r / 2, p1[i].z);
            }

            for(int i = 0; i < p1.size() - 1; i++){
                getNormal3p(p1[i].x, p1[i].y, p1[i].z, p1[i].x, p1[i].y + r / 2, p1[i].z, p1[i + 1].x, p1[i + 1].y + r / 2, p1[i + 1].z);

                glVertex3f(p1[i].x, p1[i].y, p1[i].z);
                glVertex3f(p1[i].x, p1[i].y + r / 2, p1[i].z);
                glVertex3f(p1[i + 1].x, p1[i + 1].y + r / 2, p1[i + 1].z);
                glVertex3f(p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);
            }

            getNormal3p(p1[0].x, p1[0].y, p1[0].z, p1.back().x, p1.back().y, p1.back().z, p1.back().x, p1.back().y + r / 2, p1.back().z);
            glVertex3f(p1[0].x, p1[0].y, p1[0].z);
            glVertex3f(p1.back().x, p1.back().y, p1.back().z);
            glVertex3f(p1.back().x, p1.back().y + r / 2, p1.back().z);
            glVertex3f(p1[0].x, p1[0].y + r / 2, p1[0].z);
        glEnd();
    glPopMatrix();

    p1.clear();

    for(float theta = 0; theta < 360; theta += 3){
        p1.push_back(point(r * cos(theta * M_PI / 180), r * sin(theta * M_PI / 180), 0));
    }

    glPushMatrix();
        glTranslatef(len / 2 - r / 4, 0, -r - r / 4);
        teleWire(r, r / 4, 7, 90 - 7); // Handle Helper
        drawHalfSphere(r, 60, 60, getColor(163, 143, 15));
        teleJali(p1);
        teleChaka(r, 0.05, 0.05, p1);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(-len / 2 + r / 4, 0, -r - r / 4);
        teleWire(r, r, 3, 25); // Chira wire
        teleWire(r, r / 4, 7, 90 - 7); // Handle Helper
        drawHalfSphere(r, 60, 60, getColor(163, 143, 15));
        teleJali(p1);
        teleChaka(r, 0.05, 0.05, p1);
    glPopMatrix();

    glPopMatrix();
}

void squareWithCircleGapHelper(float r){
    vector<point> p1;

    for(float theta = 0; theta <= 90; theta += 3){
        p1.push_back(point(r * cos(theta * M_PI / 180), 0, r * sin(theta * M_PI / 180)));
    }

    point head(r, 0, r);

    glPushMatrix();
        glBegin(GL_TRIANGLES);
            for(int i = 0; i < p1.size() - 1; i++){
                getNormal3p(p1[i].x, p1[i].y, p1[i].z, head.x, head.y, head.z, p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);

                glVertex3f(p1[i].x, p1[i].y, p1[i].z);
                glVertex3f(head.x, head.y, head.z);
                glVertex3f(p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);
            }
        glEnd();
    glPopMatrix();
}

void square_with_circle_gap(float r){
    glPushMatrix();
        squareWithCircleGapHelper(r);
        glRotatef(-90, 0, 1, 0);
        squareWithCircleGapHelper(r);
        glRotatef(-90, 0, 1, 0);
        squareWithCircleGapHelper(r);
        glRotatef(-90, 0, 1, 0);
        squareWithCircleGapHelper(r);
    glPopMatrix();
}

void drawQuadXZ(float x, float z){
    glPushMatrix();
        glBegin(GL_QUADS);
            getNormal3p(0, 0, 0, x, 0, 0, x, 0, z);
            glVertex3f(0, 0, 0);
            glVertex3f(x, 0, 0);
            glVertex3f(x, 0, z);
            glVertex3f(0, 0, z);
        glEnd();
    glPopMatrix();
}

void drawQuadXZBack(float x, float z){
    glPushMatrix();
        glBegin(GL_QUADS);
            getNormal3p(0, 0, z, x, 0, z, x, 0, 0);
            glVertex3f(0, 0, z);
            glVertex3f(x, 0, z);
            glVertex3f(x, 0, 0);
            glVertex3f(0, 0, 0);
        glEnd();
    glPopMatrix();
}

void drawQuadYZ(float y, float z){
    glPushMatrix();
        glBegin(GL_QUADS);
            getNormal3p(0, 0, 0, 0, -y, 0, 0, -y, z);
            glVertex3f(0, 0, 0);
            glVertex3f(0, -y, 0);
            glVertex3f(0, -y, z);
            glVertex3f(0, 0, z);
        glEnd();
    glPopMatrix();
}

void drawQuadYZRight(float y, float z){
    glPushMatrix();
        glBegin(GL_QUADS);
            getNormal3p(0, 0, z, 0, -y, z, 0, -y, 0);

            glVertex3f(0, 0, z);
            glVertex3f(0, -y, z);
            glVertex3f(0, -y, 0);
            glVertex3f(0, 0, 0);
        glEnd();
    glPopMatrix();
}

void drawQuadXY(float x, float y){
    glPushMatrix();
        glBegin(GL_QUADS);
            getNormal3p(0, 0, 0, x, 0, 0, x, y, 0);

            glVertex3f(0, 0, 0);
            glVertex3f(x, 0, 0);
            glVertex3f(x, y, 0);
            glVertex3f(0, y, 0);
        glEnd();
    glPopMatrix();
}

void drawQuadXYDown(float x, float y){
    glPushMatrix();
        glBegin(GL_QUADS);
            getNormal3p(0, y, 0, x, y, 0, x, 0, 0);

            glVertex3f(0, y, 0);
            glVertex3f(x, y, 0);
            glVertex3f(x, 0, 0);
            glVertex3f(0, 0, 0);
        glEnd();
    glPopMatrix();
}

void drawButton(float x, float y, float z, float gap, vector<float> color){
    colorMaterial(color);
    glPushMatrix();
        glPushMatrix();
            glTranslatef(0, -y, 0);
            drawCube(x, y, z, color);
        glPopMatrix();

        glPushMatrix();
            colorMaterial({0,0,0});
            glTranslatef(-gap / 2, 0, -gap / 2);
            drawQuadXZ(x + gap, z +  gap);
        glPopMatrix();
    glPopMatrix();
}

void teleButtons(float buttonSize, float buttonGap, float height){

    glPushMatrix();
        glTranslatef(buttonGap, 0, 0);
        for(int i = 0; i < 3; i++){
            glTranslatef((i > 0) * (buttonSize + buttonGap), 0, 0);
            glPushMatrix();
            for(int j = 0; j < 5; j++){
                glTranslatef(0, 0, (j > 0) * (buttonSize + buttonGap));
                drawQuadXZ(buttonSize, buttonGap);
            }
            glPopMatrix();
        }
    glPopMatrix();

    glPushMatrix();
        for(int i = 0; i < 4; i++){
            glTranslatef((i > 0) * (buttonSize + buttonGap), 0, 0);
            drawQuadXZ(buttonGap, 5 * buttonGap + 4 * buttonSize);
        }
    glPopMatrix();

    float exGap = buttonSize / 15;

    glPushMatrix();
        for(int i = 0; i < 3; i++){
            glPushMatrix();
            glTranslatef((i + 1) * buttonGap + i * buttonSize + exGap, 0, 0);
            for(int j = 0; j < 4; j++){
                glPushMatrix();
                glTranslatef(0, height / 3, j * buttonSize + (j + 1) * buttonGap + exGap);
                drawButton(buttonSize - 2 * exGap, height, buttonSize - 2 * exGap, buttonGap, {!(i == 2 && j == 0) * 1, !(i == 0 && j == 0) * 1, !((i == 0 || i == 2) && (j == 0)) * 1});
                glPopMatrix();
            }
            glPopMatrix();
        }
    glPopMatrix();
}

void roundPlate(float r, float h){
    vector<vector<float>> p1;

    for(float theta = 0; theta <= 360; theta += 3){
        p1.push_back({r * cos(theta * M_PI / 180), 0, r * sin(theta * M_PI / 180)});
    }

    //plate
    glPushMatrix();
        glPushMatrix();
            glBegin(GL_TRIANGLES);
                for(int i = 0; i < p1.size() - 1; i++){
                    getNormal3p(0, 0, 0, p1[i][0], p1[i][1], p1[i][2], p1[i + 1][0], p1[i + 1][1], p1[i + 1][2]);

                    glVertex3f(0, 0, 0);
                    glVertex3fv(&p1[i][0]);
                    glVertex3fv(&p1[i + 1][0]);
                }
            glEnd();
        glPopMatrix();

        glPushMatrix();
            glBegin(GL_QUADS);
                for(int i = 0; i < p1.size() - 1; i++){
                    getNormal3p(p1[i][0], p1[i][1] - h, p1[i][2], p1[i][0], p1[i][1], p1[i][2], p1[i + 1][0], p1[i + 1][1], p1[i + 1][2]);

                    glVertex3f(p1[i][0], p1[i][1] - h, p1[i][2]);
                    glVertex3f(p1[i][0], p1[i][1], p1[i][2]);
                    glVertex3f(p1[i + 1][0], p1[i + 1][1], p1[i + 1][2]);
                    glVertex3f(p1[i + 1][0], p1[i + 1][1] - h, p1[i + 1][2]);
                }
            glEnd();
        glPopMatrix();
    glPopMatrix();
}

void draw_telephone_stand(float len, float r, float height, float ex){
    colorMaterial(getColor(171, 160, 19));
    vector<vector<float>> p1 = {
        {0, 0, 0}, {ex, 0, 0}, {ex, 0, len + 2 * ex + 2 * r}, {0, 0, len + 2 * ex + 2 * r},
        {ex, 0, 0}, {2 * r + ex, 0, 0}, {2 * r + ex, 0, ex}, {ex, 0, ex},
        {ex, 0, ex + 2*r}, {ex + 2*r, 0, ex + 2*r}, {ex + 2*r, 0, ex + 4*r}, {ex, 0, ex + 4*r},
        {ex, 0, ex + 6 * r}, {ex + 2 * r, 0, ex + 6 * r}, {ex + 2 * r, 0, 2 * ex + 6 * r}, {ex, 0, 2 * ex + 6 * r},
        {ex + 2 * r, 0, 0}, {2 * ex + 2 * r, 0, 0}, {2 * ex + 2 * r, 0, 2 * ex + 6 * r}, {ex + 2 * r, 0, 2 * ex + 6 * r}
    };

    glPushMatrix();

    glPushMatrix();
        glTranslatef(0, -height / 2, 0);
        glBegin(GL_QUADS);
            for(int i = 0; i < p1.size(); i += 4){
                getNormal3p(p1[i][0], p1[i][1], p1[i][2],
                            p1[i + 1][0], p1[i + 1][1], p1[i + 1][2],
                            p1[i + 2][0], p1[i + 2][1], p1[i + 2][2]);

                for(int j = 0; j < 4; j++){
                    glVertex3fv(&p1[i + j][0]);
                }
            }
        glEnd();

        glPushMatrix();
            glTranslatef(ex + r, 0, ex + r);
            square_with_circle_gap(r);
            glTranslatef(0, 0, 4 * r);
            square_with_circle_gap(r);
        glPopMatrix();
    glPopMatrix();

    float bs = r / 2, bg = r / 4, totalHeight = 2 * ex + 6 * r;

    glPushMatrix();
        glTranslatef(0, -height, 0);
        glPushMatrix();
            glTranslatef(2 * ex + 2 * r, 0, 0);
            drawQuadXZ(ex, 2 * ex + 6 * r);
        glPopMatrix();

        glPushMatrix();
            glTranslatef(3 * ex + 2 * r + 4 * bg + 3 * bs, 0, 0);
            drawQuadXZ(ex, 2 * ex + 6 * r);
        glPopMatrix();

        glPushMatrix();
            glTranslatef(3 * ex + 2 * r, 0, (totalHeight - (4 * bs  + 5 * bg)) / 2);
            teleButtons(bs, bg, bs / 4);
        glPopMatrix();

        colorMaterial(getColor(171, 160, 19));

        glPushMatrix();
            glTranslatef(3 * ex + 2 * r, 0, 0);
            drawQuadXZ(4 * bg + 3 * bs, (totalHeight - (4 * bs  + 5 * bg)) / 2);
        glPopMatrix();

        glPushMatrix();
            glTranslatef(3 * ex + 2 * r, 0, (totalHeight - (4 * bs  + 5 * bg)) / 2 + 4 * bs  + 5 * bg);
            drawQuadXZ(4 * bg + 3 * bs, (totalHeight - (4 * bs  + 5 * bg)) / 2);
        glPopMatrix();
    glPopMatrix();

    glPushMatrix();
        drawQuadXZBack(4 * ex + 2 * r + 4 * bg + 3 * bs, 2 * ex + 6 * r);
        drawQuadYZ(height / 2, 2 * ex + 6 * r);
        glPushMatrix();
            glTranslatef(2 * r + 2 * ex, - height / 2, 0);
            drawQuadYZ(height / 2, 2 * ex + 6 * r);
        glPopMatrix();
        glPushMatrix();
            glTranslatef(2 * r + 2 * ex + 4 * bg + 3 * bs + 2 * ex, 0, 0);
            drawQuadYZRight(height, 2 * ex + 6 * r);
        glPopMatrix();
    glPopMatrix();

    colorMaterial(getColor(237, 224, 57));
    glPushMatrix();
        glTranslatef(ex + r, -height / 4, ex + r);
        roundPlate(r, height / 4);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(ex + r, -height / 4, ex + 5 * r);
        roundPlate(r, height / 4);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0, -height / 2, 2 * ex + 6 * r);
        drawQuadXY(2 * ex + 2 * r, height / 2);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(2 * ex + 2 * r, -height, 2 * ex + 6 * r);
        drawQuadXY(2 * ex + 4 * bg + 3 * bs, height);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0, -height / 2, 0);
        drawQuadXYDown(2 * ex + 2 * r, height / 2);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(2 * ex + 2 * r, -height, 0);
        drawQuadXYDown(2 * ex + 4 * bg + 3 * bs, height);
    glPopMatrix();

    glPushMatrix();
        glTranslatef(ex + r, - height / 4, - 2.5 * len);
        drawCylinder(r * sin(3 * M_PI / 180), 2.5 * len, 10, 60, {0.2, 0.2, 0.2});
    glPopMatrix();

    glPopMatrix();
}

void draw_only_stair(GLfloat x, GLfloat y, GLfloat z, GLfloat ex, int n){
    GLfloat perX = (x - ex) / n;
    GLfloat perZ = z / n;

    colorMaterial(getColor(255, 255, 255));

    glPushMatrix();
        glEnable(GL_TEXTURE_2D);
        glPushMatrix();
        glBindTexture(GL_TEXTURE_2D, ID[0]);
        glBegin(GL_QUADS);
            for(int i = 0; i < n; i++){ // top Floor
                getNormal3p((i != 0) * ex + i * perX, 0, z - i * perZ,
                            ex + (i + 1) * perX + 0.05, 0, z - i * perZ,
                            ex + (i + 1) * perX + 0.05, y, z - i * perZ);

                glVertex3f((i != 0) * ex + i * perX, 0, z - i * perZ);                glTexCoord2f(0, 1);
                glVertex3f(ex + (i + 1) * perX + 0.05, 0, z - i * perZ);               glTexCoord2f(1, 1);
                glVertex3f(ex + (i + 1) * perX + 0.05, y, z - i * perZ);               glTexCoord2f(1, 0);
                glVertex3f((i != 0) * ex + i * perX, y, z - i * perZ);                glTexCoord2f(0, 0);
            }

            for(int i = 0; i < n; i++){ //Side front
                getNormal3p((i != 0) * ex + i * perX, 0, 0,
                            ex + (i + 1) * perX, 0, 0,
                            ex + (i + 1) * perX, 0, z - (i * perZ));

                glVertex3f((i != 0) * ex + i * perX, 0, 0);                 glTexCoord2f(0, 1);
                glVertex3f(ex + (i + 1) * perX, 0, 0);                      glTexCoord2f(1, 1);
                glVertex3f(ex + (i + 1) * perX, 0, z - (i * perZ));         glTexCoord2f(1, 0);
                glVertex3f((i != 0) * ex + i * perX, 0, z - (i * perZ));    glTexCoord2f(0, 0);
            }

            for(int i = 0; i < n; i++){ //Side back
                getNormal3p((i != 0) * ex + i * perX, y, z - (i * perZ),
                            ex + (i + 1) * perX, y, z - (i * perZ),
                            ex + (i + 1) * perX, y, 0);

                glVertex3f((i != 0) * ex + i * perX, y, z - (i * perZ));    glTexCoord2f(0, 1);
                glVertex3f(ex + (i + 1) * perX, y, z - (i * perZ));         glTexCoord2f(1, 1);
                glVertex3f(ex + (i + 1) * perX, y, 0);                      glTexCoord2f(1, 0);
                glVertex3f((i != 0) * ex + i * perX, y, 0);                 glTexCoord2f(0, 0);
            }
        glEnd();
        glPopMatrix();
        glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    colorMaterial(getColor(8, 4, 0));

    glPushMatrix();
        glBegin(GL_QUADS);
        for(int i = 0; i < n; i++){
            getNormal3p(ex + (i + 1) * perX, 0, z - (i + 1) * perZ,
                        ex + (i + 1) * perX, y, z - (i + 1) * perZ,
                        ex  + (i + 1) * perX, y, z - i * perZ);

            glVertex3f(ex + (i + 1) * perX, 0, z - (i + 1) * perZ);
            glVertex3f(ex + (i + 1) * perX, y, z - (i + 1) * perZ);
            glVertex3f(ex  + (i + 1) * perX, y, z - i * perZ);
            glVertex3f(ex  + (i + 1) * perX, 0, z - i * perZ);
        }
        glEnd();
    glPopMatrix();
}

void shelf(int n, float x, float y, float z){
    glPushMatrix();
    glTranslatef(n * x, y, 0);
    glRotatef(180, 0, 0, 1);
    for(int i = 0; i < n; i++){
        glPushMatrix();
            glTranslatef(i * x, 0, 0);
            drawBoxThatCanOpen(x, y, z, 0);
        glPopMatrix();
    }
    glPopMatrix();
}

void plate(float R, float r, float d, vector<float> color){
    vector<vector<float>> UpOut, UpIn, DownOut, DownIn;
    colorMaterial(color);

    float h = (R - r) * tan(20 * M_PI / 180);
    float Dx = d / sin(20 * M_PI / 180);
    float dx = d * cos(20 * M_PI / 180);
    float dz = d * sin(20 * M_PI / 180);

    float URout = R, URin = R - Dx, DRout = r, DRin = r - dx;

    for(float theta = 0; theta <= 360; theta++){
        UpOut.push_back({URout * cos(theta * M_PI / 180), URout * sin(theta * M_PI / 180), h});
        UpIn.push_back({URin * cos(theta * M_PI / 180), URin * sin(theta * M_PI / 180), h});
        DownOut.push_back({DRout * cos(theta * M_PI / 180), DRout * sin(theta * M_PI / 180), 0});
        DownIn.push_back({DRin * cos(theta * M_PI / 180), DRin * sin(theta * M_PI / 180), dz});
    }

    glPushMatrix();
        glBegin(GL_TRIANGLES);
            for(int i = 0; i < DownOut.size()-1; i++){
                getNormal3p(0, 0, 0, DownOut[i + 1][0], DownOut[i + 1][1], DownOut[i + 1][2], DownOut[i][0], DownOut[i][1], DownOut[i][2]);
                glVertex3f(0, 0, 0);
                glVertex3f(DownOut[i + 1][0], DownOut[i + 1][1], DownOut[i + 1][2]);
                glVertex3f(DownOut[i][0], DownOut[i][1], DownOut[i][2]);
            }

            for(int i = 0; i < DownIn.size()-1; i++){
                getNormal3p(0, 0, 0, DownIn[i][0], DownIn[i][1], DownIn[i][2], DownIn[i + 1][0], DownIn[i + 1][1], DownIn[i + 1][2]);
                glVertex3f(0, 0, dz);
                glVertex3f(DownIn[i][0], DownIn[i][1], DownIn[i][2]);
                glVertex3f(DownIn[i + 1][0], DownIn[i + 1][1], DownIn[i + 1][2]);
            }
        glEnd();
    glPopMatrix();

    glPushMatrix();
        glBegin(GL_QUADS);
            for(int i = 0; i < DownOut.size() - 1; i++){
                getNormal3p(DownOut[i][0], DownOut[i][1], DownOut[i][2], DownOut[i + 1][0], DownOut[i + 1][1], DownOut[i + 1][2], UpOut[i + 1][0], UpOut[i + 1][1], UpOut[i + 1][2]);

                glVertex3f(DownOut[i][0], DownOut[i][1], DownOut[i][2]);
                glVertex3f(DownOut[i + 1][0], DownOut[i + 1][1], DownOut[i + 1][2]);
                glVertex3f(UpOut[i + 1][0], UpOut[i + 1][1], UpOut[i + 1][2]);
                glVertex3f(UpOut[i][0], UpOut[i][1], UpOut[i][2]);
            }

            for(int i = 0; i < DownIn.size() - 1; i++){
                getNormal3p(UpIn[i][0], UpIn[i][1], UpIn[i][2], UpIn[i + 1][0], UpIn[i + 1][1], UpIn[i + 1][2], DownIn[i + 1][0], DownIn[i + 1][1], DownIn[i + 1][2]);

                glVertex3f(UpIn[i][0], UpIn[i][1], UpIn[i][2]);
                glVertex3f(UpIn[i + 1][0], UpIn[i + 1][1], UpIn[i + 1][2]);
                glVertex3f(DownIn[i + 1][0], DownIn[i + 1][1], DownIn[i + 1][2]);
                glVertex3f(DownIn[i][0], DownIn[i][1], DownIn[i][2]);
            }

            for(int i = 0; i < UpIn.size() - 1; i++){
                getNormal3p(UpIn[i][0], UpIn[i][1], UpIn[i][2], UpOut[i][0], UpOut[i][1], UpOut[i][2], UpOut[i + 1][0], UpOut[i + 1][1], UpOut[i + 1][2]);

                glVertex3f(UpIn[i][0], UpIn[i][1], UpIn[i][2]);
                glVertex3f(UpOut[i][0], UpOut[i][1], UpOut[i][2]);
                glVertex3f(UpOut[i + 1][0], UpOut[i + 1][1], UpOut[i + 1][2]);
                glVertex3f(UpIn[i + 1][0], UpIn[i + 1][1], UpIn[i + 1][2]);
            }
        glEnd();
    glPopMatrix();
}

void forkErKata(float R, float dy, float dz){
    vector<vector<float>> p1;

    for(float theta = 90; theta <= 90 + 60; theta++){
        p1.push_back({R * cos(theta * M_PI / 180), 0, R * sin(theta * M_PI / 180)});
    }

    glPushMatrix();
        glTranslatef(0, 0, -R);
        glBegin(GL_QUADS);
            for(int i = 0; i < p1.size() - 1; i++){
                getNormal3p(p1[i][0], p1[i][1], p1[i][2], p1[i + 1][0], p1[i + 1][1], p1[i + 1][2], p1[i + 1][0], p1[i + 1][1] + dy, p1[i + 1][2]);

                glVertex3f(p1[i][0], p1[i][1], p1[i][2]);
                glVertex3f(p1[i + 1][0], p1[i + 1][1], p1[i + 1][2]);
                glVertex3f(p1[i + 1][0], p1[i + 1][1] + dy, p1[i + 1][2]);
                glVertex3f(p1[i][0], p1[i][1] + dy, p1[i][2]);
            }

            for(int i = 0; i < p1.size() - 1; i++){
                getNormal3p(p1[i][0], p1[i][1] + dy, p1[i][2] + dz, p1[i + 1][0], p1[i + 1][1] + dy, p1[i + 1][2] + dz, p1[i + 1][0], p1[i + 1][1], p1[i + 1][2] + dz);

                glVertex3f(p1[i][0], p1[i][1] + dy, p1[i][2] + dz);
                glVertex3f(p1[i + 1][0], p1[i + 1][1] + dy, p1[i + 1][2] + dz);
                glVertex3f(p1[i + 1][0], p1[i + 1][1], p1[i + 1][2] + dz);
                glVertex3f(p1[i][0], p1[i][1], p1[i][2] + dz);
            }

            for(int i = 0; i < p1.size() - 1; i++){
                getNormal3p(p1[i][0], p1[i][1], p1[i][2], p1[i][0], p1[i][1], p1[i][2] + dz, p1[i + 1][0], p1[i + 1][1], p1[i + 1][2] + dz);

                glVertex3f(p1[i][0], p1[i][1], p1[i][2]);
                glVertex3f(p1[i][0], p1[i][1], p1[i][2] + dz);
                glVertex3f(p1[i + 1][0], p1[i + 1][1], p1[i + 1][2] + dz);
                glVertex3f(p1[i + 1][0], p1[i + 1][1], p1[i + 1][2]);
            }

            for(int i = 0; i < p1.size() - 1; i++){
                getNormal3p(p1[i + 1][0], p1[i + 1][1] + dy, p1[i + 1][2], p1[i + 1][0], p1[i + 1][1] + dy, p1[i + 1][2] + dz, p1[i][0], p1[i][1] + dy, p1[i][2] + dz);

                glVertex3f(p1[i + 1][0], p1[i + 1][1] + dy, p1[i + 1][2]);
                glVertex3f(p1[i + 1][0], p1[i + 1][1] + dy, p1[i + 1][2] + dz);
                glVertex3f(p1[i][0], p1[i][1] + dy, p1[i][2] + dz);
                glVertex3f(p1[i][0], p1[i][1] + dy, p1[i][2]);
            }

            getNormal3p(p1.back()[0], p1.back()[1], p1.back()[2], p1.back()[0], p1.back()[1], p1.back()[2] + dz, p1.back()[0], p1.back()[1] + dy, p1.back()[2] + dz);
            glVertex3f(p1.back()[0], p1.back()[1], p1.back()[2]);
            glVertex3f(p1.back()[0], p1.back()[1], p1.back()[2] + dz);
            glVertex3f(p1.back()[0], p1.back()[1] + dy, p1.back()[2] + dz);
            glVertex3f(p1.back()[0], p1.back()[1] + dy, p1.back()[2]);

            getNormal3p(p1[0][0], p1[0][1] + dy, p1[0][2], p1[0][0], p1[0][1] + dy, p1[0][2] + dz, p1[0][0], p1[0][1], p1[0][2] + dz);
            glVertex3f(p1[0][0], p1[0][1] + dy, p1[0][2]);
            glVertex3f(p1[0][0], p1[0][1] + dy, p1[0][2] + dz);
            glVertex3f(p1[0][0], p1[0][1], p1[0][2] + dz);
            glVertex3f(p1[0][0], p1[0][1], p1[0][2]);
        glEnd();
    glPopMatrix();
}

void forkPlate(float r, float dz){
    vector<vector<float>> p1;
    for(float theta = -90; theta <= 90; theta++){
        p1.push_back({r * cos(theta * M_PI / 180), r * sin(theta * M_PI / 180), 0});
    }

    glPushMatrix();
        glBegin(GL_TRIANGLES);
            for(int i = 0; i < p1.size() - 1; i++){
                getNormal3p(0, 0, 0, p1[i + 1][0], p1[i + 1][1], p1[i + 1][2], p1[i][0], p1[i][1], p1[i][2]);
                glVertex3f(0, 0, 0);
                glVertex3f(p1[i + 1][0], p1[i + 1][1], p1[i + 1][2]);
                glVertex3f(p1[i][0], p1[i][1], p1[i][2]);
            }
            for(int i = 0; i < p1.size() - 1; i++){
                getNormal3p(0, 0, 0, p1[i][0], p1[i][1], p1[i][2] + dz, p1[i + 1][0], p1[i + 1][1], p1[i + 1][2] + dz);
                glVertex3f(0, 0, 0);
                glVertex3f(p1[i][0], p1[i][1], p1[i][2] + dz);
                glVertex3f(p1[i + 1][0], p1[i + 1][1], p1[i + 1][2] + dz);
            }
        glEnd();

        glBegin(GL_QUADS);
            for(int i = 0; i < p1.size() - 1; i++){
                getNormal3p(p1[i][0], p1[i][1], p1[i][2], p1[i + 1][0], p1[i + 1][1], p1[i + 1][2], p1[i + 1][0], p1[i + 1][1], p1[i + 1][2] + dz);
                glVertex3f(p1[i][0], p1[i][1], p1[i][2]);
                glVertex3f(p1[i + 1][0], p1[i + 1][1], p1[i + 1][2]);
                glVertex3f(p1[i + 1][0], p1[i + 1][1], p1[i + 1][2] + dz);
                glVertex3f(p1[i][0], p1[i][1], p1[i][2] + dz);
            }

            getNormal3p(p1[0][0], p1[0][1], p1[0][2], p1[0][0], p1[0][1], p1[0][2] + dz, p1.back()[0], p1.back()[1], p1.back()[2] + dz);
            glVertex3f(p1[0][0], p1[0][1], p1[0][2]);
            glVertex3f(p1[0][0], p1[0][1], p1[0][2] + dz);
            glVertex3f(p1.back()[0], p1.back()[1], p1.back()[2] + dz);
            glVertex3f(p1.back()[0], p1.back()[1], p1.back()[2]);
        glEnd();
    glPopMatrix();
}

void draw_fork(float R, float dy, float dz, int nKata, float gap, float hanW, float hanH, vector<float> color){
    colorMaterial(color);
    glPushMatrix();
    glTranslatef(0, 0, hanH + (nKata * dy + (nKata - 1) * gap));
    glRotatef(90, 0, 1, 0);

    glPushMatrix();
        for(int i = 0; i < nKata; i++){
            glPushMatrix();
                glTranslatef(0, i * gap + i * dy, 0);
                forkErKata(R, dy, dz);
            glPopMatrix();
        }
    glPopMatrix();

    float r = (nKata * dy + (nKata - 1) * gap);
    glPushMatrix();
        glBegin(GL_QUADS);
            getNormal3p(0, 0, 0, 0, r, 0, r, (r - hanW) / 2 + hanW, 0);
            glVertex3f(0, 0, 0);
            glVertex3f(0, r, 0);
            glVertex3f(r, (r - hanW) / 2 + hanW, 0);
            glVertex3f(r, (r - hanW) / 2, 0);

            getNormal3p(r, (r - hanW) / 2, dz, r, (r - hanW) / 2 + hanW, dz, 0, r, dz);
            glVertex3f(r, (r - hanW) / 2, dz);
            glVertex3f(r, (r - hanW) / 2 + hanW, dz);
            glVertex3f(0, r, dz);
            glVertex3f(0, 0, dz);

            getNormal3p(0, 0, 0, r, (r - hanW) / 2, 0, r, (r - hanW) / 2, dz);
            glVertex3f(0, 0, 0);
            glVertex3f(r, (r - hanW) / 2, 0);
            glVertex3f(r, (r - hanW) / 2, dz);
            glVertex3f(0, 0, dz);

            getNormal3p(0, r, dz, r, (r - hanW) / 2 + hanW, dz, r, (r - hanW) / 2 + hanW, 0);
            glVertex3f(0, r, dz);
            glVertex3f(r, (r - hanW) / 2 + hanW, dz);
            glVertex3f(r, (r - hanW) / 2 + hanW, 0);
            glVertex3f(0, r, 0);

            getNormal3p(r, (r -hanW) / 2, dz, r, (r - hanW) / 2, 0, r, (r - hanW) / 2 + hanW, 0);
            glVertex3f(r, (r -hanW) / 2, dz);
            glVertex3f(r, (r - hanW) / 2, 0);
            glVertex3f(r, (r - hanW) / 2 + hanW, 0);
            glVertex3f(r, (r - hanW) / 2 + hanW, dz);

            getNormal3p(0, 0, dz, 0, r, dz, 0, r, 0);
            glVertex3f(0, 0, dz);
            glVertex3f(0, r, dz);
            glVertex3f(0, r, 0);
            glVertex3f(0, 0, 0);
        glEnd();
    glPopMatrix();

    glPushMatrix();
        glTranslatef(r, (r - hanW)/2, 0);
        drawCube(hanH, hanW, dz, color);
    glPopMatrix();

    glPopMatrix();
}

void drawUp(float x, float y, float z, int ti = 0){
    colorMaterial({1, 1, 1});
    vector<vector<GLfloat>> points = generatePoints8p(x, y, z);
    vector<vector<GLubyte>> ind = getMyIndicesForCube();

    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glBindTexture(GL_TEXTURE_2D, ID[ti]);
    glBegin(GL_QUADS);

        for(int i = 0; i < 6; i++){

            getNormal3p(
                        points[ind[i][0]][0], points[ind[i][0]][1], points[ind[i][0]][2],
                        points[ind[i][1]][0], points[ind[i][1]][1], points[ind[i][1]][2],
                        points[ind[i][2]][0], points[ind[i][2]][1], points[ind[i][2]][2]
            );

            for(int j = 0; j < 4; j++){
                glVertex3fv(&points[ind[i][j]][0]);
                if(j == 0) glTexCoord2f(0,0);
                if(j == 1) glTexCoord2f(0,15.0 / 27.0);
                if(j == 2) glTexCoord2f(1, 15.0 / 27.0);
                if(j == 3) glTexCoord2f(1,0);
            }
        }

    glEnd();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
}

void drawUpPart2(float x, float y, float z, int ti = 0){
    colorMaterial({1, 1, 1});
    vector<vector<GLfloat>> points = generatePoints8p(x, y, z);
    vector<vector<GLubyte>> ind = getMyIndicesForCube();

    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    glBindTexture(GL_TEXTURE_2D, ID[ti]);
    glBegin(GL_QUADS);

        for(int i = 0; i < 6; i++){

            getNormal3p(
                        points[ind[i][0]][0], points[ind[i][0]][1], points[ind[i][0]][2],
                        points[ind[i][1]][0], points[ind[i][1]][1], points[ind[i][1]][2],
                        points[ind[i][2]][0], points[ind[i][2]][1], points[ind[i][2]][2]
            );

            for(int j = 0; j < 4; j++){
                glVertex3fv(&points[ind[i][j]][0]);
                if(j == 0) glTexCoord2f(0,15.0 / 27);
                if(j == 1) glTexCoord2f(0, 1);
                if(j == 2) glTexCoord2f(14.0 / 37, 1);
                if(j == 3) glTexCoord2f(14.0 / 37, 15.0 / 27);
            }
        }

    glEnd();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
}

void draw_first_floor(float x, float y, float z){
    glPushMatrix();
        draw_cube_textured(x, 0.5, z, {0, 0, 0}, 0); // left side
        glPushMatrix();
            glTranslatef(0, 0.5, 0);
            draw_cube_textured(0.5, y - 0.5, z, {0,0,0}, 0); // back side
        glPopMatrix();
        glPushMatrix();
            glTranslatef(0.5, y - 0.5, 0);
            draw_cube_textured(x - 0.5, 0.5, z, {0, 0, 0}, 0); // right side
        glPopMatrix();

        glPushMatrix(); // front wall of two room
            glTranslatef(15, 0.5, 0);
            draw_cube_textured(0.5, 14, z, {0, 0, 0}, 0); // Be alert here fixed size provided
            glTranslatef(0, y - 15, 0);
            draw_cube_textured(0.5, 14, z, {0, 0, 0}, 0);
        glPopMatrix();

        glPushMatrix(); // front wall
            glTranslatef(27 - 0.5, 0.5, 0);
            draw_cube_textured(0.5, 14, z, {0, 0, 0}, 0);
        glPopMatrix();

        glPushMatrix();
            glTranslatef(0, 0, z);
            drawUp(15, 37, 0.5);
            glTranslatef(15, 0, 0);
            drawUpPart2(12, 14.5, 0.5);
        glPopMatrix();

        glPushMatrix();
            glTranslatef(15 + (27 - 15) / 2, 7, z);
            myLight(0.4, 0.3, {1, 1, 1}, Mainlight);
        glPopMatrix();

        glPushMatrix(); // stair
            glTranslatef(15 + 0.5, 36 + 0.5 - 5, 0);
            draw_only_stair(8, 5, z, 0, 8);
        glPopMatrix();

        glPushMatrix(); // telephone
            glTranslatef(16, 36 - 12, 0);
            glRotatef(125, 0, 0, 1);
            draw_telephone(0.75, 0.75 / 4);
        glPopMatrix();

        glPushMatrix(); // telephone stand
            glTranslatef(15.5, 36 - 13, 4.5);
            glRotatef(90, 0, 0, 1);
            draw_telephone_stand(0.70, 0.70 / 4, 0.4, 0.1);
        glPopMatrix();

        glPushMatrix(); // Bed
            glTranslatef(0.5 + 9 + 0.75, 0.5, 0);
            glRotatef(90, 0, 0, 1);
            drawBed(9, 6, 1);
        glPopMatrix();

        glPushMatrix();
            glTranslatef(0.5, 0.5 + 7, 0);
            glPushMatrix(); //shelf
                glTranslatef(2, 0, 0);
                glRotatef(90, 0, 0, 1);
                drawBoxThatCanOpen(2, 2, 3, 1);
            glPopMatrix();

            glPushMatrix(); // lamp
                glTranslatef(1, 1, 3);
                drawLamp(3, 1, 0.4, 1, 0.03, 0.75, 0.05);
            glPopMatrix();
        glPopMatrix();

        glPushMatrix(); // Cooking shelf
            glTranslatef(15.5, 0.5, 0);
            shelf(3, 3, 3, 3.5);
        glPopMatrix();

        glPushMatrix();
            glTranslatef(20, 2, 3.5);
            //glRotatef(-80, 1, 0, 0);
            plate(0.5, 0.3, 0.02, getColor(242, 216, 240));
            glTranslatef(2, 0, 0);
            plate(0.5, 0.3, 0.02, getColor(242, 216, 240));
        glPopMatrix();

        glPushMatrix();
            glTranslatef(19, 2, 3.5);
            glRotatef(90, 1, 0, 0);
            draw_fork(0.2, 0.02, 0.01, 4, 0.02, 0.03, 0.4, getColor(48, 48, 48));
            glTranslatef(2, 0, 0);
            draw_fork(0.2, 0.02, 0.01, 4, 0.02, 0.03, 0.4, getColor(48, 48, 48));
        glPopMatrix();

        glPushMatrix();
            glTranslatef(0.5 + 0.1, 18.5, 6);
            glRotatef(90, 0, 1, 0);
            drawClock(2);
        glPopMatrix();

//        glPushMatrix();
//            glTranslatef(10 + 0.5 + 0.1, 18.5, 6);
//            glRotatef(90, 0, 1, 0);
//            drawPhoto(2);
//        glPopMatrix();
    glPopMatrix();
}

void play(){
    glPushMatrix();
        draw_area(53, 50, {0, 0.4, 0});
    glPopMatrix();

    glPushMatrix();
        glTranslatef(3, 3, 0);
        draw_basement();
        glTranslatef(4, 4, 3);
        draw_first_floor(27, 37, 8);
    glPopMatrix();
}

class Monster{
public:
    float R, H, r, h, tx, ty, tz;
    float eyeAliveRat = 1.0 / 40, eyeDeadRat = 1.0 / 20;
    float rotTail;
    bool alive;
    bool shoot;
    bool shooting;
    float shootingStep;
    float currStep;
    bool inc, dec;
    float oR;
    point top;
    point bottom;
    point face;
    float ms = 0.3;
    float health;

    myVector up, look, left;
    // R - radius of monster, H = height of Monster (cylindrical part), r - radius of arms, h - height of arms
    Monster(){
        R = 0, H = 0, r = 0, h = 0;
        look = myVector(0, -1, 0);
        up = myVector(0, 0, 1);
        left = myVector(1, 0, 0);
        alive = true;
        rotTail = 0;
        shoot = false;
        shooting = false;
    }

    Monster(float R, float H, float r, float h, float tx, float ty, float tz){
        this->R = R;
        this->r = r;
        this->H = H;
        this->h = h;
        this->tx = tx;
        this->ty = ty;
        this->tz = tz;
        look = myVector(0, -1, 0);
        up = myVector(0, 0, 1);
        left = myVector(1, 0, 0);
        alive = true;
        rotTail = 0;
        shoot = false;
        shooting = false;
        shootingStep = h / 3;
        currStep = 0;
        inc = 0, dec = 0;
        face = point(0 + tx, 0 + ty - R, H + tz);
        top = point(0 + tx, 0 + ty, H + tz);
        bottom = point(0 + tx, 0 + ty, 0 + tz);
        oR = R - 2 * R / 15;
        health = 100;
    }

    void drawCircularLine(float rad, float from, float to, vector<float> color, bool glow = false){
        if(glow) colorMaterialEmit(color);
        else colorMaterial(color);

        glPushMatrix();
            glTranslatef(-rad * cos(to * M_PI / 180), 0, 0);
            glBegin(GL_LINES);
                for(float theta = from; theta + 1 <= to; theta++){
                    glVertex3f(rad * cos(theta * M_PI / 180), 0, rad * sin(theta * M_PI / 180));
                    glVertex3f(rad * cos((theta + 1) * M_PI / 180), 0, rad * sin((theta + 1) * M_PI / 180));
                }
            glEnd();
        glPopMatrix();
    }

    void drawBody(){
        colorMaterial(getColor(11, 1, 33));

        vector<point> p1;
        for(float theta = 0; theta <= 360; theta++){
            p1.push_back(point(R * cos(theta * M_PI / 180), R * sin(theta * M_PI / 180), H));
        }

        float eyeH = H * (alive ? eyeAliveRat : eyeDeadRat);
        float eyePos = H - H / 6;
        float oR = R - 2 * R / 15;

        vector<point> p2;
        for(float theta = 0; theta <= 360; theta++){
            p2.push_back(point(oR * cos(theta * M_PI / 180), oR * sin(theta * M_PI / 180), H));
        }

        glPushMatrix();
            glBegin(GL_QUADS);
                for(int i = 0; i < p1.size() - 1; i++){
                    getNormal3p(p1[i].x, p1[i].y, p1[i].z, p1[i].x, p1[i].y, eyePos, p1[i + 1].x, p1[i + 1].y, eyePos);
                    glVertex3f(p1[i].x, p1[i].y, p1[i].z);
                    glVertex3f(p1[i].x, p1[i].y, eyePos);
                    glVertex3f(p1[i + 1].x, p1[i + 1].y, eyePos);
                    glVertex3f(p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);
                }

                for(int i = 0; i < p1.size(); i++){
                    if((i >= 235 && i + 1 <= 265) || (i >= 275 && i + 1 <= 305)){
                        glPushMatrix();
                            colorMaterialEmit(getColor(242, 255, 0));
                            getNormal3p(p1[i].x, p1[i].y, eyePos, p1[i].x, p1[i].y, eyePos - eyeH, p1[i + 1].x, p1[i + 1].y, eyePos - eyeH);

                            glVertex3f(p1[i].x, p1[i].y, eyePos);
                            glVertex3f(p1[i].x, p1[i].y, eyePos - eyeH);
                            glVertex3f(p1[i + 1].x, p1[i + 1].y, eyePos - eyeH);
                            glVertex3f(p1[i + 1].x, p1[i + 1].y, eyePos);
                        glPopMatrix();
                    }
                    else{
                        glPushMatrix();
                            colorMaterial(getColor(11, 1, 33));
                            getNormal3p(p1[i].x, p1[i].y, eyePos, p1[i].x, p1[i].y, eyePos - eyeH, p1[i + 1].x, p1[i + 1].y, eyePos - eyeH);

                            glVertex3f(p1[i].x, p1[i].y, eyePos);
                            glVertex3f(p1[i].x, p1[i].y, eyePos - eyeH);
                            glVertex3f(p1[i + 1].x, p1[i + 1].y, eyePos - eyeH);
                            glVertex3f(p1[i + 1].x, p1[i + 1].y, eyePos);
                        glPopMatrix();
                    }
                }

                colorMaterial(getColor(11, 1, 33));
                for(int i = 0; i < p1.size() - 1; i++){
                    getNormal3p(p1[i].x, p1[i].y, eyePos - eyeH, p1[i].x, p1[i].y, 0, p1[i + 1].x, p1[i + 1].y, 0);
                    glVertex3f(p1[i].x, p1[i].y, eyePos - eyeH);
                    glVertex3f(p1[i].x, p1[i].y, 0);
                    glVertex3f(p1[i + 1].x, p1[i + 1].y, 0);
                    glVertex3f(p1[i + 1].x, p1[i + 1].y, eyePos - eyeH);
                }

                for(int i = 0; i < p1.size() - 1; i++){
                    getNormal3p(p1[i].x, p1[i].y, p1[i].z, p1[i + 1].x, p1[i + 1].y, p1[i + 1].z, p2[i + 1].x, p2[i + 1].y, p2[i + 1].z);
                    glVertex3f(p1[i].x, p1[i].y, p1[i].z);
                    glVertex3f(p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);
                    glVertex3f(p2[i + 1].x, p2[i + 1].y, p2[i + 1].z);
                    glVertex3f(p2[i].x, p2[i].y, p2[i].z);
                }
            glEnd();
        glPopMatrix();

        glPushMatrix();
            glBegin(GL_TRIANGLES);
                for(int i = 0; i < p1.size() - 1; i++){
                    getNormal3p(0, 0, 0, p1[i + 1].x, p1[i + 1].y, 0, p1[i].x, p1[i].y, 0);
                    glVertex3f(0, 0, 0);
                    glVertex3f(p1[i + 1].x, p1[i + 1].y, 0);
                    glVertex3f(p1[i].x, p1[i].y, 0);
                }
            glEnd();
        glPopMatrix();

        glPushMatrix();
            glTranslatef(0, 0, H);
            drawHalfSphere(oR, 80, 80, getColor(242, 255, 0), true);
        glPopMatrix();

        glPushMatrix();
            glRotatef(rotTail, 0, 0, 1);
            glTranslatef(0, 0, -oR * sin(M_PI / 3));
            for(float i = 0; i < 360; i += 24){
                glPushMatrix();
                    glRotatef(i, 0, 0, 1);
                    drawCircularLine(oR, 0, 60, getColor(242, 255, 0), true);
                glPopMatrix();
            }
        glPopMatrix();
    }

    vector<point> generatePXZ(float rad){
        vector<point> p1;
        for(float theta = 0; theta <= 360; theta += 3){
            p1.push_back(point(rad * cos(theta * M_PI / 180), 0, rad * sin(theta * M_PI / 180)));
        }

        return p1;
    }

    void drawCylinderXZOut(vector<point> p1, vector<float> color, bool glow = false){
        if(glow) colorMaterialEmit(color);
        else colorMaterial(color);

        glPushMatrix();
        glBegin(GL_QUADS);
            for(int i = 0; i < p1.size() - 1; i++){
                getNormal3p(p1[i].x, p1[i].y - h, p1[i].z, p1[i].x, p1[i].y, p1[i].z, p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);

                glVertex3f(p1[i].x, p1[i].y - h, p1[i].z);
                glVertex3f(p1[i].x, p1[i].y, p1[i].z);
                glVertex3f(p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);
                glVertex3f(p1[i + 1].x, p1[i + 1].y - h, p1[i + 1].z);
            }
        glEnd();
        glPopMatrix();
    }

    void drawCylinderXZIn(vector<point> p1, vector<float> color, bool glow = false){
        if(glow) colorMaterialEmit(color);
        else colorMaterial(color);

        glPushMatrix();
        glBegin(GL_QUADS);
            for(int i = 0; i < p1.size() - 1; i++){
                getNormal3p(p1[i + 1].x, p1[i + 1].y - h, p1[i + 1].z, p1[i + 1].x, p1[i + 1].y, p1[i + 1].z, p1[i].x, p1[i].y, p1[i].z);

                glVertex3f(p1[i + 1].x, p1[i + 1].y - h, p1[i + 1].z);
                glVertex3f(p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);
                glVertex3f(p1[i].x, p1[i].y, p1[i].z);
                glVertex3f(p1[i].x, p1[i].y - h, p1[i].z);
            }
        glEnd();
        glPopMatrix();
    }

    void drawCircleXZFront(vector<point> p1, vector<float> color, bool glow = false){
        if(glow) colorMaterialEmit(color);
        else colorMaterial(color);

        glPushMatrix();
            glBegin(GL_TRIANGLES);
                for(int i = 0; i < p1.size() - 1; i++){
                    getNormal3p(0, -h, 0, p1[i].x, p1[i].y - h, p1[i].z, p1[i + 1].x, p1[i + 1].y - h, p1[i + 1].z);
                    glVertex3f(0, -h, 0);
                    glVertex3f(p1[i].x, p1[i].y - h, p1[i].z);
                    glVertex3f(p1[i + 1].x, p1[i + 1].y - h, p1[i + 1].z);
                }
            glEnd();
        glPopMatrix();
    }

    void drawCircleXZBack(vector<point> p1, vector<float> color, bool glow = false){
        if(glow) colorMaterialEmit(color);
        else colorMaterial(color);

        glPushMatrix();
            glBegin(GL_TRIANGLES);
                for(int i = 0; i < p1.size() - 1; i++){
                    getNormal3p(0, 0, 0, p1[i + 1].x, p1[i + 1].y, p1[i + 1].z, p1[i].x, p1[i].y, p1[i].z);
                    glVertex3f(0, 0, 0);
                    glVertex3f(p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);
                    glVertex3f(p1[i].x, p1[i].y, p1[i].z);
                }
            glEnd();
        glPopMatrix();
    }
    void drawConeXZ(vector<point> p1, float y, vector<float> color, bool glow = false){
        if(glow) colorMaterialEmit(color);
        else colorMaterial(color);

        glPushMatrix();
            glBegin(GL_TRIANGLES);
                for(int i = 0; i < p1.size() - 1; i++){
                    getNormal3p(0, -y, 0, p1[i].x, p1[i].y, p1[i].z, p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);
                    glVertex3f(0, -y, 0);
                    glVertex3f(p1[i].x, p1[i].y, p1[i].z);
                    glVertex3f(p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);
                }
            glEnd();
        glPopMatrix();
    }
    void drawPartialCircleFront(vector<point>p1, vector<point> p2, vector<float> color, bool glow = false){
        if(glow) colorMaterialEmit(color);
        else colorMaterial(color);

        glPushMatrix();
            glBegin(GL_QUADS);
                for(int i = 0; i < p1.size() - 1; i++){
                    getNormal3p(p1[i].x, p1[i].y, p1[i].z, p2[i].x, p2[i].y, p2[i].z, p2[i + 1].x, p2[i + 1].y, p2[i + 1].z);
                    glVertex3f(p1[i].x, p1[i].y, p1[i].z);
                    glVertex3f(p2[i].x, p2[i].y, p2[i].z);
                    glVertex3f(p2[i + 1].x, p2[i + 1].y, p2[i + 1].z);
                    glVertex3f(p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);
                }
            glEnd();
        glPopMatrix();
    }

    void drawPartialCircleBack(vector<point>p1, vector<point> p2, vector<float> color, bool glow = false){
        if(glow) colorMaterialEmit(color);
        else colorMaterial(color);

        glPushMatrix();
            glBegin(GL_QUADS);
                for(int i = 0; i < p1.size() - 1; i++){
                    getNormal3p(p1[i + 1].x, p1[i + 1].y, p1[i + 1].z, p2[i + 1].x, p2[i + 1].y, p2[i + 1].z, p2[i].x, p2[i].y, p2[i].z);

                    glVertex3f(p1[i + 1].x, p1[i + 1].y, p1[i + 1].z);
                    glVertex3f(p2[i + 1].x, p2[i + 1].y, p2[i + 1].z);
                    glVertex3f(p2[i].x, p2[i].y, p2[i].z);
                    glVertex3f(p1[i].x, p1[i].y, p1[i].z);
                }
            glEnd();
        glPopMatrix();
    }

    void arm(){
        float curr = r / 4;
        vector<point> p1;
        p1 = generatePXZ(curr);

        glPushMatrix(); //Layer 0
            if(alive) glTranslatef(0, -(shooting * currStep * 3), 0);
            drawCircleXZBack(p1, getColor(11, 1, 33), false);
            drawCircleXZFront(p1, getColor(11, 1, 33), false);
            drawCylinderXZOut(p1, getColor(242, 255, 0), true);
            glPushMatrix();
                glTranslatef(0, -h, 0);
                drawConeXZ(p1, h / 2, getColor(242, 255, 0), true);
            glPopMatrix();
        glPopMatrix();

        curr = r / 2;
        vector<point> p2 = generatePXZ(curr); //Layer 1

        glPushMatrix();
            if(alive) glTranslatef(0, -(shooting * currStep * 2), 0);
            drawCylinderXZOut(p2, getColor(242, 255, 0), true);
            drawPartialCircleBack(p1, p2, getColor(11, 1, 33), false);
            glPushMatrix();
            glTranslatef(0, -h, 0);
            drawPartialCircleFront(p1, p2, getColor(11, 1, 33), false);
            glPopMatrix();
        glPopMatrix();

        p1.clear();
        p1 = p2;
        curr = 3 * r / 4;
        p2 = generatePXZ(curr);

        glPushMatrix(); // Layer 2
            if(alive) glTranslatef(0, -(shooting * currStep), 0);
            drawCylinderXZOut(p2, getColor(242, 255, 0), true);
            drawPartialCircleBack(p1, p2, getColor(11, 1, 33), false);
            glPushMatrix();
            glTranslatef(0, -h, 0);
            drawPartialCircleFront(p1, p2, getColor(11, 1, 33), false);
            glPopMatrix();
        glPopMatrix();

        p1.clear();
        p1 = p2;
        curr = r;
        p2 = generatePXZ(curr);

        glPushMatrix(); // Layer 3
            drawCylinderXZOut(p2, getColor(242, 255, 0), true);
            drawPartialCircleBack(p1, p2, getColor(11, 1, 33), false);
            glPushMatrix();
            glTranslatef(0, -h, 0);
            drawPartialCircleFront(p1, p2, getColor(11, 1, 33), false);
            glPopMatrix();
        glPopMatrix();

        if(shooting){
            if(inc){
                if(currStep < h) currStep += shootingStep;
                else{
                    currStep -= shootingStep;
                    dec = true;
                    inc = false;
                }
            }
            else{
                currStep -= shootingStep;
                if(currStep < 0){
                    currStep = 0;
                    dec = false;
                    shooting = false;
                }
            }
        }
    }

    void draw(){
        glPushMatrix();
        //glTranslatef(tx, ty, tz);
        rotTail += alive * 5;
        if(rotTail >= 360) rotTail -= 360;
        glPushMatrix();
            drawBody();
        glPopMatrix();
        glPushMatrix();
            glTranslatef(-R-r, 0, 2 * H / 3);
            arm();
        glPopMatrix();

        glPushMatrix();
            glTranslatef(R + r, 0, 2 * H / 3);
            arm();
        glPopMatrix();

        glPopMatrix();
    }

    void shootIt(){
        if(shooting) return;

        shooting  = true;
        inc = true;
        currStep += shootingStep;
    }

    bool check(){
        myVector L = LOOK;
        L.makeUnitVec();
        point A = point(EYE.x, EYE.y, EYE.z);
        float projDist = L.projectionDist(A, top); // projection of top point
        point projPoint = A.moveIn({L.x, L.y, L.z}, projDist);
        float projRad = L.projectionRad(A, top);
        if(projRad <= oR) {
            sndPlaySound("lagse.wav", SND_ASYNC);
            health -= 50;
            if(health <= 0) alive = false;
            return true;
        }

        myVector ownUp = myVector(top.x - bottom.x, top.y - bottom.y, top.z - bottom.z);
        ownUp.makeUnitVec();

        float FinalProjDistFromBottom = ownUp.projectionDist(bottom, projPoint);

        if(FinalProjDistFromBottom < 0) {
            //sndPlaySound("shot.wav", SND_ASYNC);
            return false;
        }

        point FinalProj = bottom.moveIn({ownUp.x, ownUp.y, ownUp.z}, FinalProjDistFromBottom);

        float distInUp = FinalProjDistFromBottom;
        float FinalProjDist = ownUp.projectionRad(bottom, projPoint);

        if(distInUp > H + oR) {
            //sndPlaySound("shot.wav", SND_ASYNC);
            return false;
        }

        if(distInUp <= H){
            if(FinalProjDist <= R){
                sndPlaySound("lagse.wav", SND_ASYNC);
                //drawBlood(A);
                //cout << "Shooted" << endl;
                health -= 40;
                if(health <= 0) alive = false;
                return true;
            }
        }

//        float ex = H + oR - distInUp;
//        float expected = sqrt(oR * oR - ex * ex);
//
//        cout << FinalProjDist << " " << expected << endl;
//
//        if(FinalProjDist > expected) return false;
//
//        //cout << "SHooted" << endl;
        //sndPlaySound("shot.wav", SND_ASYNC);
        return false;
    }
};

vector<vector<bool>> stars;

void sky(double radius,int slices = 100,int stacks = 100){
    colorMaterial(getColor(8, 27, 89));

	struct point points[slices + 1][stacks + 1];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(M_PI/2));
		r=radius*cos(((double)i/(double)stacks)*(M_PI/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*M_PI);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*M_PI);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		//glColor3f(1, 1, 1);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);

			    getNormal3p(points[i][j].x,points[i][j].y,points[i][j].z,
                    points[i][j+1].x,points[i][j+1].y,points[i][j+1].z,
                    points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);

			    //upper hemisphere
			    glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
			    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
			    glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
			glEnd();
		}
	}

	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		//glColor3f(1, 1, 1);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_POINTS);
                colorMaterialEmit({1, 1, 1});
				if(stars[i][j]) glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
			glEnd();
		}
	}
}

vector<Monster> currM;
vector<float> rotProp;
vector<myVector> tarM;
vector<float> lx, rx, ly, ry;
int iM = 0;
bool gameOn = false;
bool P = false;

void drawHelperTar(){
    colorMaterial({1, 0, 0});

    point A = point(EYE.x, EYE.y, EYE.z).moveIn({LOOK.x, LOOK.y, LOOK.z}, 1);

    glPushMatrix();
    glBegin(GL_LINES);
        glVertex3f(A.x - 0.05 * RIGHT.x, A.y - 0.05 * RIGHT.y, A.z - 0.05 * RIGHT.z);
        glVertex3f(A.x + 0.05 * RIGHT.x, A.y + 0.05 * RIGHT.y, A.z + 0.05 * RIGHT.z);

        glVertex3f(A.x - 0.05 * UP.x, A.y - 0.05 * UP.y, A.z - 0.05 * UP.z);
        glVertex3f(A.x + 0.05 * UP.x, A.y + 0.05 * UP.y, A.z + 0.05 * UP.z);
    glEnd();
    glPopMatrix();
}

void display(void){
 glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    //glFrustum(-1,1,-1,1,0.8,50);
    gluPerspective(45, ar, 0.8, 100);
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    gluLookAt(EYE.x, EYE.y, EYE.z,   EYE.x + LOOK.x, EYE.y + LOOK.y, EYE.z + LOOK.z,   UP.x,UP.y,UP.z);
    light();

    drawHelperTar();

//Main Code
    play();
    glPushMatrix();
        glTranslatef(25, 25, 0);
        sky(60);
    glPopMatrix();

    glPushMatrix();
        drawCube(53, 0.5, 6, getColor(38, 39, 41));
        glTranslatef(0, 50 - 0.5, 0);
        drawCube(53, 0.5, 6, getColor(38, 39, 41));
    glPopMatrix();

    glPushMatrix();
        glTranslatef(0, 0.5, 0);
        drawCube(0.5, 49, 6, getColor(38, 39, 41));
        glTranslatef(53 - 0.5, 0, 0);
        drawCube(0.5, 49, 6, getColor(38, 39, 41));
    glPopMatrix();

if(P){
    glPushMatrix();
        for(int i = 0; i < currM.size(); i++){
            glPushMatrix();
                if(EYE.x >= lx[i] && EYE.x <= rx[i] && EYE.y >= ly[i] && EYE.y <= ry[i] && currM[i].alive){
                    myVector lk = myVector(currM[i].face.x - currM[i].top.x,
                                           currM[i].face.y - currM[i].top.y,
                                           0);
                    myVector ly = myVector(EYE.x - currM[i].face.x, EYE.y - currM[i].face.y, 0);
                    ly.makeUnitVec();
                    float ang = lk.angleWith(ly);

                    point newP(currM[i].face.x + currM[i].ms * ly.x, currM[i].face.y + currM[i].ms * ly.y, currM[i].face.z + currM[i].ms * ly.z);

                    if(newP.dist(point(EYE.x, EYE.y, EYE.z)) > 2.5){
                        currM[i].tx += currM[i].ms * ly.x;  currM[i].top.x += currM[i].ms * ly.x;
                        currM[i].ty += currM[i].ms * ly.y;  currM[i].top.y += currM[i].ms * ly.y;
                        currM[i].tz += currM[i].ms * ly.z;  currM[i].top.z += currM[i].ms * ly.z;

                        currM[i].bottom.x += currM[i].ms * ly.x;  currM[i].face.x += currM[i].ms * ly.x;
                        currM[i].bottom.y += currM[i].ms * ly.y;  currM[i].face.y += currM[i].ms * ly.y;
                        currM[i].bottom.z += currM[i].ms * ly.z;  currM[i].face.z += currM[i].ms * ly.z;
                    }

                    glTranslatef(currM[i].tx, currM[i].ty, currM[i].tz);
                    glRotatef(ang, 0, 0, 1);
                    currM[i].draw();
                    currM[i].shootIt();
                }
                else{
                    glTranslatef(currM[i].tx, currM[i].ty, currM[i].tz);
                    currM[i].draw();
                }

                if(shoot){
                    if(currM[i].check()){
                        shoot = false;
                    }
                }
            glPopMatrix();
        }
    glPopMatrix();

    if(shoot){
        sndPlaySound("shot.wav", SND_ASYNC);
        shoot = false;
    }

}


//MainCode End


    glFlush();
    glutSwapBuffers();
}

void movewitharrow( int key, int x, int y ){
	switch ( key )
	{
		case GLUT_KEY_DOWN:
		    EYE.Move(LOOK, -step);
			break;
		case GLUT_KEY_UP:
		    EYE.Move(LOOK, step);
			break;
		case GLUT_KEY_RIGHT:
			EYE.Move(RIGHT, step);
			break;
		case GLUT_KEY_LEFT:
			EYE.Move(RIGHT, -step);
			break;
		case GLUT_KEY_PAGE_UP:
			EYE.Move(UP, step);
			break;
		case GLUT_KEY_PAGE_DOWN:
			EYE.Move(UP, -step);
			break;

		case GLUT_KEY_END:	// Escape key
			exit(1);
	}

	glutPostRedisplay();
}

void myKeyboardFunc( unsigned char key, int x, int y ){
	switch ( key )
	{
        case 'a':
            LOOK.Rotate(UP, -stepAngle);
			RIGHT.Rotate(UP, -stepAngle);
            break;
        case 'd':
            LOOK.Rotate(UP, stepAngle);
			RIGHT.Rotate(UP, stepAngle);
            break;
        case 'S':
            shoot = !shoot;
            break;
        case 'w':
            LOOK.Rotate(RIGHT, -stepAngle);
			UP.Rotate(RIGHT, -stepAngle);
            break;
        case 'z':
            LOOK.Rotate(RIGHT, stepAngle);
			UP.Rotate(RIGHT, stepAngle);
            break;
        case 'P':
            P = !P;
            EYE.x = 50, EYE.y = 1, EYE.z = 4.5;
            break;
        case 'q':
            //UP.Rotate(LOOK, stepAngle);
			//RIGHT.Rotate(LOOK, stepAngle);
			break;
        case 'x':
            //UP.Rotate(LOOK, -stepAngle);
			//RIGHT.Rotate(LOOK, -stepAngle);
			break;
        case 'L':
            lamp = !lamp;
            break;
        case 'l':
            Mainlight = !Mainlight;
            break;
		case 27:	// Escape key
			exit(1);
	}
	glutPostRedisplay();
}

void LoadTexture(const char*filename, unsigned int& ID){
    glGenTextures(1, &ID);
    glBindTexture(GL_TEXTURE_2D, ID);
    glPixelStorei(GL_UNPACK_ALIGNMENT, ID);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    BmpLoader bl(filename);
    gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGB, bl.iWidth, bl.iHeight, GL_RGB, GL_UNSIGNED_BYTE, bl.textureData );
}

int main (int argc, char **argv){

    cout << "Press P to start" << endl;
    cout << "Use arrows to move backward forward left & right" << endl;
    cout << "a - d to rotate left and right" << endl;
    cout << "w - z to rotate up & down" << endl;
    cout << "S to shoot" << endl;

    currM = {
        Monster(0.5, 1.5, 0.2, 0.2, 42, 12, 3),
        Monster(0.5, 1.5, 0.2, 0.2, 50, 45, 3),

        Monster(0.5, 1.5, 0.2, 0.2, 25, 36, 3 + 3),

        Monster(0.5, 1.5, 0.2, 0.2, 24, 12, 3 + 3),
        Monster(0.5, 1.5, 0.2, 0.2, 32, 12, 3 + 3),

        Monster(0.5, 1.5, 0.2, 0.2, 7 + 15 - 2, 7 + 2, 3 + 3),

        Monster(0.5, 1.5, 0.2, 0.2, 9, 50 - 7 - 2, 3 + 3),
        Monster(0.5, 1.5, 0.2, 0.2, 12, 50 - 7 - 2, 3 + 3),
        Monster(0.5, 1.5, 0.2, 0.2, 15, 50 - 7 - 2, 3 + 3),
        Monster(0.5, 1.5, 0.2, 0.2, 18, 50 - 7 - 2, 3 + 3)
    };

    lx = {38,  38,  22,  22,  22,  7,  7,  7, 7, 7};
    rx = {53,  53,  40,  34,  34,  22, 22, 22, 22, 22};

    ly = {0,   0,  20,  7,  7, 7,  22, 22, 22, 22};
    ry = {18, 53,  40, 26, 26, 21, 40, 40, 40, 40};

    nCr.resize(101, vector<long long> (101));
    calnCr(100, 100);

    stars.resize(100, vector<bool> (100));
    for(int i = 0; i < 100; i++){
        for(int j = 0; j < 100; j++){
            int k = uniform_int_distribution<int>(1, 50)(rng);
            stars[i][j] = k == 5;
        }
    }

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

    glutInitWindowPosition(300,80);
    glutInitWindowSize(windowWidth, windowHeight);
    glutCreateWindow("Shoot 'Em All (ver 3)");


    glShadeModel( GL_SMOOTH );
    glEnable( GL_DEPTH_TEST );
    glEnable( GL_NORMALIZE );
    glEnable( GL_LIGHTING );

    ID.resize(2);
    ID[0] = 0;
    ID[1] = 1;
    LoadTexture("F:\\Shoot 'Em All V3\\wall.bmp", ID[0]);
    //LoadTexture("F:\\Shoot 'Em All V3\\slash.bmp", ID[1]);

//    auto finish = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed = finish - start;

    glutReshapeFunc(resize);
    glutDisplayFunc(display);
    glutKeyboardFunc(myKeyboardFunc);
    glutSpecialFunc(movewitharrow);
    glutIdleFunc(animate);

    glutMainLoop();

    return 0;
}
