// ========================================================================
// COSC 363  Computer Graphics  Lab07
// A simple ray tracer
// The Geany compile commands:
// g++ -Wall -c "%f" Color.cpp CheckedPlane.cpp Object.cpp Plane.cpp Sphere.cpp Vector.cpp
// ========================================================================

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include "Vector.h"
#include "Sphere.h"
#include "Color.h"
#include "Object.h"
#include "Plane.h"
#include "CheckedPlane.h"
#include "Cylinder.h"
#include <GL/glut.h>
#include <omp.h>
using namespace std;


#define PIXELS_HORZ 600
#define PIXELS_VERT 600
#define PI 3.14159265358979

// Antialiasing. Note that if n > 1, it will take n*n times as long as 1
#define ANTI_ALIAS_FACTOR 3

// These cannot be applied together; invertColor takes priority if true.
const bool invertColor = false;
const bool monochrone = false;
// Change this to the color of the light you want to use.
// Eg: Color(1, 0.5, 0) will take 100% of the red channel and 50% of
// of the green channel and 0% of the blue channel.
// This above example will show only the red and green channels.
const Color monoChromeColor = Color(1, 0.5, 0);

const float WIDTH = 20.0;  
const float HEIGHT = 20.0;
const float EDIST = 40.0;
const int PPU = 30;
const float PIXEL_SIZE = 1.0/PPU;
const int MAX_STEPS = 5;
const float XMIN = -WIDTH * 0.5;
const float XMAX =  WIDTH * 0.5;
const float YMIN = -HEIGHT * 0.5;
const float YMAX =  HEIGHT * 0.5;

// Progress used within the progress bar
static int progress = 0;

GLuint txId[1];

// Rotation matrix for PI/4 radians (45 degrees).
float rotationMatrix[4][4] = {{cos(PI/4), 0, -sin(PI/4), 0},
							   {0, 1, 0, 0},
							   {sin(PI/4), 0, cos(PI/4), 0},
							   {0, 0, 0, 1}};

vector<Object*> sceneObjects;

// Light array
Vector light[3] = {Vector(10.0, 40.0, -5.0),
                   Vector(0, 20.0, -35.0), 
                   Vector(-10, -5, 0.0)};
float lightAngle[3] = {360, PI/50, 360};
Vector lightAim[3] = {Vector(0, 0, 0), Vector(0,-10,-35.0), Vector(0, 0, -50)};
Color lightColor[3] = {Color::WHITE, Color(0, 1, 0), Color::WHITE};

Color shadowColor = Color(0, 0, 0);
Color fogColor = Color(0, 0, 0);

//A useful struct
struct PointBundle
{
	Vector point;
	int index;
	float dist;
};

/**
* This function compares the given ray with all objects in the scene
* and computes the closest point  of intersection.
*/
PointBundle closestPt(Vector pos, Vector dir)
{
    Vector  point(0, 0, 0);
	float min = 10000.0;

	PointBundle out = {point, -1, 0.0};

    for(unsigned int i = 0;  i < sceneObjects.size();  i++)
	{
        float t = sceneObjects[i]->intersect(pos, dir);
		if(t > 0)        //Intersects the object
		{
			point = pos + dir*t;
			if(t < min)
			{
				out.point = point;
				out.index = i;
				out.dist = t;
				min = t;
			}
		}
	}

	return out;
}

/**
 * Computes the colour value obtained by tracing a ray.
 * If reflections and refractions are to be included, then secondary rays will 
 * have to be traced from the point, by converting this method to a recursive
 * procedure.
 */
Color trace(Vector pos, Vector dir, int step, int lightIndex)
{
	Color colorSum;
	Color temp;
	
	Vector lightPosition = light[lightIndex];
	bool isSpotight = lightAngle[lightIndex] <= 180;
	float spotlightAngle = lightAngle[lightIndex];
	Vector lightNormal = lightPosition - lightAim[lightIndex];// - lightPosition;
	Color lightCol = lightColor[lightIndex];
	
    PointBundle q = closestPt(pos, dir);
   
    q.dist = (q.point - pos).length();
    
	Vector n;
    if (q.index != 0) {
		if(q.index == -1) {
			// Intersects with no objects in the scene
			colorSum = fogColor;
			q.dist = 200; // Maximum fog.
		}
		else {
			n = sceneObjects[q.index]->normal(q.point);
			Vector l = light[lightIndex] - q.point;
			l.normalise();
			float lDotn = l.dot(n);
			
			PointBundle shadowRay = closestPt(q.point, l);

			Color col = sceneObjects[q.index]->getColor(q.point); //Object's colour
			
			float reflectCoeff = 0.8;
			
			Vector r = ((n * 2) * lDotn) - l;
			r.normalise();
			Vector v(-dir.x, -dir.y, -dir.z);
			
			//View vector;
			float rDotv = r.dot(v);
			float spec;
			if(rDotv < 0) spec = 0.0;
			else spec = pow(rDotv, 10);
			
			colorSum = col.phongLight(shadowColor, lDotn, spec);


			if 	(lDotn <= 0) {
				// In the 'dark side'; removes the sharp cresent-shape on
				// the spheres
				temp = shadowColor;
				temp.combineColor(colorSum, 0.1);
				colorSum.combineColor(shadowColor, 0.1);
			}
			
			if (shadowRay.index != -1) {
				// In shadow
				colorSum = col.phongLight(shadowColor, 0.0, 0.0);
				Color softShadow = col.phongLight(shadowColor, lDotn, spec);
				colorSum.combineColor(softShadow, 0.5);
			}
			else  {
				// If not spotlight
				if (!isSpotight) {
					colorSum = col.phongLight(lightCol, lDotn, spec);
				}
				// Else if spotlight
				else if (acos(((lightPosition.dot(l))/(l.length()*lightNormal.length())) > spotlightAngle)) {
					colorSum = col.phongLight(lightCol, 0.0, 0.0);
				}
			}

			
			//Generate reflection ray
			if((col.r + col.b + col.g < 2.9) && step < MAX_STEPS) {
				float nDotv = n.dot(v);
				Vector reflectionVector = ((n*2)*(nDotv) - v);
				reflectionVector.normalise();
				Color reflectionCol = trace(q.point, reflectionVector, step+1, lightIndex);
				colorSum.combineColor(reflectionCol, reflectCoeff);
			}
			
		}
	}
	else {
		//Refraction
		if (step <  MAX_STEPS) //refractive sphere
		{
			float n1, n2;
			n = sceneObjects[q.index]->normal(q.point);
			//center and radius of refractive sphere
			Vector c = Vector(5, -1, -50);
			float r = 4.0;
			
			float dist = ((pos-c).length()) - r;
			if (fabs(dist) < 0.001) //ray leaving sphere
			{
				n1 = 1.004;
				n2 = 1.0;
				n.scale(-1.0);
			} else { //ray entering sphere
				n1 = 1.0;
				n2 = 1.004;
			}
			
			float cosAngle = sqrt(1 - (pow((n1/n2), 2) * (1 - pow((dir.dot(n)), 2))));
			Vector refractionVector = (dir * (n1/n2)) - (n * (((n1/n2) * (dir.dot(n))) + cosAngle));
			
			Color refractionCol = trace(q.point, refractionVector, step+1, lightIndex);
			colorSum = refractionCol;
		}
	}
	
	
	// Fog
	if (q.dist >= 100) {
		colorSum = fogColor;
	}
	else if (step == 1) {
		temp = fogColor;
		float fogCoeff = (100 - q.dist)/100.0f;
		temp.combineColor(colorSum, fogCoeff);
		colorSum = temp;
	}
	
    return colorSum;

}


// AA factor must be a multiple of 4!!!!!
Color runAntiAliasing(Vector eye, float x1, float y1, int aaFactor)
{
	if ((aaFactor % 4) != 0)
	{
		return Color::RED; // undesirable color if paramaters were bad.
	}	
	
	float redPart = 0;
	float greenPart = 0;
	float bluePart = 0;

	int X__BIT = 1;
	int Y__BIT = 2;

	int innerAAFactor = aaFactor / 4;
	float fractionOfPixel = PIXEL_SIZE / float(innerAAFactor);

	int i = 0;
	for (i = 0; i < 4; i ++)
	{
		int j = 1;
		for (j = 1; j <= innerAAFactor; j++)
		{
			Vector direction(x1 + (i & X__BIT ? (j * fractionOfPixel) : 0),
						 y1 + (i & Y__BIT ? j * (fractionOfPixel) : 0),
						 -EDIST);
			direction.normalise();
			Color col = trace (eye, direction, 1, 0);
			Color col2 = trace(eye, direction, 1, 1);
			Color col3 = trace(eye, direction, 1, 2);
			if(col.r > 1 || col.g > 1 || col.b > 1) {
				float colMax = col.r > col.g ? col.b > col.r ?
					col.b : col.r : col.b > col.g ? col.b : col.g;
				col.r = col.r/colMax;
				col.g = col.g/colMax;
				col.b = col.b/colMax;
			}
			if(col2.r > 1 || col2.g > 1 || col2.b > 1) {
				float col2Max = col2.r > col2.g ? col2.b > col2.r ?
					col2.b : col2.r : col2.b > col2.g ? col2.b : col2.g;
				col2.r = col2.r/col2Max;
				col2.g = col2.g/col2Max;
				col2.b = col2.b/col2Max;
			}
			if(col3.r > 1 || col3.g > 1 || col3.b > 1) {
				float col3Max = col3.r > col3.g ? col3.b > col3.r ?
					col3.b : col3.r : col3.b > col3.g ? col3.b : col3.g;
				col3.r = col3.r/col3Max;
				col3.g = col3.g/col3Max;
				col3.b = col3.b/col3Max;
			}
			redPart += (col.r + col2.r + col3.r)/3;
			greenPart += (col.g + col2.g+col3.g)/3;
			bluePart += (col.b + col2.b+col3.b)/3;
		}
		
	}

	redPart /= aaFactor;
	bluePart /= aaFactor;
	greenPart /= aaFactor;

	Color anti_aliasedColor = Color(redPart, greenPart, bluePart);

	return anti_aliasedColor;
}

void progressBar(int percent) {
	if (percent > progress) {
		cout << "\r" << setw(3) << percent << "%  [";
		cout << string(percent/3, '-') << "(━┳━ _ ━┳━)" << string(100/3-percent/3, ' ') << "]" << flush;

		progress = percent;
	}
}

//---The main display module -----------------------------------------------------------
// In a ray tracing application, it just displays the ray traced image by drawing
// each pixel as quads.
//---------------------------------------------------------------------------------------
void display()
{
	int widthInPixels = (int)(WIDTH * PPU);
	int heightInPixels = (int)(HEIGHT * PPU);
	float x1, y1;
	Vector eye(0., 0., 0.);
	
	float screenBuffer[widthInPixels][heightInPixels][3];
	
	for(int i = 0; i < widthInPixels; i++)	//Scan every "pixel"
	{
		x1 = XMIN + i*PIXEL_SIZE;
		progressBar(100 * (i + 1) / widthInPixels);
		
		// Used for enabling multi-threading.
		#pragma omp parallel for
		for(int j = 0; j < heightInPixels; j++)
		{
			y1 = YMIN + j*PIXEL_SIZE;
			
			Color anti_aliasedColor = runAntiAliasing(eye, x1, y1, ANTI_ALIAS_FACTOR * 4);
			screenBuffer[i][j][0] = anti_aliasedColor.r;
			screenBuffer[i][j][1] = anti_aliasedColor.g;
			screenBuffer[i][j][2] = anti_aliasedColor.b;
        }
    }
    
	glClear(GL_COLOR_BUFFER_BIT);
	glBegin(GL_QUADS);  //Each pixel is a quad.
	for(int i = 0; i < widthInPixels; i++)	//Scan every "pixel"
	{
		x1 = XMIN + i*PIXEL_SIZE;
		for(int j = 0; j < heightInPixels; j++) {
			y1 = YMIN + j*PIXEL_SIZE;
			if (invertColor) {
				glColor3f(1-screenBuffer[i][j][0],
			          1-screenBuffer[i][j][1],
			          1-screenBuffer[i][j][2]);
		    }
			else if (monochrone) {
				
				glColor3f(screenBuffer[i][j][0]*monoChromeColor.r,
						  screenBuffer[i][j][1]*monoChromeColor.g,
						  screenBuffer[i][j][2]*monoChromeColor.b);
			}
			else {
			glColor3f(screenBuffer[i][j][0],
			          screenBuffer[i][j][1],
			          screenBuffer[i][j][2]);
			}
			glVertex2f(x1, y1);				//Draw each pixel with its color value
			glVertex2f(x1 + PIXEL_SIZE, y1);
			glVertex2f(x1 + PIXEL_SIZE, y1 + PIXEL_SIZE);
			glVertex2f(x1, y1 + PIXEL_SIZE);
		}
	}
    glEnd();
    glFlush();
}

/**
 * Adds spheres to the scene.
 */
void makeBalls() {
    Sphere *sphere1 = new Sphere(Vector(4, 8, -60), 1.5, Color(1, 1, 0));
	sceneObjects.push_back(sphere1);
	Sphere *sphere2 = new Sphere(Vector(0, 6, -75), 8.0, Color(0.9,0.9,0.9));
	sceneObjects.push_back(sphere2);
	Sphere *sphere3 = new Sphere(Vector(0, -5, -65), 2.0, Color::GREEN);
	sceneObjects.push_back(sphere3);
	Sphere *sphere4 = new Sphere(Vector(-3, -5, -50), 3.0, Color::RED);
	sceneObjects.push_back(sphere4);
	Cylinder *cylinder = new Cylinder(Vector(-7, -5, -55), 6.0, 1.0, Color::GREEN);
	sceneObjects.push_back(cylinder);
}


void drawParallelpiped(Vector center, Vector xAxis, Vector yAxis, Vector zAxis, Vector lengths, Color color) {
	
	xAxis.normalise();
	yAxis.normalise();
	zAxis.normalise();
	
	Vector v000 = center;
	Vector v001 = center + zAxis*lengths.z;
	Vector v010 = center + yAxis*lengths.y;
	Vector v011 = v010 + zAxis*lengths.z;
	Vector v100 = center + xAxis*lengths.x;
	Vector v101 = v100 + zAxis*lengths.z;
	Vector v110 = v100 + yAxis*lengths.y;
	Vector v111 = v110 + zAxis*lengths.z;
	
	Plane *back = new Plane(v101, v001, v011, v111, color);
	sceneObjects.push_back(back);
	Plane *right = new Plane(v100, v101, v111, v110, color);
	sceneObjects.push_back(right);
	Plane *front = new Plane(v000, v100, v110, v010, color);
	sceneObjects.push_back(front);
	Plane *left = new Plane(v001, v000, v010, v011, color);
	sceneObjects.push_back(left);
	Plane *top = new Plane(v010, v110, v111, v011, color);
	sceneObjects.push_back(top);
	Plane *bottom = new Plane(v100, v000, v001, v101, color);
	sceneObjects.push_back(bottom);
	

}


/**
 * Adds a floor to the scene.
 */
void makeFloor() {
	int floorX = 300;
	int floorY = 300;
	int numTiles = 39;
	
	
	Vector frontLeftBottom = Vector(-floorX/2, -10, 0);
	Vector frontRightBottom = Vector(+floorX/2, -10, 0);
	Vector backRightBottom = Vector(+floorX/2, -10, -floorY/2);
	Vector backLeftBottom = Vector(-floorX/2, -10, -floorY/2);
	
	Vector backRightTop = Vector(+floorX/2, 20, -floorY/2);
	Vector backLeftTop = Vector(-floorX/2, 20, -floorY/2);
	
	CheckedPlane *floor = new CheckedPlane(frontLeftBottom,
							 frontRightBottom,
							 backRightBottom,
							 backLeftBottom,
							 Color::BLACK, Color::WHITE, numTiles);
	sceneObjects.push_back(floor);
	
	Vector bringForwardLeft = Vector(100, 0, 0);
	
	Plane *plane = new Plane(
							 backLeftBottom + bringForwardLeft,
							 backRightBottom,
							 backRightTop,
							 backLeftTop + bringForwardLeft,
							 Color::BLACK);
	sceneObjects.push_back(plane);
}

/**
 * Initialises the scene and populates it with spheres and stuff.
 */
void initialize() {
    glMatrixMode(GL_PROJECTION);
    gluOrtho2D(XMIN, XMAX, YMIN, YMAX);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    Sphere* refractiveSphere = new Sphere(Vector(5, -1, -50), 4, Color::WHITE);
	sceneObjects.push_back(refractiveSphere);
    
    
    makeBalls();
    
    Vector center = Vector(5, -5, -50);
    Vector xAxis = Vector(cos(PI/5), 0, sin(PI/5));
    Vector yAxis = Vector(0, -1, 0);
    Vector zAxis = Vector(-sin(PI/5), 0, cos(PI/5));
    Vector lengths = Vector(5, 5, 3);
    drawParallelpiped(center, xAxis, yAxis, zAxis, lengths, Color::BLUE);
    
    makeFloor();
    
	glClearColor(0, 0, 0, 1);
	
}


int main(int argc, char *argv[]) 
{
	omp_set_num_threads(1000);
	cout << "Starting" << endl;
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(20, 20);
    glutCreateWindow("Raytracing");

    glutDisplayFunc(display);
    initialize();

    glutMainLoop();
    return 0;
}
