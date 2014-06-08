/*----------------------------------------------------------
* COSC363  Ray Tracer
*
*  The Cylinder class
*  This is a subclass of Object, and hence implements the
*  methods intersect() and normal().
-------------------------------------------------------------*/

#include "Cylinder.h"
#include <math.h>

/**
* Cylinder's intersection method.  The input is a ray (pos, dir). 
*/
float Cylinder::intersect(Vector pos, Vector dir)
{
	
    float a = (dir.x*dir.x + dir.z*dir.z);
    float b = 2 * (dir.x * (pos.x-base.x) + dir.z * (pos.z - base.z));
    float c = pow((pos.x - base.x),2.0) + pow((pos.z - base.z), 2.0) - pow(radius, 2);
    float delta = pow(b,2) - 4*a*c;
   
	if(fabs(delta) < 0.001) return -1.0; 
    if(delta < 0.0) return -1.0;

    float t1 = (-b - sqrt(delta)) / (2 * a);
    float t2 = (-b + sqrt(delta)) / (2 * a);
    
   
    if(fabs(t1) < 0.001 )
    {
        if (t2 > 0) return t2;
        else t1 = -1.0;
    }

    if(fabs(t2) < 0.001 ) t2 = -1.0;
	
	float t = (t1 < t2) ? t1 : t2;
	
	float intersectionHeight = pos.y + (dir.y*t);
	float insideCoeff = (intersectionHeight - base.y);
	if (insideCoeff >= base.y && insideCoeff <= height)
	{
		return t;
	}
	
	// Doesn't intersect
	return -1;
}

/**
* Returns the unit normal vector at a given point.
* Assumption: The input point p lies on the sphere.
*/
Vector Cylinder::normal(Vector p)
{
    Vector n = Vector(p.x - base.x, 0, p.z - base.z);
    n.normalise();
    return n;
}
