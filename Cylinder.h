/*----------------------------------------------------------
* COSC363  Ray Tracer
*
*  The cylinder class
*  This is a subclass of Object, and hence implements the
*  methods intersect() and normal().
-------------------------------------------------------------*/

#ifndef H_CYLINDER
#define H_CYLINDER

#include "Object.h"

/**
 * Defines a simple Cylinder located at 'center' 
 * with the specified radius
 */
class Cylinder : public Object
{

private:
    Vector base;
    float radius;
    float height;

public:	
	Cylinder()
		: base(Vector()), radius(1), height(3)  //Default constructor creates a generic cylinder
	{
		color = Color::WHITE;
	};
	
    Cylinder(Vector c, float r, float h, Color col)
		: base(c), radius(r), height(h)
	{
		color = col;
	};

	float intersect(Vector pos, Vector dir);

	Vector normal(Vector p);

};

#endif //!H_CYLINDER
