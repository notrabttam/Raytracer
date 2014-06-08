/*----------------------------------------------------------
* COSC363  Ray Tracer
*
*  The Plane class
*  This is a subclass of Object, and hence implements the
*  methods intersect() and normal().
-------------------------------------------------------------*/

#include "Plane.h"
#include "CheckedPlane.h"
#include "Color.h"
#include "Vector.h"
#include <stdlib.h>
#include <math.h>


Color CheckedPlane::getColor(Vector point)
{
	int X = tileCount*(point.x - a.x)/(b.x-a.x);
	int Y = tileCount*(point.z - a.z)/(d.z-a.z);
	return (X+Y)%2 == 0 ? color1 : color2;
}
