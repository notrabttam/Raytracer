/*----------------------------------------------------------
* COSC363  Ray Tracer
*
*  The Plane class
*  This is a subclass of Object, and hence implements the
*  methods intersect() and normal().
-------------------------------------------------------------*/

#ifndef H_CHECKEDPLANE
#define H_CHECKEDPLANE

#include "Vector.h"
#include "Object.h"
#include "Plane.h"
#include "Color.h"

class CheckedPlane : public Plane
{
private:
	int tileCount;
	Color color1, color2;
	
public:	
	CheckedPlane(void);
	
    CheckedPlane(Vector pa, Vector pb, Vector pc, Vector pd, Color col1, Color col2, int numTiles)
		: color1(col1), color2(col2)
	{
		tileCount = numTiles;
		
		a = pa;
		b = pb;
		c = pc;
		d = pd;
	};
	
	virtual Color getColor(Vector point);
};

#endif //!H_CHECKEDPLANE
