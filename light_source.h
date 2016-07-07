/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		   light source classes

***********************************************************/

#include "util.h"
#include "stdlib.h"

// Base class for a light source.  You could define different types
// of lights here, but point light is sufficient for most scenes you
// might want to render.  Different light sources shade the ray 
// differently.
class LightSource {
public:
	virtual void shade( Ray3D& ) = 0;
	virtual Point3D get_position() const = 0; 
};

// A point light is defined by its position in world space and its
// colour.
class PointLight : public LightSource {
public:
	PointLight( Point3D pos, Colour col ) : _pos(pos), _col_ambient(col), 
	_col_diffuse(col), _col_specular(col) {}
	PointLight( Point3D pos, Colour ambient, Colour diffuse, Colour specular ) 
	: _pos(pos), _col_ambient(ambient), _col_diffuse(diffuse), 
	_col_specular(specular) {}
	void shade( Ray3D& ray );
	Point3D get_position() const { return _pos; }
	
private:
	Point3D _pos;
	Colour _col_ambient;
	Colour _col_diffuse; 
	Colour _col_specular; 
};

// A area light is defined by a origin and two range vectors in world 
// space and its colour
class AreaLight : public LightSource {
public:
    AreaLight(Point3D origin, Vector3D v1, Vector3D v2, Colour col)
        : _origin(origin), _v1(v1), _v2(v2), _n(v1.cross(v2)), _col_ambient(col),
            _col_diffuse(col), _col_specular(col) {}
    AreaLight(Point3D origin, Vector3D v1, Vector3D v2, Colour ambient, Colour diffuse, Colour specular)
        : _origin(origin), _v1(v1), _v2(v2), _n(v1.cross(v2)), _col_ambient(ambient), _col_diffuse(diffuse),
            _col_specular(specular) {}

    void shade(Ray3D &ray);
    // randomly return a postion in this area
    Point3D get_position() const {return _origin + (rand() % 100/ 100.0)* _v1 + ((rand() % 100/ 100.0)* _v2);}; 
    Vector3D get_v1() const { return _v1; }
    Vector3D get_v2() const { return _v2; }
    //Vector3D get_v3() const { return _v3; }
    Vector3D get_normal() const { return _n; }
private:
    Point3D _origin;
    Vector3D _v1;
    Vector3D _v2;
    //Vector3D _v3;
    Vector3D _n;
    Colour _col_ambient;
    Colour _col_diffuse;
    Colour _col_specular;
};
