/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include "light_source.h"

void PointLight::shade( Ray3D& ray ) {

	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  
	
	// scene signature
	/* ray.col = ray.intersection.mat->diffuse; */

	//scene with diffuse plus ambient shading
	Vector3D light_vector = get_position() - ray.intersection.point;
	light_vector.normalize();

	// Dot product
	double n_dot_l = ray.intersection.normal.dot(light_vector);
	if (n_dot_l < 0)
		n_dot_l = 0;

	// Colour colour = _col_ambient * ray.intersection.mat->ambient
	// 	+ n_dot_l * ray.intersection.mat-> diffuse * _col_diffuse;

	// colour.clamp();
	// ray.col = colour;

	// scene with the full Phong model
	Vector3D view_vector = - ray.dir;
	view_vector.normalize();

	double specular = (2 * ray.intersection.normal.dot(light_vector) * ray.intersection.normal 
					 - light_vector).dot(view_vector);

	if (specular < 0)
		specular = 0;

	specular = pow(specular, ray.intersection.mat->specular_exp);

	Colour colour = _col_ambient * ray.intersection.mat->ambient
		+ n_dot_l * ray.intersection.mat-> diffuse * _col_diffuse
		+ specular * ray.intersection.mat->specular * _col_specular;

	colour.clamp();
	ray.col = colour;
}

void AreaLight::shade(Ray3D &ray) {
   // to realise soft shadows
   

    Vector3D N = ray.intersection.normal;
    Vector3D L = this->get_position() - ray.intersection.point;
    Vector3D V = -ray.dir;
    Vector3D R = 2.* (L.dot(N) * N) - L;

    N.normalize();
    L.normalize();
    V.normalize();
    R.normalize();

    Colour Ia = _col_ambient;
    Colour Id = _col_diffuse;
    Colour Is = _col_specular;
    Colour Ka = ray.intersection.mat->ambient;
    Colour Kd = ray.intersection.mat->diffuse;
    Colour Ks = ray.intersection.mat->specular;

    Colour ambient = (1) ? (Ka * Ia) : Colour(0,0,0);
    Colour diffuse = (1) ? (fmax(0, N.dot(L)) * Kd * Id) : Colour(0,0,0);
    Colour specular = (1) ? (pow(fmax(0, V.dot(R)), ray.intersection.mat->specular_exp) * Ks * Is) : Colour(0,0,0);

    if(ray.in_shadow){
		ray.col = ambient;
		return;
	}
	
    ray.col = ambient + diffuse + specular;
}
