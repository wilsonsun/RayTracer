/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"

bool UnitSquare::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSquare, which is
	// defined on the xy-plane, with vertices (0.5, 0.5, 0), 
	// (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.

	// Transform the ray into object space
	Ray3D rayOS;
	rayOS.origin = worldToModel * ray.origin;
	rayOS.dir = worldToModel * ray.dir;

	// t = (p1-q).n /(r.n).
	Vector3D normal(0.0, 0.0, 1.0);
	Point3D p1(0.5, 0.5, 0.0);
	Vector3D p1_minus_q = p1 - rayOS.origin;

	double t = ((p1 - rayOS.origin).dot(normal))/(rayOS.dir.dot(normal));

	// Check t value
	if (t < 0)
		return false;

	// This is used when intersecting multiply objects and only 
    // want to keep the nearest intersection.
    if (ray.intersection.none == false && t > ray.intersection.t_value)
    	return false;

	// Plug-in t to compute the intersection point
	rayOS.intersection.point = rayOS.origin + t * rayOS.dir;

	// Check if the intersection point fall within the plane
	if (rayOS.intersection.point[0] <= 0.5 && rayOS.intersection.point[0] >= -0.5
		&& rayOS.intersection.point[1] <= 0.5 && rayOS.intersection.point[1] >= -0.5) {

		// Convert intersection point to world space
		ray.intersection.point = modelToWorld * rayOS.intersection.point;
		// Convert normal vector to world space
		ray.intersection.normal = transNorm(worldToModel, normal);
		ray.intersection.normal.normalize();
		ray.intersection.none = false;
		ray.intersection.t_value = t;

		return true;
	}

	return false;
}

bool UnitSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSphere, which is centred 
	// on the origin.  
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.
	
	// Transform the ray into object space
	Point3D center(0.0, 0.0, 0.0);
	double r = 1.0;

	Ray3D rayOS;
	rayOS.origin = worldToModel * ray.origin;
	rayOS.dir = worldToModel * ray.dir;

	double a = rayOS.dir.dot(rayOS.dir);
	double b = 2 * (rayOS.dir.dot(rayOS.origin - center));
	double c = (rayOS.origin - center).dot(rayOS.origin - center) - r * r;
	double t;

	// t^2 = b^2 - 4ac
	double alpha = b * b - 4 * a * c;

	if (alpha < 0)
		return false;

	double t1 = (-b + sqrt(alpha))/(2 * a);
	double t2 = (-b - sqrt(alpha))/(2 * a);

	if (t1 > 0 && t2 > 0 && t1 < t2)
		t = t1;
	else if (t2 > 0)
		t = t2;
	else if (t1 > 0)
		t = t1;
	else
		return false;

	if (ray.intersection.none == false && t > ray.intersection.t_value)
		return false;

	rayOS.intersection.point = rayOS.origin + t * rayOS.dir;

	// Convert intersection point to world space
	ray.intersection.point = modelToWorld * rayOS.intersection.point;

	Vector3D normal = rayOS.intersection.point - center;

	// Convert normal vector to world space
	ray.intersection.normal = transNorm(worldToModel, normal);
	ray.intersection.normal.normalize();
	ray.intersection.none = false;
	ray.intersection.t_value = t;
	
	return true;
}

bool UnitCylinder::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {

	// Transform the ray into object space
	Ray3D rayOS;
	rayOS.origin = worldToModel * ray.origin;
	rayOS.dir = worldToModel * ray.dir;

	// Define radius of the cylinder
	double r = 1.0;

	// Calculaion of intersections between line and cylinder

	double a = pow(rayOS.dir[0], 2) + pow(rayOS.dir[2], 2);
	double b = 2 * (rayOS.origin[0] * rayOS.dir[0] + rayOS.origin[2] * rayOS.dir[2]);
	double c = pow(rayOS.origin[0], 2) + pow(rayOS.origin[2], 2) - pow(r, 2);

	// Alpha calculation
	double alpha = b * b - 4 * a * c;

	if (alpha < 0)
		return false;

	double t1 = (-b + sqrt(alpha))/(2 * a);
	double t2 = (-b - sqrt(alpha))/(2 * a);

	// Check if t1 and t2 fall with the height
	Point3D p1 = rayOS.origin + t1 * rayOS.dir;
	Point3D p2 = rayOS.origin + t2 * rayOS.dir;
	if (p1[1] > 1 || p1[1] < -1) {
		// If t1 doesn't fall with the height
		// Set it to 0
		t1 = -1;
	}
        
    if (p2[1] > 1 || p2[1] < -1) {
    	// If t2 doesn't fall with the height
		// Set it to 0
        t2 = -1;
    }

    // Find the intersection between the line and the plane
    double t3 = (Vector3D(0, 1, 0).dot(Vector3D(0, 1, 0)) - rayOS.origin[1]) 
    			/ rayOS.dir.dot(Vector3D(0, 1, 0));
    Point3D p3 = rayOS.origin + t3 * rayOS.dir;

    double t4 = (Vector3D(0, -1, 0).dot(Vector3D(0, -1, 0)) + rayOS.origin[1])
    			 / rayOS.dir.dot(Vector3D(0, -1, 0));
    Point3D p4 = rayOS.origin + t4 * rayOS.dir;

    if (pow(p3[0], 2) + pow(p3[2], 2) > 1) {
    	// If t3 doesn't fall with the height
		// Set it to 0
        t3 = -1;
    }

    if (pow(p4[0], 2) + pow(p4[2], 2) > 1) {
    	// If t3 doesn't fall with the height
    	// Set it to 0
         t4 = -1;
    }

    if (t1 < 0 && t2 < 0 && t3 < 0 && t4 < 0) {
    	// If No intersection
    	return false;
    }

    // Find the closest t
    double res;
    res = t1;
    if (res < 0 || (t2 < res && t2 > 0))
          res = t2;
    if (res < 0 || (t3 < res && t3 > 0))
          res = t3;
    if (res < 0 || (t4 < res && t4 > 0))
          res = t4;

    if (!ray.intersection.none && ray.intersection.t_value < res)
    	return false;

    rayOS.intersection.point = rayOS.origin + res * rayOS.dir;

    // Calculate the normal for the intersection
    Vector3D normal;
    if (res == t3)
    	normal = Vector3D(0, 1, 0);
    else if (res == t4)
    	normal = Vector3D(0, -1, 0);
    else
    	normal = rayOS.intersection.point - Point3D(0, rayOS.intersection.point[1], 0);

    ray.intersection.normal = transNorm(worldToModel, normal);
    ray.intersection.normal.normalize();
    ray.intersection.point = modelToWorld * rayOS.intersection.point;
    ray.intersection.t_value = res;
    ray.intersection.none = false;

    return true;

}

