/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		Implementations of functions in raytracer.h, 
		and the main function which specifies the 
		scene to be rendered..

***********************************************************/


#include "raytracer.h"
#include "bmp_io.h"
#include <cmath>
#include <iostream>
#include <cstdlib>

bool ANTI_ALISING_ENABLED = false;
bool GLOSSY_REFLECTION_ENABLED = false;
bool MOTION_BLUR_ENABLED = false;

Raytracer::Raytracer() : _lightSource(NULL) {
	_root = new SceneDagNode();
}

Raytracer::~Raytracer() {
	delete _root;
}

SceneDagNode* Raytracer::addObject( SceneDagNode* parent, 
		SceneObject* obj, Material* mat ) {
	SceneDagNode* node = new SceneDagNode( obj, mat );
	node->parent = parent;
	node->next = NULL;
	node->child = NULL;
	
	// Add the object to the parent's child list, this means
	// whatever transformation applied to the parent will also
	// be applied to the child..
	if (parent->child == NULL) {
		parent->child = node;
	}
	else {
		parent = parent->child;
		while (parent->next != NULL) {
			parent = parent->next;
		}
		parent->next = node;
	}
	
	return node;;
}

LightListNode* Raytracer::addLightSource( LightSource* light ) {
	LightListNode* tmp = _lightSource;
	_lightSource = new LightListNode( light, tmp );
	return _lightSource;
}

void Raytracer::rotate( SceneDagNode* node, char axis, double angle ) {
	Matrix4x4 rotation;
	double toRadian = 2*M_PI/360.0;
	int i;
	
	for (i = 0; i < 2; i++) {
		switch(axis) {
			case 'x':
				rotation[0][0] = 1;
				rotation[1][1] = cos(angle*toRadian);
				rotation[1][2] = -sin(angle*toRadian);
				rotation[2][1] = sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'y':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][2] = sin(angle*toRadian);
				rotation[1][1] = 1;
				rotation[2][0] = -sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'z':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][1] = -sin(angle*toRadian);
				rotation[1][0] = sin(angle*toRadian);
				rotation[1][1] = cos(angle*toRadian);
				rotation[2][2] = 1;
				rotation[3][3] = 1;
			break;
		}
		if (i == 0) {
		    node->trans = node->trans*rotation; 	
			angle = -angle;
		} 
		else {
			node->invtrans = rotation*node->invtrans; 
		}	
	}
}

void Raytracer::translate( SceneDagNode* node, Vector3D trans ) {
	Matrix4x4 translation;
	
	translation[0][3] = trans[0];
	translation[1][3] = trans[1];
	translation[2][3] = trans[2];
	node->trans = node->trans*translation; 	
	translation[0][3] = -trans[0];
	translation[1][3] = -trans[1];
	translation[2][3] = -trans[2];
	node->invtrans = translation*node->invtrans; 
}

void Raytracer::scale( SceneDagNode* node, Point3D origin, double factor[3] ) {
	Matrix4x4 scale;
	
	scale[0][0] = factor[0];
	scale[0][3] = origin[0] - factor[0] * origin[0];
	scale[1][1] = factor[1];
	scale[1][3] = origin[1] - factor[1] * origin[1];
	scale[2][2] = factor[2];
	scale[2][3] = origin[2] - factor[2] * origin[2];
	node->trans = node->trans*scale; 	
	scale[0][0] = 1/factor[0];
	scale[0][3] = origin[0] - 1/factor[0] * origin[0];
	scale[1][1] = 1/factor[1];
	scale[1][3] = origin[1] - 1/factor[1] * origin[1];
	scale[2][2] = 1/factor[2];
	scale[2][3] = origin[2] - 1/factor[2] * origin[2];
	node->invtrans = scale*node->invtrans; 
}

Matrix4x4 Raytracer::initInvViewMatrix( Point3D eye, Vector3D view, 
		Vector3D up ) {
	Matrix4x4 mat; 
	Vector3D w;
	view.normalize();
	up = up - up.dot(view)*view;
	up.normalize();
	w = view.cross(up);

	mat[0][0] = w[0];
	mat[1][0] = w[1];
	mat[2][0] = w[2];
	mat[0][1] = up[0];
	mat[1][1] = up[1];
	mat[2][1] = up[2];
	mat[0][2] = -view[0];
	mat[1][2] = -view[1];
	mat[2][2] = -view[2];
	mat[0][3] = eye[0];
	mat[1][3] = eye[1];
	mat[2][3] = eye[2];

	return mat; 
}

void Raytracer::traverseScene( SceneDagNode* node, Ray3D& ray ) {
	SceneDagNode *childPtr;

	// Applies transformation of the current node to the global
	// transformation matrices.
	_modelToWorld = _modelToWorld*node->trans;
	_worldToModel = node->invtrans*_worldToModel; 
	if (node->obj) {
		// Perform intersection.
		if (node->obj->intersect(ray, _worldToModel, _modelToWorld)) {
			ray.intersection.mat = node->mat;
		}
	}
	// Traverse the children.
	childPtr = node->child;
	while (childPtr != NULL) {
		traverseScene(childPtr, ray);
		childPtr = childPtr->next;
	}

	// Removes transformation of the current node from the global
	// transformation matrices.
	_worldToModel = node->trans*_worldToModel;
	_modelToWorld = _modelToWorld*node->invtrans;
}

void Raytracer::computeShading( Ray3D& ray ) {
   // Implement softshadow here
   LightListNode* curLight = _lightSource;
	for (curLight; curLight != NULL; curLight = curLight -> next) {
		if (curLight == NULL) break;
		// Each lightSource provides its own shading function.
         Colour previous_colour(ray.col);
         Vector3D shadowd;
         Ray3D shadow;
         int i;
		// Implement shadows here if needed.
         for(i = 0; i < 5; i++){
            shadowd = curLight->light->get_position() - ray.intersection.point;
            shadowd.normalize();
            shadow.origin = ray.intersection.point + 0.01 * shadowd;
            shadow.dir = shadowd;

            traverseScene(_root, shadow);
            
            ray.in_shadow = (!shadow.intersection.none) && (shadow.intersection.point - ray.intersection.point).length() < shadowd.length();
			
			
            curLight->light->shade(ray);
            
               
	}
	if(i > 0 && !shadow.intersection.none){
		ray.col = (1.0 / 5) * ray.col;
}
    ray.col = ray.col + previous_colour;
    ray.col.clamp();
}
    }

void Raytracer::initPixelBuffer() {
	int numbytes = _scrWidth * _scrHeight * sizeof(unsigned char);
	_rbuffer = new unsigned char[numbytes];
	_gbuffer = new unsigned char[numbytes];
	_bbuffer = new unsigned char[numbytes];
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			_rbuffer[i*_scrWidth+j] = 0;
			_gbuffer[i*_scrWidth+j] = 0;
			_bbuffer[i*_scrWidth+j] = 0;
		}
	}
}

void Raytracer::flushPixelBuffer( char *file_name ) {
	bmp_write( file_name, _scrWidth, _scrHeight, _rbuffer, _gbuffer, _bbuffer );
	delete _rbuffer;
	delete _gbuffer;
	delete _bbuffer;
}


Colour Raytracer::shadeRay( Ray3D& ray ) {
	Colour col(0.0, 0.0, 0.0); 
	traverseScene(_root, ray); 
	
	// Don't bother shading if the ray didn't hit 
	// anything.
	if (!ray.intersection.none) {
		computeShading(ray); 
		col = ray.col;  
		// You'll want to call shadeRay recursively (with a different ray, 
		// of course) here to implement reflection/refraction effects.  

		// Reflection effect
		//if (ray.intersection.mat->reflectivity > 0.0 && ray.depth < 5) {
			Ray3D rayFL;

            rayFL.depth = ray.depth + 1;
            // Reflective ray's origin will be intersection point
            ray.dir.normalize();

            // Calculate the reflective ray's direction
            rayFL.dir = ray.dir - 2 * ray.intersection.normal.dot(ray.dir) * ray.intersection.normal;
            rayFL.dir.normalize();
            rayFL.origin = ray.intersection.point + 0.01 * rayFL.dir;

            Colour fl_col;
            int blur_factor = 16;
            int num_of_samples = 50;

            if (GLOSSY_REFLECTION_ENABLED) {
             	Colour tc(0.0, 0.0, 0.0);

            	for (int i = 0; i < num_of_samples; i++) {
                	Ray3D rayT;
                	double t, p, x, y;

                  	t = acos(pow((1 - (rand() / (double) RAND_MAX)), ray.intersection.mat->reflectivity));
                   	p = 2 * M_PI * (rand() / (double) RAND_MAX);
                    x = sin(p) * cos(t) / blur_factor;
                  	y = sin(p) * sin(t) / blur_factor;

                 	Vector3D v1 = rayFL.dir.cross(ray.intersection.normal);
                   	Vector3D v2 = rayFL.dir.cross(v1);

                   	rayT.dir = x * v1 + y * v2 + rayFL.dir;
                   	rayT.dir.normalize();
                	rayT.origin = rayFL.origin;
                 	rayT.depth = rayFL.depth;
               		tc = tc + ray.intersection.mat->reflectivity *shadeRay(rayT);
   				}
				col = col + tc / num_of_samples;
            } else {

	            fl_col = shadeRay(rayFL);
	            //col = col + (ray.intersection.mat->reflectivity) * fl_col;
	            col = col + 1.0 * fl_col;
	        }

            col.clamp();
		//}

		// Refraction effect
		if (ray.intersection.mat->refractiveIndex > 1.0 && ray.depth < 10) {

			Vector3D rayOSDir;
			Vector3D refraDir;
			bool TIR = false;	// Total internal refraction
			Colour tc;
			double oddotn;

			rayOSDir = ray.dir - (2 * (ray.dir.dot(ray.intersection.normal)) * ray.intersection.normal);
			rayOSDir.normalize();
			

			Ray3D rayFL = Ray3D(ray.intersection.point + 0.001 * rayOSDir, rayOSDir);
			double ddotn = ray.dir.dot(ray.intersection.normal);

			// If the ray is coming in to the medium
			if (ddotn < 0) {

				// Calculate the refractive dir
				double temp = sqrt(1.0 - ((pow(ray.intersection.mat->refractiveIndex, 2) 
								* (1 - pow(ddotn, 2))) / pow(1.0, 2)));

				refraDir = (ray.intersection.mat->refractiveIndex * 
							(ray.dir - (ddotn * ray.intersection.normal))) - temp * ray.intersection.normal;
				refraDir.normalize();
				oddotn = -ray.dir.dot(ray.intersection.normal);
				tc = Colour(1.0, 1.0, 1.0);
			} else {
				// The ray is going out from the medium
				tc = Colour(exp(-0.15 * ray.intersection.mat->diffuse[0] * ray.intersection.t_value),
							exp(-0.15 * ray.intersection.mat->diffuse[1] * ray.intersection.t_value),
							exp(-0.15 * ray.intersection.mat->diffuse[2] * ray.intersection.t_value));

				double mddotn = ray.dir.dot(-ray.intersection.normal);
				double temp = sqrt(1.0 - ((pow(ray.intersection.mat->refractiveIndex, 2) 
								* (1 - pow(mddotn, 2))) / pow(1.0, 2)));

				refraDir = ((1/ray.intersection.mat->refractiveIndex) * 
							(ray.dir - (mddotn * ray.intersection.normal))) - temp * ray.intersection.normal;
				refraDir.normalize();

				if (temp >= 0)
					oddotn = refraDir.dot(ray.intersection.normal);
				else{
					rayFL.depth = ray.depth + 1;
					col = col + tc * shadeRay(rayFL);
					TIR = true;
				}

			}
						
			if (!TIR) {
				Ray3D rayTS = Ray3D(ray.intersection.point + 0.001 * refraDir, refraDir);
				double i = pow(ray.intersection.mat->refractiveIndex - 1, 2) / pow(ray.intersection.mat->refractiveIndex + 1, 2)
						 + (1 - pow(ray.intersection.mat->refractiveIndex - 1, 2) / pow(ray.intersection.mat->refractiveIndex + 1, 2)) 
						 * pow(1 - oddotn, 5);

				rayFL.depth = ray.depth + 1;
				Colour RLC = i * shadeRay(rayFL); 

				rayTS.depth = ray.depth + 1;
				Colour RC = (1 - i) * shadeRay(rayTS);

				col = col + tc * (RC + RLC);
			}

		}
			
	}

	col.clamp();
		
	return col; 
}	


void Raytracer::render( int width, int height, Point3D eye, Vector3D view, 
		Vector3D up, double fov, char* fileName, SceneDagNode* object, Raytracer* raytracer ) {
	Matrix4x4 viewToWorld;
	_scrWidth = width;
	_scrHeight = height;
	double factor = (double(height)/2)/tan(fov*M_PI/360.0);

	initPixelBuffer();
	viewToWorld = initInvViewMatrix(eye, view, up);
	if (MOTION_BLUR_ENABLED) {
		for (double t = 0.0; t < 1; t = t + 0.1 ) {
			raytracer->translate(object, Vector3D(t, 0, 0));

			// Construct a ray for each pixel.
			for (int i = 0; i < _scrHeight; i++) {
				for (int j = 0; j < _scrWidth; j++) {
					// Sets up ray origin and direction in view space, 
					// image plane is at z = -1.
					Point3D origin(0, 0, 0);
					Point3D imagePlane;

					// TODO: Convert ray to world space and call 
					// shadeRay(ray) to generate pixel colour. 
					Ray3D ray;
					// Point on plane in world space	
					int num_of_samples = 10;
					Colour col(0.0, 0.0, 0.0);

					// Implement Anti-alising feature
					if (ANTI_ALISING_ENABLED) {
						for (int k = 0; k < num_of_samples; k++) {
							double tx = rand() / (double) RAND_MAX - 0.5;
							double ty = rand() / (double) RAND_MAX - 0.5;

							// Shoot out more random ray
							imagePlane[0] = (-double(width)/2 + 0.5 + j + tx)/factor;
							imagePlane[1] = (-double(height)/2 + 0.5 + i + ty)/factor;
							imagePlane[2] = -1;

							ray.origin = viewToWorld * origin;
							ray.dir = viewToWorld * 
									Point3D(imagePlane[0], imagePlane[1], imagePlane[2]) - ray.origin;
							ray.dir.normalize();
		                    col = col + shadeRay(ray);
						}
						col = (1 / (double)num_of_samples) * col;
					} else {
						imagePlane[0] = (-double(width)/2 + 0.5 + j)/factor;
						imagePlane[1] = (-double(height)/2 + 0.5 + i)/factor;
						imagePlane[2] = -1;
						ray.origin = viewToWorld * origin;
						ray.dir = viewToWorld * 
									Point3D(imagePlane[0], imagePlane[1], imagePlane[2]) - ray.origin;
						ray.dir.normalize();

						col = col + shadeRay(ray); 
					}
					col.clamp();

					_rbuffer[i*width+j] = _rbuffer[i*width+j] + 0.1 * int(col[0]*255);
					_gbuffer[i*width+j] = _gbuffer[i*width+j] + 0.1 * int(col[1]*255);
					_bbuffer[i*width+j] = _bbuffer[i*width+j] + 0.1 * int(col[2]*255);
				}
			}
		}
	} else {
			// Construct a ray for each pixel.
			for (int i = 0; i < _scrHeight; i++) {
				for (int j = 0; j < _scrWidth; j++) {
					// Sets up ray origin and direction in view space, 
					// image plane is at z = -1.
					Point3D origin(0, 0, 0);
					Point3D imagePlane;

					// TODO: Convert ray to world space and call 
					// shadeRay(ray) to generate pixel colour. 
					Ray3D ray;
					// Point on plane in world space	
					int num_of_samples = 10;
					Colour col(0.0, 0.0, 0.0);

					// Implement Anti-alising feature
					if (ANTI_ALISING_ENABLED) {
						for (int k = 0; k < num_of_samples; k++) {
							double tx = rand() / (double) RAND_MAX - 0.5;
							double ty = rand() / (double) RAND_MAX - 0.5;

							// Shoot out more random ray
							imagePlane[0] = (-double(width)/2 + 0.5 + j + tx)/factor;
							imagePlane[1] = (-double(height)/2 + 0.5 + i + ty)/factor;
							imagePlane[2] = -1;

							ray.origin = viewToWorld * origin;
							ray.dir = viewToWorld * 
									Point3D(imagePlane[0], imagePlane[1], imagePlane[2]) - ray.origin;
							ray.dir.normalize();
		                    col = col + shadeRay(ray);
						}
						col = (1 / (double)num_of_samples) * col;
					} else {
						imagePlane[0] = (-double(width)/2 + 0.5 + j)/factor;
						imagePlane[1] = (-double(height)/2 + 0.5 + i)/factor;
						imagePlane[2] = -1;
						ray.origin = viewToWorld * origin;
						ray.dir = viewToWorld * 
									Point3D(imagePlane[0], imagePlane[1], imagePlane[2]) - ray.origin;
						ray.dir.normalize();

						col = col + shadeRay(ray); 
					}
					col.clamp();

					_rbuffer[i*width+j] = int(col[0]*255);
					_gbuffer[i*width+j] = int(col[1]*255);
					_bbuffer[i*width+j] = int(col[2]*255);
				}
			}

	}

	flushPixelBuffer(fileName);
}

int main(int argc, char* argv[])
{	
	// Build your scene and setup your camera here, by calling 
	// functions from Raytracer.  The code here sets up an example
	// scene and renders it from two different view points, DO NOT
	// change this if you're just implementing part one of the 
	// assignment.  
	Raytracer raytracer;
	int width = 320; 
	int height = 240; 

	if (argc == 3) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}
    //raytracer.enableShadows();
	// Camera parameters.
	Point3D eye(0, 0, 1);
	Vector3D view(0, 0, -1);
	Vector3D up(0, 1, 0);
	double fov = 60;

	// Defines a material for shading.
	Material gold( Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648), 
			Colour(0.628281, 0.555802, 0.366065), 
			51.2, 1.0, 0.0, 0.0 );
	Material jade( Colour(0, 0, 0), Colour(0.54, 0.89, 0.63), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8, 1.0, 0.0, 0.0 );

	Material glass( Colour(0.0, 0.0, 0.0), Colour(0.0, 0.0, 0.0),        
        Colour(0.0, 0.0, 0.0), 10.1, 1.05, 0.9, 0.0 );
	glass.refractiveIndex = 1.05;

	Material bronze( Colour(0.2125, 0.1275, 0.054), Colour(0.714, 0.4284, 0.18144), 
			Colour(0.393548, 0.271906, 0.166721), 
			0.0, 1.0, 0.0, 0.6 );

	Material red_rubber( Colour(0.05, 0.0, 0.0), Colour(0.5, 0.4, 0.4), 
			Colour(0.7, 0.04, 0.04), 
			0.0, 1.0, 0.0, 0.0 );

	Material copper( Colour(0.19125, 0.0735, 0.0225), Colour(0.7038, 0.27048, 0.0828), 
			Colour(0.256777, 0.137622, 0.086014), 
			0.1, 1.0, 0.0, 0.0 );

	Material chrome (Colour(0.25, 0.25, 0.25),Colour(0.4, 0.4, 0.4), 
			Colour(0.774597, 0.774597, 0.774597),16.8, 1.0, 0.0, 1.0);

	Material green (Colour(0.0, 0.25, 0.0),Colour(0.0, 1.0, 0.0), 
			Colour(0.0225,0.0225,0.0225),12.8, 1.0, 0.0, 0.7);

	// Defines a point light source.
	
	raytracer.addLightSource( new PointLight(Point3D(0, 0, 5), 
				Colour(0.9, 0.9, 0.9) ) );
	
	raytracer.addLightSource( new PointLight(Point3D(0, 0, -5), 
				Colour(0.9, 0.9, 0.9) ) );
			
	
       // raytracer.addLightSource( new AreaLight(Point3D(0, 0, 5), Vector3D(0, 0, 5),Vector3D(5, 0, 0),  Colour(0.9, 0.9, 0.9)));
	
	// Add a unit square into the scene with material mat.
	SceneDagNode* sphere = raytracer.addObject( new UnitSphere(), &red_rubber );
	SceneDagNode* cylinder = raytracer.addObject( new UnitCylinder(), &jade );
	SceneDagNode* plane1 = raytracer.addObject( new UnitSquare(), &glass );
	SceneDagNode* plane2 = raytracer.addObject( new UnitSquare(), &glass );
	SceneDagNode* plane3 = raytracer.addObject( new UnitSquare(), &glass );
	
	// Apply some transformations to the unit square.
	double factor1[3] = { 1.0, 1.0, 1.0 };
	double factor2[3] = { 10.0, 10.0, 10.0 };
	raytracer.translate(sphere, Vector3D(0, 0, -3));	
	raytracer.rotate(sphere, 'x', -45); 
	raytracer.rotate(sphere, 'z', 45); 
	raytracer.scale(sphere, Point3D(0, 0, 0), factor1);

	raytracer.translate(cylinder, Vector3D(1, 0, -5));	


	raytracer.translate(plane1, Vector3D(0, 0, -10));	
	raytracer.scale(plane1, Point3D(0, 0, 0), factor2);

	raytracer.translate(plane2, Vector3D(0, -5, -5));	
	raytracer.rotate(plane2, 'x', -90); 
	raytracer.scale(plane2, Point3D(0, 0, 0), factor2);

	raytracer.translate(plane3, Vector3D(-5, 0, -5));	
	raytracer.rotate(plane3, 'y', 90); 
	raytracer.scale(plane3, Point3D(0, 0, 0), factor2);

	// Render the scene, feel free to make the image smaller for
	// testing purposes.	
	if (ANTI_ALISING_ENABLED)
		raytracer.render(width, height, eye, view, up, fov, "ANTI_ALIASING_ENABLED-view1.bmp", sphere, &raytracer);
	else
		raytracer.render(width, height, eye, view, up, fov, "view1.bmp", sphere, &raytracer);
	
	// Render it from a different point of view.
	Point3D eye2(4, 2, 1);
	Vector3D view2(-4, -2, -6);
	if (ANTI_ALISING_ENABLED)
		raytracer.render(width, height, eye2, view2, up, fov, "ANTI_ALIASING_ENABLED-view2.bmp", sphere, &raytracer);
	else
		raytracer.render(width, height, eye2, view2, up, fov, "view2.bmp", sphere, &raytracer);
	
	return 0;
}

