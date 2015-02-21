#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <limits>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#ifdef _WIN32
static DWORD lastTime;
#else
static struct timeval lastTime;
#endif

#undef max
#undef min

#define PI 3.14159265  // Should be used from mathlib
inline float sqr(float x) { return x*x; }

using namespace std;

//****************************************************
// Helper Classes
//****************************************************

class Viewport {
public:
	int w, h; // width and height
};

class Point {
public:
	GLfloat x, y, z;
};

class Color {
public:
	GLfloat r, g, b;
};

class Patch {
public:
	Point p1, p2, p3, p4; //counterclockwise -- facing toward audience
	Color emission, reflectance, incident, excident;
};

class Line{
public:
	Point p1, p2, p3, p4;
};

class Quad {
public:
	Line l1, l2, l3, l4;
};

class Ray{
public:
	Point direction;
	Point position;
	GLfloat min;
	GLfloat max;
};

class Node {
public:
	Patch* patch;
	std::vector<Node*> children;
	std::vector<Node*> potentialChildren;
};

//****************************************************
// Global Variables
//****************************************************
Viewport	viewport;
Color emission, reflectance, incident, excident;
std::string filename;

std::vector<Node*> nodes;
std::vector<Patch*> patches;
std::vector<Patch*> hemiPatches;
std::vector<Point> vertices;
int maxPass;

GLfloat stepSize;
GLfloat hemiResolution;
bool outfile = false;
bool adaptiveBool = false;
// Wired Mode or Filled Mode
bool wired = false;
bool smooth = true;
/*
GLfloat cameraX = 0.0;
GLfloat cameraY = 0.0;
GLfloat cameraZ = 10.0;

GLfloat lx = 0.0;
GLfloat ly = 0.0;
GLfloat lz = -1.0;

GLfloat cameraAngle = 0.0;
*/
GLfloat yRot = 0.0;
GLfloat xRot = 0.0;
GLfloat xTran = 0.0;
GLfloat yTran = 0.0;

GLfloat scaleValue = 1.0;
GLfloat maxX = 10;
GLfloat maxY = 10;


//****************************************************
// reshape viewport if the window is resized
//****************************************************
void myReshape(int w, int h) {
	viewport.w = w;
	viewport.h = h;

	glViewport (0,0,viewport.w,viewport.h);
	glMatrixMode(GL_PROJECTION);

	glLoadIdentity();
	gluPerspective(90, 1, 0.01, 1000);
	gluLookAt(0.0, 0.0, -10.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

	glMatrixMode(GL_MODELVIEW);
}

//****************************************************
// Simple init function
//****************************************************

void initScene(){
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f); // Clear to black, fully transparent
	myReshape(viewport.w,viewport.h);

	glEnable(GL_NORMALIZE);
	glEnable(GL_DEPTH_TEST);
}

//*********************************************w
// Helper Methods
//********************************************* 

Point multiplyPoint(GLfloat s, Point p){
	Point r;
	r.x = p.x * s;
	r.y = p.y * s;
	r.z = p.z * s;
	return r;
}

Point dividePoint(GLfloat s, Point p){
	// s cannot be zero
	Point r;         
	r.x = p.x / s;
	r.y = p.y / s;
	r.z = p.z / s;
	return r;
}

Point addPoint(Point p1, Point p2){
	Point r;
	r.x = p1.x + p2.x;
	r.y = p1.y + p2.y;
	r.z = p1.z + p2.z;
	return r;
}

Point subtractPoint(Point p1, Point p2){
	Point r;
	r.x = p1.x - p2.x;
	r.y = p1.y - p2.y;
	r.z = p1.z - p2.z;
	return r;
}

Color subtractColorAbs(Color color1, Color color2){
    Color color;
    color.r = abs(color1.r - color2.r);
    color.g = abs(color1.g - color2.g);
    color.b = abs(color1.b - color2.b);
    return color;
}

bool colorLessThan(GLfloat threshold, Color color){
    return color.r < threshold && color.g < threshold && color.b < threshold;
}

Patch quadToPatch(Quad quad, Color emmi, Color ref, Color in, Color ex){
	Patch patch;
	patch.p1 = quad.l1.p1;
	patch.p2 = quad.l1.p4;
	patch.p3 = quad.l4.p1;
	patch.p4 = quad.l4.p4;
	patch.emission = emmi;
	patch.reflectance = ref;
	patch.incident = in;
	patch.excident = ex;

	return patch;
}

Quad patchToQuad(Patch patch){
	Quad quad;
	Line l1, l2, l3, l4;

	Point ul, ur, ll, lr;
	ul = patch.p1;
	ll = patch.p2;
	lr = patch.p3;
	ur = patch.p4;

	Point pp1, pp2, pp3, pp4, pp5, pp6, pp7, pp8, pp9, pp10, pp11, pp12, pp13, pp14, pp15, pp16;
	pp1 = ul;
	pp4 = ur;
	pp13 = ll;
	pp16 = lr;

	Point hori1 = subtractPoint(pp4, pp1);
	hori1 = dividePoint(3, hori1);
	pp2 = addPoint(pp1, hori1);
	pp3 = addPoint(pp1, multiplyPoint(2, hori1));

	Point hori2 = subtractPoint(pp16, pp13);
	hori2 = dividePoint(3, hori2);
	pp14 = addPoint(pp13, hori2);
	pp15 = addPoint(pp13, multiplyPoint(2, hori2));

	Point vert1 = subtractPoint(pp1, pp13);
	vert1 = dividePoint(3, vert1);
	pp9 = addPoint(pp13, vert1);
	pp5 = addPoint(pp13, multiplyPoint(2, vert1));

	Point vert2 = subtractPoint(pp4, pp16);
	vert2 = dividePoint(3, vert2);
	pp12 = addPoint(pp16, vert2);
	pp8 = addPoint(pp16, multiplyPoint(2, vert2));

	Point mid1 = subtractPoint(pp8, pp5);
	mid1 = dividePoint(3, mid1);
	pp6 = addPoint(pp5, mid1);
	pp7 = addPoint(pp5, multiplyPoint(2, mid1));

	Point mid2 = subtractPoint(pp12, pp9);
	mid2 = dividePoint(3, mid2);
	pp10 = addPoint(pp9, mid2);
	pp11 = addPoint(pp9, multiplyPoint(2, mid2));

	l1.p1 = pp1;
	l1.p2 = pp2;
	l1.p3 = pp3;
	l1.p4 = pp4;

	l2.p1 = pp5;
	l2.p2 = pp6;
	l2.p3 = pp7;
	l2.p4 = pp8;

	l3.p1 = pp9;
	l3.p2 = pp10;
	l3.p3 = pp11;
	l3.p4 = pp12;

	l4.p1 = pp13;
	l4.p2 = pp14;
	l4.p3 = pp15;
	l4.p4 = pp16;

	quad.l1 = l1;
	quad.l2 = l2;
	quad.l3 = l3;
	quad.l4 = l4;

	return quad; 
}

Point normalize(Point p){
	Point r;
	GLfloat len = sqrt(pow(p.x, 2.0) + pow(p.y, 2.0) + pow(p.z, 2.0));
	r.x = p.x / len;
	r.y = p.y / len;
	r.z = p.z / len;
	return r;
}

Point crossProduct(Point p1, Point p2){
	Point r;
	r.x = p1.y * p2.z - p1.z * p2.y;
	r.y = p1.z * p2.x - p1.x * p2.z;
	r.z = p1.x * p2.y - p1.y * p2.x;
	return r;
}

GLfloat dotProduct(Point p1, Point p2){
	return (p1.x * p2.x) + (p1.y * p2.y) + (p1.z * p2.z);
}

Point getNormal(Point p1, Point p2, Point p3){
	return crossProduct(subtractPoint(p3, p1), subtractPoint(p2, p1));
}

Point midPoint(Point p1, Point p2) {
	Point r;
	r.x = (p1.x + p2.x)/2.0;
	r.y = (p1.y + p2.y)/2.0;
	r.z = (p1.z + p2.z)/2.0;
	return r;
}

GLfloat distancePoint(Point p1, Point p2) {
	return sqrt(pow(p1.x-p2.x, 2.0) + pow(p1.y-p2.y, 2.0) + pow(p1.z-p2.z, 2.0));
}

Point bernstein(GLfloat u, Line line){
	Point a, b, c, d, e, p;

	a = addPoint(multiplyPoint((1.0-u), line.p1), multiplyPoint(u, line.p2));
	b = addPoint(multiplyPoint((1.0-u), line.p2), multiplyPoint(u, line.p3));
	c = addPoint(multiplyPoint((1.0-u), line.p3), multiplyPoint(u, line.p4));

	d = addPoint(multiplyPoint((1.0-u), a), multiplyPoint(u, b));
	e = addPoint(multiplyPoint((1.0-u), b), multiplyPoint(u, c));

	p = addPoint(multiplyPoint((1.0-u), d), multiplyPoint(u, e));

	return p;
}

Point quadPoint(GLfloat u, GLfloat v, Quad quad) {
	// assume all input patches are quadrilaterial
	Line vline;
	Point p;

	vline.p1 = bernstein(u, quad.l1);
	vline.p2 = bernstein(u, quad.l2);
	vline.p3 = bernstein(u, quad.l3);
	vline.p4 = bernstein(u, quad.l4);

	p = bernstein(v, vline);

	return p;
}

void drawPolygon(Patch patch){
	glBegin(GL_POLYGON);
	glColor3f(patch.excident.r, patch.excident.g, patch.excident.b);
	glVertex3f(patch.p1.x, patch.p1.y, patch.p1.z);
	glVertex3f(patch.p2.x, patch.p2.y, patch.p2.z);
	glVertex3f(patch.p3.x, patch.p3.y, patch.p3.z);
	glVertex3f(patch.p4.x, patch.p4.y, patch.p4.z);
	glEnd();
}
/*
void drawPolygon2(Patch patch){
glBegin(GL_POLYGON);
glColor3f(1, 1, 1);
glVertex3f(patch.p1.x, patch.p1.y, patch.p1.z);
glVertex3f(patch.p2.x, patch.p2.y, patch.p2.z);
glVertex3f(patch.p3.x, patch.p3.y, patch.p3.z);
glVertex3f(patch.p4.x, patch.p4.y, patch.p4.z);
glEnd();
}
*/
Point patchNormal(Patch patch){
	Point normal = normalize(crossProduct(subtractPoint(patch.p1, patch.p4), subtractPoint(patch.p3, patch.p4)));
	return normal;
}

Point patchCenter(Patch patch){
	Point center = dividePoint(4, addPoint(patch.p1, addPoint(patch.p2, addPoint(patch.p3, patch.p4))));
	return center;
}

Quad getHemiQuad(Patch patch, GLfloat side){ 
	Point normal = dividePoint(2, patchNormal(patch)); 
	Point horizontal = dividePoint(2, normalize(subtractPoint(patch.p4, patch.p1))); 
	Point vertical = dividePoint(2, normalize(crossProduct(normal, horizontal))); 
	Point center = patchCenter(patch); 

	Patch* hemiPatch = new Patch; 
	if (side == 0.0){           // front scene 
		hemiPatch->p1 = subtractPoint(addPoint(vertical, addPoint(center, normal)), horizontal); 
		hemiPatch->p2 = subtractPoint(subtractPoint(addPoint(center, normal), horizontal), vertical); 
		hemiPatch->p3 = subtractPoint(addPoint(addPoint(center, normal), horizontal), vertical); 
		hemiPatch->p4 = addPoint(addPoint(addPoint(center, normal), vertical), horizontal); 
	}else if (side == 1.0){     // right scene 
		hemiPatch->p1 = subtractPoint(addPoint(vertical, addPoint(center, normal)), horizontal); 
		hemiPatch->p2 = subtractPoint(subtractPoint(addPoint(center, normal), horizontal), vertical); 
		hemiPatch->p3 = subtractPoint(subtractPoint(subtractPoint(center, normal), horizontal), vertical); 
		hemiPatch->p4 = subtractPoint(addPoint(vertical, subtractPoint(center, normal)), horizontal); 
	}else if (side == 2.0){     // left scene 
		hemiPatch->p1 = addPoint(addPoint(addPoint(center, normal), vertical), horizontal); 
		hemiPatch->p2 = subtractPoint(addPoint(addPoint(center, normal), horizontal), vertical); 
		hemiPatch->p3 = subtractPoint(addPoint(subtractPoint(center, normal), horizontal), vertical); 
		hemiPatch->p4 = addPoint(addPoint(subtractPoint(center, normal), vertical), horizontal); 
	}else if (side == 3.0){     // top scene 
		hemiPatch->p1 = subtractPoint(addPoint(vertical, addPoint(center, normal)), horizontal);  
		hemiPatch->p2 = subtractPoint(addPoint(vertical, subtractPoint(center, normal)), horizontal); 
		hemiPatch->p3 = addPoint(addPoint(subtractPoint(center, normal), vertical), horizontal); 
		hemiPatch->p4 = addPoint(addPoint(addPoint(center, normal), vertical), horizontal); 
	}else if (side == 4.0){     // bottom scene 
		hemiPatch->p1 = subtractPoint(subtractPoint(addPoint(center, normal), horizontal), vertical); 
		hemiPatch->p2 = subtractPoint(subtractPoint(subtractPoint(center, normal), horizontal), vertical); 
		hemiPatch->p3 = subtractPoint(addPoint(subtractPoint(center, normal), horizontal), vertical); 
		hemiPatch->p4 = subtractPoint(addPoint(addPoint(center, normal), horizontal), vertical);
	} 

	hemiPatches.push_back(hemiPatch); 

	Quad hemiQuad = patchToQuad(*hemiPatch); 
	return hemiQuad; 
}

bool arePointsEquivalent(Point p1, Point p2){
	return (p1.x == p2.x) && (p1.y == p2.y) && (p1.z == p2.z);
}

Ray generateRay(Quad hemiQuad, Point camera, GLfloat u, GLfloat v){
	Ray ray;
	ray.direction = subtractPoint(quadPoint(u+hemiResolution/2, v+hemiResolution/2, hemiQuad), camera);
	ray.position = camera;
	ray.min = 0.001;
	ray.max = std::numeric_limits<float>::max();
	return ray;
}

bool arePatchesEquivalent(Patch patch1, Patch patch2){
	return arePointsEquivalent(patch1.p1, patch2.p1) && arePointsEquivalent(patch1.p2, patch2.p2) && arePointsEquivalent(patch1.p3, patch2.p3) && arePointsEquivalent(patch1.p4, patch2.p4);
}

std::vector<Node*> leafNodes(vector<Node*> nodes) {
	std::vector<Node*> leaves;
	for (int i = 0; i < nodes.size(); i++) {
		if (nodes[i]->children.size() > 0) {
			std::vector<Node*> temp = leafNodes(nodes[i]->children);
			//leaves.reserve(leaves.size() + temp.size());
			leaves.insert(leaves.end(), temp.begin(), temp.end());
		} else {
			leaves.push_back(nodes[i]);
		}
	}
	return leaves;
}


GLfloat intersection(Ray ray, Patch patch) {

	Point normal = patchNormal(patch);
	float denom = dotProduct(ray.direction, normal);

	if (denom == 0.0){      // the direction 
		return std::numeric_limits<float>::max();
	}

	float num = dotProduct(subtractPoint(patch.p1, ray.position), normal);

	float timeHit = num / denom; 

	if (timeHit < ray.min || timeHit > ray.max){  // hit time outside of boundary
		return std::numeric_limits<float>::max();
	}

	Point hitPoint = addPoint(ray.position, multiplyPoint(timeHit, ray.direction));

	Point edge21 = subtractPoint(patch.p2, patch.p1);
	Point edge32 = subtractPoint(patch.p3, patch.p2);
	Point edge43 = subtractPoint(patch.p4, patch.p3);
	Point edge14 = subtractPoint(patch.p1, patch.p4);

	Point edge12 = subtractPoint(patch.p1, patch.p2);
	Point edge23 = subtractPoint(patch.p2, patch.p3);
	Point edge34 = subtractPoint(patch.p3, patch.p4);
	Point edge41 = subtractPoint(patch.p4, patch.p1);

	Point hitPointToP1 = subtractPoint(hitPoint, patch.p1);
	Point hitPointToP2 = subtractPoint(hitPoint, patch.p2);
	Point hitPointToP3 = subtractPoint(hitPoint, patch.p3);
	Point hitPointToP4 = subtractPoint(hitPoint, patch.p4);

	float leftOfEdge1 = dotProduct(crossProduct(edge21, edge41), crossProduct(edge21, hitPointToP1)); 
	float leftOfEdge2 = dotProduct(crossProduct(edge32, edge12), crossProduct(edge32, hitPointToP2));
	float leftOfEdge3 = dotProduct(crossProduct(edge43, edge23), crossProduct(edge43, hitPointToP3));
	float leftOfEdge4 = dotProduct(crossProduct(edge14, edge34), crossProduct(edge14, hitPointToP4));


	if (leftOfEdge1 < 0.0 || leftOfEdge2 < 0.0 || leftOfEdge3 < 0.0 || leftOfEdge4 < 0.0){   // hit point not inside the patch 
		return std::numeric_limits<float>::max();    
	}   

	return timeHit;
}

Node* intersectChildren(Ray ray, Patch referencePatch, std::vector<Node*> children){ 
	float minTimeHit = std::numeric_limits<float>::max(); 
	float timeHit;
	Node* bestChild = NULL;
	for (int i = 0; i < children.size(); i++){ 
		if (!arePatchesEquivalent(referencePatch, *(children[i]->patch))){ 
			timeHit = intersection(ray, *(children[i]->patch)); 
			if (timeHit < std::numeric_limits<float>::max()){ 
				if (timeHit < minTimeHit){ 
					minTimeHit = timeHit; 
					bestChild = children[i];
				} 
			} 
		}
	} 
	
	if (minTimeHit == std::numeric_limits<float>::max()) { 
		bestChild = NULL;
	}
	
	return bestChild; 
}

Color intersectAllPatch(Ray ray, Patch referencePatch){
	Node* bestParentNode = intersectChildren(ray, referencePatch, nodes);

	Color color;
	color.r = 0.0;
	color.g = 0.0;
	color.b = 0.0;
	
	if (bestParentNode != NULL) {
		while (bestParentNode != NULL && bestParentNode->children.size() > 0) {
			bestParentNode = intersectChildren(ray, referencePatch, bestParentNode->children);
		}

		if (bestParentNode != NULL) {
			color = bestParentNode->patch->excident;
		}
	}
	return color;
}

std::vector<Node*> uniformSubdivisionSinglePatch(Patch patch, GLfloat stepSize){
	
	Quad quad = patchToQuad(patch);
	Point new1, new2, new3, new4;

	std::vector<Node*> children;

	GLfloat old_u, new_u, old_v, new_v;

	old_u = 0.0;
	for (GLfloat u = 0.0; u < 1.0; u += stepSize) {
		new_u = old_u + stepSize;
		if (new_u > 1.0){
			new_u = 1.0;
		}
		old_v = 0.0;
		for (GLfloat v = 0.0; v < 1.0; v += stepSize) {
			new_v = old_v + stepSize;
			if (new_v > 1.0){
				new_v = 1.0;
			}
			new1 = quadPoint(old_u, old_v, quad);
			new2 = quadPoint(old_u, new_v, quad);
			new3 = quadPoint(new_u, new_v, quad);
			new4 = quadPoint(new_u, old_v, quad);

			//vector of patch.add
			Patch* newPatch = new Patch;
			newPatch->p1 = new1;
			newPatch->p2 = new2;
			newPatch->p3 = new3;
			newPatch->p4 = new4;

			newPatch->emission = patch.emission;
			newPatch->reflectance = patch.reflectance;
			newPatch->incident = patch.incident;
			newPatch->excident = patch.excident;

			Node* node = new Node;
			node->patch = newPatch;

			children.push_back(node);

			old_v = new_v;
		}
		old_u = new_u;
	}

	return children;
}

void uniformTesselation(){
	for(int i = 0; i < nodes.size(); i++){
		nodes[i]->children = uniformSubdivisionSinglePatch(*(nodes[i]->patch), stepSize);
	}
}

void drawScene() {
	std::vector<Node*> allLeafNodes = leafNodes(nodes);

	for (int i = 0; i < allLeafNodes.size(); i++) {
		drawPolygon(*(allLeafNodes[i])->patch);
	}
	/*
	for(unsigned int i = 0; i < hemiPatches.size(); i++){
	drawPolygon2(*(hemiPatches[i]));
	}
	*/
}

GLfloat formFactor(Ray ray, GLfloat side, GLfloat v){
	GLfloat weight;
	GLfloat distancePowTwo = ray.direction.x * ray.direction.x + ray.direction.y * ray.direction.y + ray.direction.z * ray.direction.z;
	GLfloat z = 0.5 - v;
	if (side == 0.0){
		weight = 0.5 / (distancePowTwo * distancePowTwo * PI);
	} else if (side == 1.0){
		weight = z / (distancePowTwo * distancePowTwo * PI);
	} else if (side == 2.0){
		weight = z / (distancePowTwo * distancePowTwo * PI);
	} else if (side == 3.0){
		weight = z / (distancePowTwo * distancePowTwo * PI);
	} else if (side == 4.0){
		weight = z / (distancePowTwo * distancePowTwo * PI);
	}
	return weight;
}

void getIncidentForPatch(Patch* patch) {
	std::vector<Color> visibleColors;
	GLfloat uHemiMax, vHemiMax;
	Point center = patchCenter(*patch);

	for (int i = 0; i < 5; i++) {
		
		Quad quad = getHemiQuad(*patch, i);

		if (i == 0) {
			uHemiMax = 1;
			vHemiMax = 1;
		} else {
			uHemiMax = 1;
			vHemiMax = 0.5;
		}

		for (GLfloat u = 0.0; u <= uHemiMax; u += hemiResolution) {
			for (GLfloat v = 0.0; v <= vHemiMax; v += hemiResolution) {
				Ray ray = generateRay(quad, center, u, v);
				Color color = intersectAllPatch(ray, *patch);
				GLfloat factor = formFactor(ray, i, v);
				color.r = color.r * factor;
				color.g = color.g * factor;
				color.b = color.b * factor; 
				visibleColors.push_back(color);
			}
		}
	}
	Color* avgColor = new Color;
	avgColor->r = 0;
	avgColor->g = 0;
	avgColor->b = 0;

	for (int i = 0; i < visibleColors.size(); i++) {
		avgColor->r += visibleColors[i].r;
		avgColor->g += visibleColors[i].g;
		avgColor->b += visibleColors[i].b;
	}

	avgColor->r /= visibleColors.size();
	avgColor->g /= visibleColors.size();
	avgColor->b /= visibleColors.size();

	patch->incident = *avgColor;
}

void updateExcidentForPatch(Patch* patch) {
	patch->excident.r = (patch->incident.r * patch->reflectance.r) + patch->emission.r;
	patch->excident.g = (patch->incident.g * patch->reflectance.g) + patch->emission.g;
	patch->excident.b = (patch->incident.b * patch->reflectance.b) + patch->emission.b;

	//if (patch->excident.r > 1.0) patch->excident.r = 1.0;
	//if (patch->excident.g > 1.0) patch->excident.g = 1.0;
	//if (patch->excident.b > 1.0) patch->excident.b = 1.0;
}

bool thresholdComparison(Node* node){
    Color color1 = node->potentialChildren[0]->patch->incident;
    Color color2 = node->potentialChildren[1]->patch->incident;
    Color color3 = node->potentialChildren[2]->patch->incident;
    Color color4 = node->potentialChildren[3]->patch->incident;

    Color diff1 = subtractColorAbs(color1, node->patch->incident);
    Color diff2 = subtractColorAbs(color2, node->patch->incident);
    Color diff3 = subtractColorAbs(color3, node->patch->incident);
    Color diff4 = subtractColorAbs(color4, node->patch->incident);

    return colorLessThan(stepSize, diff1) && colorLessThan(stepSize, diff2) && colorLessThan(stepSize, diff3) && colorLessThan(stepSize, diff4);
    
}

void adaptiveSubdivisionIncident(Node* node){
    std::vector<Node*> newChildren = uniformSubdivisionSinglePatch(*(node->patch), 0.5);
    for (int i = 0; i < newChildren.size(); i++) {
        getIncidentForPatch(newChildren[i]->patch);
    }
    node->potentialChildren = newChildren;
    getIncidentForPatch(node->patch);
}

void adaptiveSubdivisionExcident(Node* node){
    if (!thresholdComparison(node)){
        node->children = node->potentialChildren;
        for (int i = 0; i < node->children.size(); i++) {
            updateExcidentForPatch(node->children[i]->patch);
        }
    }else{
        updateExcidentForPatch(node->patch);
    }
}

void radiosity() {
	for(int passIndex = 0; passIndex < maxPass; passIndex++) {
		cout << "Starting pass " << passIndex+1 << " of " << maxPass << ".\n";

		std::vector<Node*> allLeafNodes = leafNodes(nodes);

		for (int i = 0; i < allLeafNodes.size(); i++) {
			cout << "Starting patch " << i+1 << " of " << allLeafNodes.size() << ".\n";
			if (adaptiveBool) {
				adaptiveSubdivisionIncident(allLeafNodes[i]);
			} else {
				getIncidentForPatch(allLeafNodes[i]->patch);
			}
		}
		for (int i = 0; i < allLeafNodes.size(); i++) {
			if (adaptiveBool) {
				adaptiveSubdivisionExcident(allLeafNodes[i]);
			} else {
				updateExcidentForPatch(allLeafNodes[i]->patch);
			}
		}
		cout << "Pass " << passIndex+1 << " of " << maxPass << " done.\n";
	}
}
//****************************************************
// File Parser
//****************************************************
void loadScene(std::string file) {
	incident.r = excident.r = 0.0;
	incident.g = excident.g = 0.0;
	incident.b = excident.b = 0.0;

	std::ifstream inpfile(file.c_str());
	if(!inpfile.is_open()) {
		std::cout << "Unable to open file" << std::endl; 
	} else {
		std::string line;

		while(inpfile.good()) {
			std::vector<std::string> splitline;
			std::string buf;

			std::getline(inpfile,line);
			std::stringstream ss(line);

			while (ss >> buf) {
				splitline.push_back(buf);
			}
			//Ignore blank lines
			if(splitline.size() == 0) {
				continue;
			}
			if (!outfile) {
				//Ignore comments
				if(splitline[0][0] == '#') {
					continue;
				}
				//max passes
				//  max # of radiosity passes (default 5)
				else if(!splitline[0].compare("maxpass")) {
					// maxdepth: atoi(splitline[1].c_str())
					maxPass = atoi(splitline[1].c_str());
				}
				//output filename
				//  output file to write patch data to 
				else if(!splitline[0].compare("output")) {
					filename = splitline[1];
				}
				//vertex x y z
				//  Deﬁnes a vertex at the given location.
				//  The vertex is put into a pile, starting to be numbered at 0.
				else if(!splitline[0].compare("vertex")) {
					// x: atof(splitline[1].c_str()),
					// y: atof(splitline[2].c_str()),
					// z: atof(splitline[3].c_str()));
					// Create a new vertex with these 3 values, store in some array
					Point* p = new Point;
					p->x = atof(splitline[1].c_str());
					p->y = atof(splitline[2].c_str());
					p->z = atof(splitline[3].c_str());
					vertices.push_back(*p);
				}
				//quad a b c d
				//	Defines a quad at the given location
				//	Points are given in counterclockwise rotation
				else if (!splitline[0].compare("quad")) {
					// a: atof(splitline[1].c_str()),
					// b: atof(splitline[2].c_str()),
					// c: atof(splitline[3].c_str()));
					// d: atof(splitline[4].c_str()));
					// Create a new quad with these 4 values, store in Patch with current emission, reflectance settings
					Patch* patch = new Patch;
					patch->p1 = vertices[atof(splitline[1].c_str())];
					patch->p2 = vertices[atof(splitline[2].c_str())];
					patch->p3 = vertices[atof(splitline[3].c_str())]; 
					patch->p4 = vertices[atof(splitline[4].c_str())];

					patch->emission = emission;
					patch->reflectance = reflectance;
					patch->incident = incident;
					patch->excident = excident;

					//patches.push_back(patch);

					Node* node = new Node;
					node->patch = patch;

					nodes.push_back(node);
				}
				//shininess s
				//  speciﬁes the shininess of the surface.
				else if(!splitline[0].compare("reflectance")) {
					// shininess: atof(splitline[1].c_str())
					// Update current properties
					reflectance.r = atof(splitline[1].c_str());
					reflectance.g = atof(splitline[2].c_str());
					reflectance.b = atof(splitline[3].c_str());
				}
				//emission r g b
				//  gives the emissive color of the surface.
				else if(!splitline[0].compare("emission")) {
					// r: atof(splitline[1].c_str())
					// g: atof(splitline[2].c_str())
					// b: atof(splitline[3].c_str())
					// Update current properties
					emission.r = atof(splitline[1].c_str());
					emission.g = atof(splitline[2].c_str());
					emission.b = atof(splitline[3].c_str());

					excident.r = atof(splitline[1].c_str());
					excident.g = atof(splitline[2].c_str());
					excident.b = atof(splitline[3].c_str());
				} else {
					std::cerr << "Unknown command: " << splitline[0] << std::endl;
				}
			} else {
				Point* p1 = new Point;
				Point* p2 = new Point;
				Point* p3 = new Point;
				Point* p4 = new Point;

				p1->x = atof(splitline[0].c_str());
				p1->y = atof(splitline[1].c_str());
				p1->z = atof(splitline[2].c_str());

				p2->x = atof(splitline[3].c_str());
				p2->y = atof(splitline[4].c_str());
				p2->z = atof(splitline[5].c_str());

				p3->x = atof(splitline[6].c_str());
				p3->y = atof(splitline[7].c_str());
				p3->z = atof(splitline[8].c_str());

				p4->x = atof(splitline[9].c_str());
				p4->y = atof(splitline[10].c_str());
				p4->z = atof(splitline[11].c_str());

				Patch* patch = new Patch;
				patch->p1 = *p1;
				patch->p2 = *p2;
				patch->p3 = *p3;
				patch->p4 = *p4;

				Color* color = new Color;

				color->r = atof(splitline[12].c_str());
				color->g = atof(splitline[13].c_str());
				color->b = atof(splitline[14].c_str());

				patch->excident = *color;

				Node* node = new Node;
				node->patch = patch;

				nodes.push_back(node);
			}
		}
		inpfile.close();
	}
}

//****************************************************
// Output Generator
//****************************************************
void generateOutput() {

	std::vector<Node*> allLeafNodes = leafNodes(nodes);

	ofstream output;
	output.open(filename.c_str());
	for (int i = 0; i < allLeafNodes.size(); i++) {
		Patch patch = *allLeafNodes[i]->patch;
		output << patch.p1.x << " ";
		output << patch.p1.y << " ";
		output << patch.p1.z << " ";

		output << patch.p2.x << " ";
		output << patch.p2.y << " ";
		output << patch.p2.z << " ";

		output << patch.p3.x << " ";
		output << patch.p3.y << " ";
		output << patch.p3.z << " ";

		output << patch.p4.x << " ";
		output << patch.p4.y << " ";
		output << patch.p4.z << " ";

		output << patch.excident.r << " ";
		output << patch.excident.g << " ";
		output << patch.excident.b << "\n";
	}
	output.close();
}

//****************************************************
// function that does the actual drawing of stuff
//***************************************************
void myDisplay() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// clear the color buffer

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	/*
	gluLookAt(cameraX,      cameraY,      cameraZ,
	cameraX + lx, cameraY + ly, cameraZ + lz,
	0.0, 1.0, 0.0);
	*/
	glScalef(scaleValue, scaleValue, scaleValue);

	glTranslatef(xTran, 0.0, 0.0);
	glTranslatef(0.0, yTran, 0.0);

	glRotatef(xRot, 1.0, 0.0, 0.0);
	glRotatef(yRot, 0.0, 1.0, 0.0);

	drawScene();

	glFlush();
	glutSwapBuffers();					// swap buffers (we earlier set double buffer)
}

void keyboard(unsigned char key, int x, int y) {
	switch (key)  {
	case 32: // Space key
		exit (0);
		break;
	case 61: //+ key
		scaleValue += 0.1;
		break;
	case 45: // - key
		scaleValue -= 0.1;
		break;
	case 115: //s key
		if (smooth){
			glShadeModel(GL_FLAT);
			smooth = false;
		}else{
			glShadeModel(GL_SMOOTH);
			smooth = true;
		}
		break;
	case 49: //1 key
		if (!wired){
			glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
			wired = true;
		}else{
			glPolygonMode( GL_FRONT_AND_BACK, GL_FILL);
			wired = false;
		}
		break;
	}
}

void SpecialKeys(int key, int x, int y) {

	GLfloat fraction = 0.2;

	switch (key) {
	case GLUT_KEY_LEFT:
		if (!(glutGetModifiers() & GLUT_ACTIVE_SHIFT)) {
			/*
			cameraAngle -= 0.01;
			lx = sin(cameraAngle);
			lz = -cos(cameraAngle);
			*/
			yRot += 10;
		} else {
			xTran -= 0.5;
		}
		break;
	case GLUT_KEY_RIGHT:
		if (!(glutGetModifiers() & GLUT_ACTIVE_SHIFT)) {
			/*
			cameraAngle += 0.01;
			lx = sin(cameraAngle);
			lz = -cos(cameraAngle);
			*/
			yRot -= 10;
		} else {
			xTran += 0.5;
		}
		break;
	case GLUT_KEY_UP:
		if (!(glutGetModifiers() & GLUT_ACTIVE_SHIFT)) {
			/*
			cameraAngle += 0.01;
			ly = sin(cameraAngle);
			lz = -cos(cameraAngle);
			*/
			xRot -= 10;
		} else {
			/*
			cameraX += lx * fraction;
			cameraY += ly * fraction;
			cameraZ += lz * fraction;
			*/
			yTran += 0.5;
		}
		break;
	case GLUT_KEY_DOWN:
		if (!(glutGetModifiers() & GLUT_ACTIVE_SHIFT)) {
			/*
			cameraAngle -= 0.01;
			ly = sin(cameraAngle);
			lz = -cos(cameraAngle);
			*/
			xRot += 10;
		} else {
			/*
			cameraX -= lx * fraction;
			cameraY -= ly * fraction;
			cameraZ -= lz * fraction;
			*/
			yTran -= 0.5;
		}
		break;
	}
}

void myFrameMove() {
	//nothing here for now
#ifdef _WIN32
	Sleep(10);                                   //give ~10ms back to OS (so as not to waste the CPU)
#endif
	glutPostRedisplay(); // forces glut to call the display function (myDisplay())
}

//****************************************************
// the usual stuff, nothing exciting here
//****************************************************
int main(int argc, char *argv[]) {

	if (argc > 5) {
		outfile = true;
	}

	std::string filename = argv[1];
	char* adaptiveUniform = argv[2];
	stepSize = 1.0/atof(argv[3]);
	hemiResolution = 1.0/atof(argv[4]);

	if (adaptiveUniform[0] == 'a') {
		adaptiveBool = true;
	}

	loadScene(filename);

	if (!outfile) {
		if(!adaptiveBool) {
			uniformTesselation();
		}
		radiosity();
		generateOutput();
	}

	//This initializes glut
	glutInit(&argc, argv);

	//This tells glut to use a double-buffered window with red, green, and blue channels 
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

	// Initalize theviewport size
	viewport.w = 1000;
	viewport.h = 1000;

	//The size and position of the window
	glutInitWindowSize(viewport.w, viewport.h);
	glutInitWindowPosition(0,0);
	glutCreateWindow(argv[0]);

	initScene();							// quick function to set up scene

	glutDisplayFunc(myDisplay);				// function to run when its time to draw something
	glutReshapeFunc(myReshape);				// function to run when the window gets resized
	glutIdleFunc(myFrameMove);	
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(SpecialKeys);

	glutMainLoop();							// infinite loop that will keep drawing and resizing
	// and whatever else

	return 0;
}

