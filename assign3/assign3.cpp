/*
CSCI 420 Computer Graphics
Assignment 3: Ray Tracer
Name: <Ching-Chih Chen>
*/

#include <windows.h>
#include <stdlib.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <stdio.h>
#include <string>

#include "opencv2/core/core.hpp"
#include "opencv2/core/utility.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/imgcodecs/imgcodecs.hpp"

const int MAX_TRIANGLES = 2000;
const int MAX_SPHERES = 10;
const int MAX_LIGHTS = 10;

char* filename = 0;

//different display modes
enum class PLOTMODE { MODE_DISPLAY = 1, MODE_JPEG = 2 };
PLOTMODE mode = PLOTMODE::MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
//const int WIDTH = 960, HEIGHT = 720;
const int WIDTH = 640, HEIGHT = 480;
//const int WIDTH = 320, HEIGHT = 240;
const float ASPECT_RATIO = static_cast<float>(WIDTH) / static_cast<float>(HEIGHT);

//the field of view of the camera
const double fov = 120.0, Pi = 3.1415926, fovPi = (fov / 360.0) * Pi;

//boundary coordinates
double yMax = tan(fovPi / 2), xMax = ASPECT_RATIO * yMax;

unsigned char buffer[HEIGHT][WIDTH][3];

struct vec3 {
	double x, y, z;

	vec3(const double x0 = 0, const double y0 = 0, const double z0 = 0) : x(x0), y(y0), z(z0) {}
	void normalize() {
		double scalar = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
		x /= scalar, y /= scalar, z /= scalar;
	}
};

struct Vertex {
	double position[3];
	double color_diffuse[3];
	double color_specular[3];
	double normal[3];
	double shininess;
};

typedef struct _Triangle {
	struct Vertex v[3];
	vec3 surfaceNormal;
} Triangle;

typedef struct _Sphere {
	double position[3];
	double color_diffuse[3];
	double color_specular[3];
	double shininess;
	double radius;
} Sphere;

typedef struct _Light {
	double position[3];
	double color[3];
} Light;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

inline double dotProduct(const vec3& u, const vec3& v) {
	return (u.x * v.x) + (u.y * v.y) + (u.z * v.z);
}

inline void crossProduct(const vec3& u, const vec3& v, vec3& result) {
	result.x = u.y * v.z - u.z * v.y, result.y = u.z * v.x - u.x * v.z, result.z = u.x * v.y - u.y * v.x;
}

inline double getScalar(const vec3& v) {
	return sqrt(pow(v.x, 2) + pow(v.y, 2) + pow(v.z, 2));
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b);

bool triangleIntersection(double x0, double y0, double z0, const vec3& ray, vec3& normal, vec3& surfaceNormal,
						  double ks[3], double kd[3], double& shininess, double& t) {
	bool intersect = false;

	for (int i = 0; i < num_triangles; i++) {
		//ray direction & surface normal
		double NdotD = dotProduct(triangles[i].surfaceNormal, ray);
		if (NdotD != 0) {
			double d = (triangles[i].surfaceNormal.x * triangles[i].v[0].position[0] +
						triangles[i].surfaceNormal.y * triangles[i].v[0].position[1] +
						triangles[i].surfaceNormal.z * triangles[i].v[0].position[2]) * -1;
			double t0 = -(triangles[i].surfaceNormal.x * x0 +
						  triangles[i].surfaceNormal.y * y0 +
						  triangles[i].surfaceNormal.z * z0 + d) / NdotD;

			if (t0 > 0.0001) {		//in front of camera, use 0.0001 instead of 0 to prevent "self block"
				// triangle vertex
				double x1 = triangles[i].v[0].position[0], y1 = triangles[i].v[0].position[1], z1 = triangles[i].v[0].position[2],
					x2 = triangles[i].v[1].position[0], y2 = triangles[i].v[1].position[1], z2 = triangles[i].v[1].position[2],
					x3 = triangles[i].v[2].position[0], y3 = triangles[i].v[2].position[1], z3 = triangles[i].v[2].position[2],
					xp = x0 + t0 * ray.x, yp = y0 + t0 * ray.y, zp = z0 + t0 * ray.z; // intersection point
				double areaTotal = abs((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1));
				double areaP1P2P = 0.0;
				double areaP2P3P = 0.0;
				double areaP1P3P = 0.0;
				if (areaTotal != 0) {	//same x or y
					areaP1P2P = abs((x2 - xp) * (y1 - yp) - (y2 - yp) * (x1 - xp));
					areaP2P3P = abs((x2 - xp) * (y3 - yp) - (y2 - yp) * (x3 - xp));
					areaP1P3P = abs((x1 - xp) * (y3 - yp) - (y1 - yp) * (x3 - xp));
				}
				else {
					areaTotal = abs((x2 - x1) * (z3 - z1) - (z2 - z1) * (x3 - x1));
					areaP1P2P = abs((x2 - xp) * (z1 - zp) - (z2 - zp) * (x1 - xp));
					areaP2P3P = abs((x2 - xp) * (z3 - zp) - (z2 - zp) * (x3 - xp));
					areaP1P3P = abs((x1 - xp) * (z3 - zp) - (z1 - zp) * (x3 - xp));
				}
				//use area difference to do inside test
				if (abs(areaP1P2P + areaP2P3P + areaP1P3P - areaTotal) < 0.000001) {
					if (t > t0) {
						t = t0;
						intersect = true;

						double alpha = areaP2P3P / areaTotal;
						double beta = areaP1P3P / areaTotal;
						double gamma = 1 - alpha - beta;

						normal.x = alpha * triangles[i].v[0].normal[0] + beta * triangles[i].v[1].normal[0] + gamma * triangles[i].v[2].normal[0];
						normal.y = alpha * triangles[i].v[0].normal[1] + beta * triangles[i].v[1].normal[1] + gamma * triangles[i].v[2].normal[1];
						normal.z = alpha * triangles[i].v[0].normal[2] + beta * triangles[i].v[1].normal[2] + gamma * triangles[i].v[2].normal[2];
						normal.normalize();

						surfaceNormal.x = triangles[i].surfaceNormal.x;
						surfaceNormal.y = triangles[i].surfaceNormal.y;
						surfaceNormal.z = triangles[i].surfaceNormal.z;

						kd[0] = alpha * triangles[i].v[0].color_diffuse[0] + beta * triangles[i].v[1].color_diffuse[0] + gamma * triangles[i].v[2].color_diffuse[0];
						kd[1] = alpha * triangles[i].v[0].color_diffuse[1] + beta * triangles[i].v[1].color_diffuse[1] + gamma * triangles[i].v[2].color_diffuse[1];
						kd[2] = alpha * triangles[i].v[0].color_diffuse[2] + beta * triangles[i].v[1].color_diffuse[2] + gamma * triangles[i].v[2].color_diffuse[2];
						ks[0] = alpha * triangles[i].v[0].color_specular[0] + beta * triangles[i].v[1].color_specular[0] + gamma * triangles[i].v[2].color_specular[0];
						ks[1] = alpha * triangles[i].v[0].color_specular[1] + beta * triangles[i].v[1].color_specular[1] + gamma * triangles[i].v[2].color_specular[1];
						ks[2] = alpha * triangles[i].v[0].color_specular[2] + beta * triangles[i].v[1].color_specular[2] + gamma * triangles[i].v[2].color_specular[2];

						shininess = alpha * triangles[i].v[0].shininess + beta * triangles[i].v[1].shininess + gamma * triangles[i].v[2].shininess;
					}
				}
			}
		}
	}

	return intersect;
}

bool sphereIntersection(double x0, double y0, double z0, const vec3& ray, vec3& normal, vec3& surfaceNormal,
						double ks[3], double kd[3], double& shininess, double& t) {
	bool intersect = false;

	for (int i = 0; i < num_spheres; i++) {
		double b = 2 * (ray.x * (x0 - spheres[i].position[0]) +
						ray.y * (y0 - spheres[i].position[1]) +
						ray.z * (z0 - spheres[i].position[2]));
		double c = pow(x0 - spheres[i].position[0], 2) +
			pow(y0 - spheres[i].position[1], 2) +
			pow(z0 - spheres[i].position[2], 2) - pow(spheres[i].radius, 2);
		//intersection
		if (pow(b, 2) - 4 * c > 0) {
			double t0 = (-b - sqrt(pow(b, 2) - 4 * c)) / 2, t1 = (-b + sqrt(pow(b, 2) - 4 * c)) / 2;
			if (t0 >= 0.0001 && t1 >= 0.0001) { //use 0.001 instead of 0 to prevent "self block"
				if (t > std::min(t0, t1)) {
					t = std::min(t0, t1);
					intersect = true;

					normal.x = x0 + t * ray.x - spheres[i].position[0];
					normal.y = y0 + t * ray.y - spheres[i].position[1];
					normal.z = z0 + t * ray.z - spheres[i].position[2];
					normal.normalize();

					surfaceNormal.x = normal.x;
					surfaceNormal.y = normal.y;
					surfaceNormal.z = normal.z;

					kd[0] = spheres[i].color_diffuse[0], kd[1] = spheres[i].color_diffuse[1], kd[2] = spheres[i].color_diffuse[2];
					ks[0] = spheres[i].color_specular[0], ks[1] = spheres[i].color_specular[1], ks[2] = spheres[i].color_specular[2];
					shininess = spheres[i].shininess;
				}
			}
		}
	}

	return intersect;
}

bool findIntersection(double x0, double y0, double z0, const vec3& ray, vec3& normal, vec3& surfaceNormal,
					  double ks[3], double kd[3], double& shininess, double& t) {
	bool triangle = triangleIntersection(x0, y0, z0, ray, normal, surfaceNormal, ks, kd, shininess, t);
	bool sphere = sphereIntersection(x0, y0, z0, ray, normal, surfaceNormal, ks, kd, shininess, t);

	return triangle || sphere;
}

//MODIFY THIS FUNCTION
void draw_scene() {
	glPointSize(2.0);


	//for every pixel
	for (unsigned int x = 0; x < WIDTH; x++) {
		glBegin(GL_POINTS);
		for (unsigned int y = 0; y < HEIGHT; y++) {
			//anti aliasing
			double pixelR = 0.0, pixelG = 0.0, pixelB = 0.0;
			for (int m = 0; m < 2; m++) {
				for (int n = 0; n < 2; n++) {
					//pixel coordinates
					double x0 = -xMax + 2 * xMax * (static_cast<float>(x * 2 + m) / static_cast<float>(WIDTH * 2));
					double y0 = -yMax + 2 * yMax * (static_cast<float>(y * 2 + n) / static_cast<float>(HEIGHT * 2));
					//double x0 = -xMax + 2 * xMax * (static_cast<float>(x) / static_cast<float>(WIDTH));
					//double y0 = -yMax + 2 * yMax * (static_cast<float>(y) / static_cast<float>(HEIGHT));

					// colors
					double r0 = 0.0, g0 = 0.0, b0 = 0.0;
					double localR = 1.0, localG = 1.0, localB = 1.0;

					// camera position
					double startX = 0.0, startY = 0.0, startZ = 0.0;
					//generate ray vector
					vec3 ray(x0, y0, -1.0);
					ray.normalize();

					double t = 1000.0;
					vec3 normal;
					vec3 surfaceNormal;
					double kd[3] = { 0.0,0.0,0.0 };
					double ks[3] = { 0.0,0.0,0.0 };
					double prevKs[3] = { 1.0,1.0,1.0 };
					double shininess = 0.0;

					int reflection = 10;		// how many reflections
					// while intersect
					while (reflection && findIntersection(startX, startY, startZ, ray, normal, surfaceNormal, ks, kd, shininess, t)) {
						--reflection;
						//ambient light
						localR = ambient_light[0], localG = ambient_light[1], localB = ambient_light[2];

						//diffuse & specular
						for (int j = 0; j < num_lights; j++) {
							//soft shadows
							vec3 lightVector;
							//point light becomes lightArea * lightArea area light, each point intensity becomes 1/lightArea * lightArea
							int lightArea = 1, lightIntensity = pow(lightArea, 2);
							for (int k = -lightArea / 2; k < (lightArea + 1) / 2; k++) {
								for (int l = -lightArea / 2; l < (lightArea + 1) / 2; l++) {
									lightVector.x = lights[j].position[0] + 0.008 * k - (startX + t * ray.x);
									lightVector.y = lights[j].position[1] + 0.008 * l - (startY + t * ray.y);
									lightVector.z = lights[j].position[2] - (startZ + t * ray.z);
									double scalarToLight = getScalar(lightVector);
									lightVector.normalize();

									vec3 normal2;
									vec3 surfaceNormal2;
									vec3 toLight(lightVector);
									double kd2[3] = { 0.0,0.0,0.0 };
									double ks2[3] = { 0.0,0.0,0.0 };
									double shininess2 = 0.0;
									double t2 = 1000.0;

									if (!(findIntersection(startX + t * ray.x, startY + t * ray.y, startZ + t * ray.z, lightVector, normal2, surfaceNormal2, ks2, kd2, shininess2, t2) &&
										  t2 < scalarToLight)) {

										double NdotL = dotProduct(normal, lightVector);
										if (NdotL > 0) {
											localR += kd[0] * lights[j].color[0] / lightIntensity * NdotL;
											localG += kd[1] * lights[j].color[1] / lightIntensity * NdotL;
											localB += kd[2] * lights[j].color[2] / lightIntensity * NdotL;
										}

										//reflected vector
										vec3 reflectVector(2 * NdotL * normal.x - lightVector.x,
														   2 * NdotL * normal.y - lightVector.y,
														   2 * NdotL * normal.z - lightVector.z);
										reflectVector.normalize();
										double VdotR = -dotProduct(ray, reflectVector);
										if (VdotR > 0) {
											localR += ks[0] * lights[j].color[0] / lightIntensity * pow(VdotR, shininess);
											localG += ks[1] * lights[j].color[1] / lightIntensity * pow(VdotR, shininess);
											localB += ks[2] * lights[j].color[2] / lightIntensity * pow(VdotR, shininess);
										}
									}
								}
							}
						}

						if (localR > 1.0) { localR = 1.0; }
						if (localG > 1.0) { localG = 1.0; }
						if (localB > 1.0) { localB = 1.0; }

						r0 += localR * (1 - ks[0]) * prevKs[0], g0 += localG * (1 - ks[1]) * prevKs[1], b0 += localB * (1 - ks[2]) * prevKs[2];
						//r0 += localR , g0 += localG , b0 += localB ;
						prevKs[0] *= ks[0], prevKs[1] *= ks[1], prevKs[2] *= ks[2];
						if (std::min(prevKs[0], std::min(prevKs[1], prevKs[2])) < 0.01)
							break;

						//if white
						if (r0 >= 1.0 && g0 >= 1.0 && b0 >= 1.0)
							break;

						startX += t * ray.x, startY += t * ray.y, startZ += t * ray.z;
						//reflect direction
						double NdotV = -dotProduct(surfaceNormal, ray);
						ray.x = 2 * NdotV * surfaceNormal.x + ray.x;
						ray.y = 2 * NdotV * surfaceNormal.y + ray.y;
						ray.z = 2 * NdotV * surfaceNormal.z + ray.z;
						ray.normalize();

						t = 999.0;
					}

					if (t == 1000.0) {
						//r0 = 0.3, g0 = 0.3, b0 = 0.3;
						r0 = 1.0, g0 = 1.0, b0 = 1.0;
					}
					else {
						//ambient reflection 
						r0 += ambient_light[0] * prevKs[0], g0 += ambient_light[1] * prevKs[1], b0 += ambient_light[2] * prevKs[2];

						if (r0 > 1) { r0 = 1; }
						if (g0 > 1) { g0 = 1; }
						if (b0 > 1) { b0 = 1; }
					}

					pixelR += r0;
					pixelG += g0;
					pixelB += b0;
				}
			}

			unsigned char r = pixelR / 4 * 255;
			unsigned char g = pixelG / 4 * 255;
			unsigned char b = pixelB / 4 * 255;

			plot_pixel(x, y, r, g, b);
		}

		glEnd();
		glFlush();
	}

	printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b) {
	glColor3f(((double)r) / 256.f, ((double)g) / 256.f, ((double)b) / 256.f);
	glVertex2i(x, y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b) {
	buffer[HEIGHT - y - 1][x][0] = r;
	buffer[HEIGHT - y - 1][x][1] = g;
	buffer[HEIGHT - y - 1][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b) {
	plot_pixel_display(x, y, r, g, b);
	if (mode == PLOTMODE::MODE_JPEG)
		plot_pixel_jpeg(x, y, r, g, b);
}

/* Write a jpg image from buffer*/
void save_jpg() {
	if (filename == NULL)
		return;

	// Allocate a picture buffer // 
	cv::Mat3b bufferBGR = cv::Mat::zeros(HEIGHT, WIDTH, CV_8UC3); //rows, cols, 3-channel 8-bit.
	printf("File to save to: %s\n", filename);

	// unsigned char buffer[HEIGHT][WIDTH][3];
	for (int r = 0; r < HEIGHT; r++) {
		for (int c = 0; c < WIDTH; c++) {
			for (int chan = 0; chan < 3; chan++) {
				unsigned char red = buffer[r][c][0];
				unsigned char green = buffer[r][c][1];
				unsigned char blue = buffer[r][c][2];
				bufferBGR.at<cv::Vec3b>(r, c) = cv::Vec3b(blue, green, red);
			}
		}
	}
	if (cv::imwrite(filename, bufferBGR)) {
		printf("File saved Successfully\n");
	}
	else {
		printf("Error in Saving\n");
	}
}

void parse_check(char* expected, char* found) {
	if (stricmp(expected, found)) {
		char error[100];
		printf("Expected '%s ' found '%s '\n", expected, found);
		printf("Parse error, abnormal abortion\n");
		exit(0);
	}

}

void parse_doubles(FILE* file, char* check, double p[3]) {
	char str[100];
	fscanf(file, "%s", str);
	parse_check(check, str);
	fscanf(file, "%lf %lf %lf", &p[0], &p[1], &p[2]);
	printf("%s %lf %lf %lf\n", check, p[0], p[1], p[2]);
}

void parse_rad(FILE* file, double* r) {
	char str[100];
	fscanf(file, "%s", str);
	parse_check("rad:", str);
	fscanf(file, "%lf", r);
	printf("rad: %f\n", *r);
}

void parse_shi(FILE* file, double* shi) {
	char s[100];
	fscanf(file, "%s", s);
	parse_check("shi:", s);
	fscanf(file, "%lf", shi);
	printf("shi: %f\n", *shi);
}

int loadScene(char* argv) {
	FILE* file = fopen(argv, "r");
	int number_of_objects;
	char type[50];
	int i;
	Triangle t;
	Sphere s;
	Light l;
	fscanf(file, "%i", &number_of_objects);

	printf("number of objects: %i\n", number_of_objects);
	char str[200];

	parse_doubles(file, "amb:", ambient_light);

	for (i = 0; i < number_of_objects; i++) {
		fscanf(file, "%s\n", type);
		printf("%s\n", type);
		if (stricmp(type, "triangle") == 0) {

			printf("found triangle\n");
			int j;

			for (j = 0; j < 3; j++) {
				parse_doubles(file, "pos:", t.v[j].position);
				parse_doubles(file, "nor:", t.v[j].normal);
				parse_doubles(file, "dif:", t.v[j].color_diffuse);
				parse_doubles(file, "spe:", t.v[j].color_specular);
				parse_shi(file, &t.v[j].shininess);
			}

			if (num_triangles == MAX_TRIANGLES) {
				printf("too many triangles, you should increase MAX_TRIANGLES!\n");
				exit(0);
			}
			triangles[num_triangles++] = t;
		}
		else if (stricmp(type, "sphere") == 0) {
			printf("found sphere\n");

			parse_doubles(file, "pos:", s.position);
			parse_rad(file, &s.radius);
			parse_doubles(file, "dif:", s.color_diffuse);
			parse_doubles(file, "spe:", s.color_specular);
			parse_shi(file, &s.shininess);

			if (num_spheres == MAX_SPHERES) {
				printf("too many spheres, you should increase MAX_SPHERES!\n");
				exit(0);
			}
			spheres[num_spheres++] = s;
		}
		else if (stricmp(type, "light") == 0) {
			printf("found light\n");
			parse_doubles(file, "pos:", l.position);
			parse_doubles(file, "col:", l.color);

			if (num_lights == MAX_LIGHTS) {
				printf("too many lights, you should increase MAX_LIGHTS!\n");
				exit(0);
			}
			lights[num_lights++] = l;
		}
		else {
			printf("unknown type in scene description:\n%s\n", type);
			exit(0);
		}
	}

	//triangle normal
	for (int i = 0; i < num_triangles; i++) {
		//get surface normal
		vec3 u(triangles[i].v[1].position[0] - triangles[i].v[0].position[0],
			   triangles[i].v[1].position[1] - triangles[i].v[0].position[1],
			   triangles[i].v[1].position[2] - triangles[i].v[0].position[2]);
		vec3 v(triangles[i].v[2].position[0] - triangles[i].v[0].position[0],
			   triangles[i].v[2].position[1] - triangles[i].v[0].position[1],
			   triangles[i].v[2].position[2] - triangles[i].v[0].position[2]);
		crossProduct(u, v, triangles[i].surfaceNormal);
		triangles[i].surfaceNormal.normalize();
	}

	return 0;
}

void display() {}

void init() {
	glMatrixMode(GL_PROJECTION);
	glOrtho(0, WIDTH, 0, HEIGHT, 1, -1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glClearColor(0, 0, 0, 0);
	glClear(GL_COLOR_BUFFER_BIT);
}

void idle() {
	//hack to make it only draw once
	static int once = 0;
	if (!once) {
		draw_scene();
		if (mode == PLOTMODE::MODE_JPEG)
			save_jpg();
	}
	once = 1;
}

int main(int argc, char** argv) {
	if (argc < 2 || argc > 3) {
		printf("usage: %s <scenefile> [jpegname]\n", argv[0]);
		exit(0);
	}
	if (argc == 3) {
		mode = PLOTMODE::MODE_JPEG;
		filename = argv[2];
	}
	else if (argc == 2)
		mode = PLOTMODE::MODE_DISPLAY;

	glutInit(&argc, argv);
	printf("loadScene: %s\n\n", argv[1]);
	loadScene(argv[1]);

	glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
	glutInitWindowPosition(300, 100);
	glutInitWindowSize(WIDTH, HEIGHT);
	int window = glutCreateWindow("Ray Tracer");
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	init();
	glutMainLoop();
}
