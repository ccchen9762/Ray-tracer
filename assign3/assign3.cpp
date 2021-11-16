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
//const int WIDTH = 640;
//const int HEIGHT = 480;
const int WIDTH = 320;
const int HEIGHT = 240;
const float ASPECT_RATIO = static_cast<float>(WIDTH) / static_cast<float>(HEIGHT);

//the field of view of the camera
const double fov = 60.0;
const double Pi = 3.1415926;
const double fovPi = (fov / 360.0) * Pi;

//boundary coordinates
double yMax = tan(fovPi / 2);
double xMax = ASPECT_RATIO * yMax;

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

double dotProduct(const vec3& v0, const vec3& v1) {
	return v0.x * v1.x + v0.y + v1.y + v0.z + v1.z;
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b);

//MODIFY THIS FUNCTION
void draw_scene() {
	glPointSize(2.0);
	glBegin(GL_POINTS);
	for (unsigned int x = 0; x < WIDTH; x++) {
		for (unsigned int y = 0; y < HEIGHT; y++) {
			//pixel coordinates
			double x0 = -xMax + 2 * xMax * (static_cast<float>(x) / static_cast<float>(WIDTH));
			double y0 = -yMax + 2 * yMax * (static_cast<float>(y) / static_cast<float>(HEIGHT));
			//double x0 = -ASPECT_RATIO * tan(fovPi / 2) +(2 * ASPECT_RATIO * tan(fovPi / 2)) * (static_cast<float>(x) / static_cast<float>(WIDTH));
			//double y0 = -tan(fovPi / 2) +(2 * tan(fovPi / 2)) * (static_cast<float>(y) / static_cast<float>(HEIGHT));
			
			//ambient light
			unsigned char r = ambient_light[0] * 255;
			unsigned char g = ambient_light[1] * 255;
			unsigned char b = ambient_light[2] * 255;

			//generate vector
			vec3 ray(x0, y0, -1.0);
			ray.normalize();
			//printf("%f, %f, %f \n", ray.x, ray.y, ray.z);
			
			//check intersections
			for (int i = 0; i < num_spheres; i++) {
				double b = 2 * (ray.x * -spheres[i].position[0] + ray.y * -spheres[i].position[1] + ray.z * -spheres[i].position[2]);
				double c = pow(-spheres[i].position[0], 2) + pow(-spheres[i].position[1], 2) + pow(-spheres[i].position[2], 2) - pow(spheres[i].radius, 2);
				//intersection
				if (pow(b, 2) - 4 * c > 0) {
					//printf("%f, %f\n", b, c);
					double t0 = (-b - sqrt(pow(b, 2) - 4 * c)) / 2, t1 = (-b + sqrt(pow(b, 2) - 4 * c)) / 2;
					if (t0 >= 0 && t1 >= 0) {
						double t = std::min(t0, t1);
						vec3 normal(t * ray.x- spheres[i].position[0], t * ray.y - spheres[i].position[1], t * ray.z - spheres[i].position[2]);
						normal.normalize();
						for (int j = 0; j < num_lights; j++) {
							vec3 lightVector(t * ray.x, t * ray.y, t * ray.z);
							double NdotL = dotProduct(normal, lightVector);
							if (NdotL > 0) {
								r += spheres[i].color_diffuse[0] * lights[j].color[0] * NdotL;
								g += spheres[i].color_diffuse[1] * lights[j].color[1] * NdotL;
								b += spheres[i].color_diffuse[2] * lights[j].color[2] * NdotL;
							}
						}
					}
				}
			}
			for (int i = 0; i < num_triangles; i++) {

			}
			plot_pixel(x, y, r, g, b);
		}
	}
	glEnd();
	glFlush();
	//simple output
	/*for (x = 0; x < WIDTH; x++) {
		glPointSize(2.0);
		glBegin(GL_POINTS);
		for (y = 0; y < HEIGHT; y++) {
			plot_pixel(x, y, x % 256, y % 256, (x + y) % 256);
		}
		glEnd();
		glFlush();
	}*/
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
	glutInitWindowPosition(500, 200);
	glutInitWindowSize(WIDTH, HEIGHT);
	int window = glutCreateWindow("Ray Tracer");
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	init();
	glutMainLoop();
}
