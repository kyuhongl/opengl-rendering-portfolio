/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: kyuhong lee
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>
#include <math.h>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

// The different display modes.
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

// While solving the homework, it is useful to make the below values smaller for debugging purposes.
// The still images that you need to submit with the homework should be at the below resolution (640x480).
// However, for your own purposes, after you have solved the homework, you can increase those values to obtain higher-resolution images.
#define WIDTH 640
#define HEIGHT 480

// The field of view of the camera, in degrees.
#define fov 60.0

// Buffer to store the image when saving it to a JPEG.
unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
  double reflectivity;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
  double reflectivity;
};

struct Light
{
  double position[3];
  double color[3];
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b);

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x, HEIGHT - 1 - y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  // Flip the y-coordinate for the JPEG output
  int flipped_y = HEIGHT - 1 - y;
  buffer[flipped_y][x][0] = r;
  buffer[flipped_y][x][1] = g;
  buffer[flipped_y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x, y, r, g, b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x, y, r, g, b);
}

// Ray structure to represent rays in the scene
struct Ray {
  double origin[3];      // Origin of the ray
  double direction[3];   // Direction of the ray (should be normalized)
};

// Hit info to store information about ray-object intersections
struct HitInfo {
  bool hit;              // Whether the ray hit something
  double t;              // Parameter value at intersection point
  double point[3];       // 3D coordinates of intersection point
  double normal[3];      // Surface normal at intersection point
  double color_diffuse[3];  // Diffuse color at intersection point
  double color_specular[3]; // Specular color at intersection point
  double shininess;      // Shininess at intersection point
};

// Vector operations
void normalize(double v[3]) {
  double len = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

void cross(const double a[3], const double b[3], double result[3]) {
  result[0] = a[1]*b[2] - a[2]*b[1];
  result[1] = a[2]*b[0] - a[0]*b[2];
  result[2] = a[0]*b[1] - a[1]*b[0];
}

double dot(const double a[3], const double b[3]) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void subtract(const double a[3], const double b[3], double result[3]) {
  result[0] = a[0] - b[0];
  result[1] = a[1] - b[1];
  result[2] = a[2] - b[2];
}

void add(const double a[3], const double b[3], double result[3]) {
  result[0] = a[0] + b[0];
  result[1] = a[1] + b[1];
  result[2] = a[2] + b[2];
}

void multiply(const double a[3], double scalar, double result[3]) {
  result[0] = a[0] * scalar;
  result[1] = a[1] * scalar;
  result[2] = a[2] * scalar;
}

// Returns the squared distance between two points
double distance_squared(const double a[3], const double b[3]) {
  double dx = a[0] - b[0];
  double dy = a[1] - b[1];
  double dz = a[2] - b[2];
  return dx*dx + dy*dy + dz*dz;
}

// Generate a ray from the camera through pixel (x,y)
void generate_ray(int x, int y, Ray &ray) {
  // Camera is at origin (0,0,0)
  ray.origin[0] = 0;
  ray.origin[1] = 0;
  ray.origin[2] = 0;
  
  // Calculate the direction of the ray using the field of view
  double aspect = (double)WIDTH / (double)HEIGHT;
  double fov_rad = fov * M_PI / 180.0;
  double tan_fov = tan(fov_rad / 2.0);
  
  // Map pixel coordinates to [-1,1]
  double pixel_x = (2.0 * ((x + 0.5) / WIDTH) - 1) * aspect * tan_fov;
  double pixel_y = (1.0 - 2.0 * ((y + 0.5) / HEIGHT)) * tan_fov;
  
  // The ray direction points from the origin through the pixel on the image plane
  ray.direction[0] = pixel_x;
  ray.direction[1] = pixel_y;
  ray.direction[2] = -1.0;  // Negative z-direction (forward)
  
  // Normalize the direction
  normalize(ray.direction);
}

// =============================
// triangle intersection 
// =============================
bool intersect_triangle(const Ray &ray, const Triangle &triangle, HitInfo &hit_info) {
  double v0[3], v1[3], v2[3];
  for (int i = 0; i < 3; i++) {
    v0[i] = triangle.v[0].position[i];
    v1[i] = triangle.v[1].position[i];
    v2[i] = triangle.v[2].position[i];
  }
  
  double edge1[3], edge2[3];
  subtract(v1, v0, edge1);
  subtract(v2, v0, edge2);
  
  double pvec[3];
  cross(ray.direction, edge2, pvec);
  
  double det = dot(edge1, pvec);
  
  // Ray is parallel to the triangle or triangle is degenerate
  if (fabs(det) < 1e-8) {
    return false;
  }
  
  double inv_det = 1.0 / det;
  
  double tvec[3];
  subtract(ray.origin, v0, tvec);
  
  double u = dot(tvec, pvec) * inv_det;
  
  // Intersection outside of triangle
  if (u < 0.0 || u > 1.0) {
    return false;
  }
  
  double qvec[3];
  cross(tvec, edge1, qvec);
  
  double v = dot(ray.direction, qvec) * inv_det;
  
  // Intersection outside of triangle
  if (v < 0.0 || u + v > 1.0) {
    return false;
  }
  
  // Calculate the distance to the intersection point
  double t = dot(edge2, qvec) * inv_det;
  
  // Only consider intersections in front of the ray
  if (t < 0) {
    return false;
  }
  
  hit_info.hit = true;
  hit_info.t = t;
  
  // Calculate intersection point
  multiply(ray.direction, t, hit_info.point);
  add(ray.origin, hit_info.point, hit_info.point);
  
  // Calculate barycentric coordinates
  double w = 1.0 - u - v;
  
  // Interpolate the normal
  for (int i = 0; i < 3; i++) {
    hit_info.normal[i] = w * triangle.v[0].normal[i] + 
                          u * triangle.v[1].normal[i] + 
                          v * triangle.v[2].normal[i];
    
    // Interpolate diffuse color
    hit_info.color_diffuse[i] = w * triangle.v[0].color_diffuse[i] + 
                                 u * triangle.v[1].color_diffuse[i] + 
                                 v * triangle.v[2].color_diffuse[i];
    
    // Interpolate specular color
    hit_info.color_specular[i] = w * triangle.v[0].color_specular[i] + 
                                  u * triangle.v[1].color_specular[i] + 
                                  v * triangle.v[2].color_specular[i];
  }
  
  // Interpolate shininess
  hit_info.shininess = w * triangle.v[0].shininess + 
                       u * triangle.v[1].shininess + 
                       v * triangle.v[2].shininess;
  
  // Normalize the interpolated normal
  normalize(hit_info.normal);
  
  return true;
}

// =============================
// sphere intersection 
// =============================
bool intersect_sphere(const Ray &ray, const Sphere &sphere, HitInfo &hit_info) {
  double oc[3];
  subtract(ray.origin, sphere.position, oc);
  
  double a = dot(ray.direction, ray.direction);
  double b = 2.0 * dot(oc, ray.direction);
  double c = dot(oc, oc) - sphere.radius * sphere.radius;
  
  double discriminant = b*b - 4*a*c;
  
  if (discriminant < 0) {
    return false;  // No intersection
  }
  
  // Find the nearest intersection point
  double t = (-b - sqrt(discriminant)) / (2.0 * a);
  
  // Check if the intersection is behind the ray origin
  if (t < 0) {
    t = (-b + sqrt(discriminant)) / (2.0 * a);
    if (t < 0) {
      return false;  // Both intersection points are behind the ray
    }
  }
  
  // Calculate the intersection point
  hit_info.hit = true;
  hit_info.t = t;
  
  // Compute the intersection point
  multiply(ray.direction, t, hit_info.point);
  add(ray.origin, hit_info.point, hit_info.point);
  
  // Compute the normal at the intersection
  subtract(hit_info.point, sphere.position, hit_info.normal);
  normalize(hit_info.normal);
  
  // Store material properties
  for (int i = 0; i < 3; i++) {
    hit_info.color_diffuse[i] = sphere.color_diffuse[i];
    hit_info.color_specular[i] = sphere.color_specular[i];
  }
  hit_info.shininess = sphere.shininess;
  
  return true;
}

// =============================
// find the closest intersection in the scene
// =============================
bool trace_ray(const Ray &ray, HitInfo &hit_info, double min_t = 0.0001, double max_t = 1e10) {
  hit_info.hit = false;
  hit_info.t = max_t;

  // Check intersections with all spheres
  for (int i = 0; i < num_spheres; i++) {
    HitInfo temp_hit;
    if (intersect_sphere(ray, spheres[i], temp_hit) && temp_hit.t > min_t && temp_hit.t < hit_info.t) {
      hit_info = temp_hit;
    }
  }

  // Check intersections with all triangles
  for (int i = 0; i < num_triangles; i++) {
    HitInfo temp_hit;
    if (intersect_triangle(ray, triangles[i], temp_hit) && temp_hit.t > min_t && temp_hit.t < hit_info.t) {
      hit_info = temp_hit;
    }
  }

  return hit_info.hit;
}

// =============================
// triangle phong shading 
// =============================
// handled in shade() when hit_info comes from a triangle

// =============================
// sphere phong shading 
// =============================
// handled in shade() when hit_info comes from a sphere

// =============================
// shadow rays 
// =============================
// handled in is_in_shadow() and shade() for direct lighting
bool is_in_shadow(const double point[3], const Light &light) {
  // Soft shadows: Sample multiple points around the light
  const int shadow_samples = 8;
  int shadow_count = 0;
  
  for (int s = 0; s < shadow_samples; s++) {
    Ray shadow_ray;
    
    // Shadow ray origin is the intersection point
    for (int i = 0; i < 3; i++) {
      shadow_ray.origin[i] = point[i];
    }
    
    // Calculate a sample position around the light source
    double light_pos_jittered[3];
    for (int i = 0; i < 3; i++) {
      light_pos_jittered[i] = light.position[i];
    }
    
    // Only apply jittering if not the first sample
    if (s > 0) {
      // Create a jittered sample in a disk perpendicular to the main direction
      // First determine the main direction from point to light
      double main_dir[3];
      subtract(light.position, point, main_dir);
      normalize(main_dir);
      
      // Find an arbitrary perpendicular vector
      double perp1[3] = {0, 0, 0};
      if (fabs(main_dir[0]) < fabs(main_dir[1])) {
        if (fabs(main_dir[0]) < fabs(main_dir[2])) {
          perp1[0] = 1.0;
        } else {
          perp1[2] = 1.0;
        }
      } else {
        if (fabs(main_dir[1]) < fabs(main_dir[2])) {
          perp1[1] = 1.0;
        } else {
          perp1[2] = 1.0;
        }
      }
      
      // Cross to get the first perpendicular vector
      double perp2[3];
      cross(main_dir, perp1, perp2);
      normalize(perp2);
      
      // And cross again to get the second perpendicular vector
      cross(main_dir, perp2, perp1);
      normalize(perp1);
      
      // Now create a random point in a disk of radius 0.2
      double radius = 0.2 * sqrt((double)rand() / RAND_MAX);
      double angle = 2.0 * M_PI * ((double)rand() / RAND_MAX);
      
      // Offset the light position
      double offset[3];
      for (int i = 0; i < 3; i++) {
        offset[i] = radius * (cos(angle) * perp1[i] + sin(angle) * perp2[i]);
        light_pos_jittered[i] += offset[i];
      }
    }
    
    // Direction is from the point to the jittered light position
    subtract(light_pos_jittered, point, shadow_ray.direction);
    double light_distance = sqrt(dot(shadow_ray.direction, shadow_ray.direction));
    normalize(shadow_ray.direction);
    
    // Check if any object blocks the light
    HitInfo shadow_hit;
    if (trace_ray(shadow_ray, shadow_hit, 0.0001, light_distance)) {
      shadow_count++;
    }
  }
  
  // Return true if majority of shadow rays hit an object
  return shadow_count > shadow_samples / 2;
}

// =============================
// recursive reflection (ray bounces)
// =============================
void shade(const Ray &ray, const HitInfo &hit_info, double color[3], int depth = 0) {
  // Set initial color to ambient light
  for (int i = 0; i < 3; i++) {
    color[i] = ambient_light[i] * hit_info.color_diffuse[i];
  }
  
  // For each light source
  for (int light_idx = 0; light_idx < num_lights; light_idx++) {
    // Skip if point is in shadow for this light
    if (is_in_shadow(hit_info.point, lights[light_idx])) {
      continue;
    }
    
    // Light direction (L)
    double light_dir[3];
    subtract(lights[light_idx].position, hit_info.point, light_dir);
    normalize(light_dir);
    
    // View direction (V)
    double view_dir[3];
    multiply(ray.direction, -1.0, view_dir);
    normalize(view_dir);
    
    // Reflection direction (R)
    double reflect_dir[3];
    double dot_ln = dot(light_dir, hit_info.normal);
    double temp[3];
    multiply(hit_info.normal, 2.0 * dot_ln, temp);
    subtract(temp, light_dir, reflect_dir);
    normalize(reflect_dir);
    
    // Diffuse term: kd * (L dot N)
    double diffuse = fmax(0.0, dot_ln);
    
    // Specular term: ks * (R dot V)^sh
    double specular = pow(fmax(0.0, dot(reflect_dir, view_dir)), hit_info.shininess);
    
    // Add contribution from this light
    for (int i = 0; i < 3; i++) {
      color[i] += lights[light_idx].color[i] * (
        hit_info.color_diffuse[i] * diffuse +
        hit_info.color_specular[i] * specular
      );
    }
  }
  
  // Calculate reflection color if not at max recursion depth
  if (depth < 5) {
    // Calculate reflection direction
    double refl_dir[3];
    double dot_vn = dot(ray.direction, hit_info.normal);
    double temp[3];
    multiply(hit_info.normal, 2.0 * dot_vn, temp);
    subtract(ray.direction, temp, refl_dir);
    normalize(refl_dir);
    
    // Create reflection ray
    Ray refl_ray;
    for (int i = 0; i < 3; i++) {
      refl_ray.origin[i] = hit_info.point[i];
      refl_ray.direction[i] = refl_dir[i];
    }
    
    // Trace reflection ray
    HitInfo refl_hit;
    if (trace_ray(refl_ray, refl_hit, 0.0001)) {
      double refl_color[3] = {0.0, 0.0, 0.0};
      shade(refl_ray, refl_hit, refl_color, depth + 1);
      
      // use reflectivity for blending reflection color with local color
      double reflectivity = 0.0;
      if (hit_info.shininess > 0.0) {
        reflectivity = hit_info.shininess / 200.0; // fallback if not set
      }
      double ks_avg = (hit_info.color_specular[0] + hit_info.color_specular[1] + hit_info.color_specular[2]) / 3.0;
      if (hit_info.shininess > 0.0 && hit_info.t > 0.0) {
        reflectivity = ks_avg; // fallback if reflectivity not set
      }
      // blend using reflectivity
      for (int i = 0; i < 3; i++) {
        color[i] = (1.0 - reflectivity) * color[i] + reflectivity * refl_color[i];
      }
    }
  }
  
  // Clamp color values to [0, 1]
  for (int i = 0; i < 3; i++) {
    if (color[i] > 1.0) color[i] = 1.0;
    if (color[i] < 0.0) color[i] = 0.0;
  }
}

// Main ray tracing function
void draw_scene() {
  for(unsigned int y = 0; y < HEIGHT; y++) {
    glPointSize(2.0);
    glBegin(GL_POINTS);
    
    for(unsigned int x = 0; x < WIDTH; x++) {
      // Default color (background - set to black)
      double color[3] = {0.0, 0.0, 0.0};
      
      // use stratified sampling for anti-aliasing to reduce noise and match reference structure
      const int sqrt_samples = 4;
      const int samples = sqrt_samples * sqrt_samples;
      for (int p = 0; p < sqrt_samples; ++p) {
        for (int q = 0; q < sqrt_samples; ++q) {
          double jitter_x = (p + ((double)rand() / RAND_MAX)) / sqrt_samples;
          double jitter_y = (q + ((double)rand() / RAND_MAX)) / sqrt_samples;

          Ray ray;
          ray.origin[0] = 0;
          ray.origin[1] = 0;
          ray.origin[2] = 0;

          double aspect = (double)WIDTH / (double)HEIGHT;
          double fov_rad = fov * M_PI / 180.0;
          double tan_fov = tan(fov_rad / 2.0);

          double pixel_x = (2.0 * ((x + jitter_x) / WIDTH) - 1) * aspect * tan_fov;
          double pixel_y = (1.0 - 2.0 * ((y + jitter_y) / HEIGHT)) * tan_fov;

          ray.direction[0] = pixel_x;
          ray.direction[1] = pixel_y;
          ray.direction[2] = -1.0;
          normalize(ray.direction);

          HitInfo hit_info;
          double sample_color[3] = {0.0, 0.0, 0.0};

          if (trace_ray(ray, hit_info)) {
            shade(ray, hit_info, sample_color, 0);
          } else {
            sample_color[0] = 1.0;
            sample_color[1] = 1.0;
            sample_color[2] = 1.0;
          }

          for (int i = 0; i < 3; i++) {
            color[i] += sample_color[i] / samples;
          }
        }
      }
      
      // Convert to unsigned char for plotting
      unsigned char r = (unsigned char)(color[0] * 255);
      unsigned char g = (unsigned char)(color[1] * 255);
      unsigned char b = (unsigned char)(color[2] * 255);
      
      // Plot the pixel
      plot_pixel(x, y, r, g, b);
    }
    
    glEnd();
    glFlush();
  }
  
  printf("Ray tracing completed.\n");
  fflush(stdout);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in saving\n");
  else 
    printf("File saved successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parsing error; abnormal program abortion.\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

void parse_reflectivity(FILE* file, double* reflectivity) {
  char s[100];
  int pos = ftell(file);
  if (fscanf(file, "%s", s) == 1 && strcasecmp("ref:", s) == 0) {
    fscanf(file, "%lf", reflectivity);
    printf("ref: %f\n", *reflectivity);
  } else {
    *reflectivity = 0.0;
    fseek(file, pos, SEEK_SET);
  }
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  if (!file)
  {
    printf("Unable to open input file %s. Program exiting.\n", argv);
    exit(0);
  }

  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
        parse_reflectivity(file, &t.v[j].reflectivity);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);
      parse_reflectivity(file, &s.reflectivity);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;
    
  // Initialize random number generator for anti-aliasing and soft shadows
  srand(time(NULL));

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

