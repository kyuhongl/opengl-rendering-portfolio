/*
  CSCI 420 Computer Graphics, Computer Science, USC
  Assignment 2: Roller Coaster Simulation.
  C/C++ starter code
  
  Student username: kyuhongl
*/

//*****************************
// include files
//*****************************
// core opengl includes
#include "openGLHeader.h"
#include "glutHeader.h"

#include "openGLMatrix.h"
#include "imageIO.h"
#include "pipelineProgram.h"
#include "vbo.h"
#include "vao.h"
#include "ebo.h"

#include <sys/stat.h>
#include <string>
#include <iostream>
#include <cstring>
#include <unistd.h> 

#if defined(WIN32) || defined(_WIN32)
#ifdef _DEBUG
#pragma comment(lib, "glew32d.lib")
#else
#pragma comment(lib, "glew32.lib")
#endif
#endif

#if defined(WIN32) || defined(_WIN32)
char shaderBasePath[1024] = SHADER_BASE_PATH;
#else
char shaderBasePath[1024] = "../openGLHelper";
#endif

using namespace std;

//*****************************
// forward declarations
//*****************************
void cleanup();
void computeLowestPoint();


struct Point
{
  float x, y, z;
};

struct Spline
{
  int numControlPoints;
  Point *points;
} spline;


// camera control modes
typedef enum
{
  ROTATE,
  TRANSLATE,
  SCALE
} CONTROL_STATE;

// bird's eye view control modes
enum BIRD_EYE_CONTROL
{
  BIRD_MOVE,
  BIRD_ZOOM
};

//*****************************
// global variables - input and camera control
//*****************************
// mouse tracking
int mousePos[2] = {0, 0};  // x,y coordinates of current mouse position
int leftMouseButton = 0;   // 1 if pressed, 0 if not
int middleMouseButton = 0; // 1 if pressed, 0 if not
int rightMouseButton = 0;  // 1 if pressed, 0 if not
CONTROL_STATE controlState = ROTATE;  // current control state for mouse interaction
BIRD_EYE_CONTROL birdEyeControl = BIRD_MOVE;  // bird's eye view control mode

// window parameters
int windowWidth = 1280;
int windowHeight = 720;
char windowTitle[512] = "CSCI 420 Homework 2";

// terrain transformation parameters
float terrainRotate[3] = {0.0f, 0.0f, 0.0f};  // rotation around x, y, z axes (in degrees)
float terrainTranslate[3] = {0.0f, 0.0f, 0.0f};  // translation along x, y, z axes
float terrainScale[3] = {1.0f, 1.0f, 1.0f};  // scaling factors for x, y, z axes

// camera parameters
float cameraP = 0.0f;  // parametric position on the spline (0 to 1)
float cameraSpeed = 0.01f;  // speed of camera movement along the spline
int currentSpline = 0;  // current spline segment being traversed
bool birdEyeView = false;  // toggle for bird's eye view
float birdEyeHeight = 15.0f;  // height of bird's eye view
const float BIRD_EYE_FOV = 90.0f;  // field of view for bird's eye view
const float NORMAL_FOV = 60.0f;  // field of view for normal camera

//*****************************
// global variables - geometry and rendering
//*****************************
// geometry counts
int numVertices;  // total number of vertices
int numIndices;  // total number of indices

// rail and tie parameters
float railWidth = 0.2f;  // width of each rail
float railHeight = 0.2f;  // height of each rail
float railSeparation = 1.0f;  // distance between the two rails
float tieWidth = 0.2f;  // width of each tie
float tieHeight = 0.05f;  // height of each tie
float tieSpacing = 0.2f;  // distance between ties

// geometry data
vector<float> pointPos, pointCol;  // positions and colors for rails
vector<unsigned int> railIndices, railColors;  // indices for rendering rails

// tie geometry
vector<float> tiePositions;  // positions for ties
vector<float> tieColors;  // colors for ties
vector<unsigned int> tieIndices;  // indices for rendering ties
int numTieVertices = 0;  // number of vertices for ties

// ground plane
vector<float> positions = {
    -100.0f, -5.0f, -100.0f,  // bottom-left
    100.0f, -5.0f, -100.0f,   // bottom-right
    100.0f, -5.0f, 100.0f,    // top-right
    -100.0f, -5.0f, 100.0f    // top-left
};

vector<float> texCoords = {
    0.0f, 0.0f,  // bottom-left
    1.0f, 0.0f,  // bottom-right
    1.0f, 1.0f,  // top-right
    0.0f, 1.0f   // top-left
};

vector<unsigned int> groundIndices = {
    0, 1, 2,  // first triangle
    2, 3, 0   // second triangle
};

// skybox vertices (cube centered at origin)
std::vector<float> skyboxVertices = {
  // Positions         // Texture Coords
  // Back face (Negative Z)
  -1.0f,  1.0f, -1.0f, 0.0f, 1.0f, // 0 - Top-left
  -1.0f, -1.0f, -1.0f, 0.0f, 0.0f, // 1 - Bottom-left
   1.0f, -1.0f, -1.0f, 1.0f, 0.0f, // 2 - Bottom-right
   1.0f,  1.0f, -1.0f, 1.0f, 1.0f, // 3 - Top-right

  // Front face (Positive Z)
  -1.0f,  1.0f, 1.0f, 1.0f, 1.0f, // 4 - Top-left
  -1.0f, -1.0f, 1.0f, 1.0f, 0.0f, // 5 - Bottom-left
   1.0f, -1.0f, 1.0f, 0.0f, 0.0f, // 6 - Bottom-right
   1.0f,  1.0f, 1.0f, 0.0f, 1.0f, // 7 - Top-right

  // Left face (Negative X)
  -1.0f,  1.0f,  1.0f, 0.0f, 1.0f, // 8 - Top-left
  -1.0f, -1.0f,  1.0f, 0.0f, 0.0f, // 9 - Bottom-left
  -1.0f, -1.0f, -1.0f, 1.0f, 0.0f, // 10 - Bottom-right
  -1.0f,  1.0f, -1.0f, 1.0f, 1.0f, // 11 - Top-right

  // Right face (Positive X)
  1.0f,  1.0f, -1.0f, 0.0f, 1.0f, // 12 - Top-left
  1.0f, -1.0f, -1.0f, 0.0f, 0.0f, // 13 - Bottom-left
  1.0f, -1.0f,  1.0f, 1.0f, 0.0f, // 14 - Bottom-right
  1.0f,  1.0f,  1.0f, 1.0f, 1.0f, // 15 - Top-right

  // Bottom face (Negative Y)
  -1.0f, -1.0f, -1.0f, 0.0f, 1.0f, // 16 - Top-left
   1.0f, -1.0f, -1.0f, 1.0f, 1.0f, // 17 - Top-right
   1.0f, -1.0f,  1.0f, 1.0f, 0.0f, // 18 - Bottom-right
  -1.0f, -1.0f,  1.0f, 0.0f, 0.0f, // 19 - Bottom-left

  // Top face (Positive Y)
  -1.0f,  1.0f, -1.0f, 0.0f, 0.0f, // 20 - Bottom-left
   1.0f,  1.0f, -1.0f, 1.0f, 0.0f, // 21 - Bottom-right
   1.0f,  1.0f,  1.0f, 1.0f, 1.0f, // 22 - Top-right
  -1.0f,  1.0f,  1.0f, 0.0f, 1.0f  // 23 - Top-left
};

// skybox texture paths - corrected orientation
vector<string> skyboxFaces = {
    "textures/skybox/right.jpg",  // positive x (right)
    "textures/skybox/left.jpg",   // negative x (left)
    "textures/skybox/bottom.jpg",    // positive y (top)
    "textures/skybox/top.jpg", // negative y (bottom)
    "textures/skybox/front.jpg",  // positive z (front)
    "textures/skybox/back.jpg"    // negative z (back)
};

// recording parameters
int frameCount = 0;
const int MAX_FRAMES = 450;      // Increased to 450 frames
const int FRAMES_PER_SECOND = 15; // Record at 15 frames per second
bool isRecording = false;
const char *outputDirectory = "frames/";

// rendering flags
bool showSkybox = true;
bool showGround = false;

//*****************************
// OpenGL objects
//*****************************
// textures
GLuint textureID;
GLuint groundTextureHandle;
GLuint skyboxTexture = 0;

// helper objects
OpenGLMatrix matrix;
PipelineProgram *pipelineProgram = nullptr;
PipelineProgram *texturePipelineProgram = nullptr;
PipelineProgram *skyboxPipeline = nullptr;

// vertex buffer objects
VAO *splineVAO = nullptr;
VBO *splineVBO = nullptr;
VBO *colorVBO = nullptr;
EBO *railEBO = nullptr;

VAO *groundVAO = nullptr;
VBO *groundVBOPosition = nullptr;
VBO *groundVBOTexCoord = nullptr;
EBO *groundEBO = nullptr;

VAO *tieVAO = nullptr;
VBO *tieVBOPosition = nullptr;
VBO *tieVBOColor = nullptr;
EBO *tieEBO = nullptr;

// Define skybox indices
std::vector<unsigned int> skyboxIndices = {
  // Back face
  0, 1, 2, 0, 2, 3,
  // Front face
  4, 5, 6, 4, 6, 7,
  // Left face
  8, 9, 10, 8, 10, 11,
  // Right face
  12, 13, 14, 12, 14, 15,
  // Bottom face
  16, 17, 18, 16, 18, 19,
  // Top face
  20, 21, 22, 20, 22, 23
};

// skybox VAO and VBO for cubemap
// add EBO for skybox
VAO* skyboxVAO = nullptr;
VBO* skyboxVBO = nullptr;
GLuint skyboxEBO = 0;

// necessary variables
Spline *splines;
int numSplines;

//*****************************
// math helper functions
//*****************************

// multiply matrices: c = a * b, where a is m×p and b is p×n
// matrices are stored in column-major order
void MultiplyMatrices(int m, int p, int n, const float *A, const float *B, float *C)
{
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      float entry = 0.0;
      for (int k = 0; k < p; k++)
        entry += A[k * m + i] * B[j * p + k];
      C[m * j + i] = entry;
    }
  }
}

// normalize a 3d vector to unit length
void normalize(float *vec)
{
  float length = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
  if (length != 0.0f)
  {
    vec[0] /= length;
    vec[1] /= length;
    vec[2] /= length;
  }
  else
  {
    // default to z-axis if vector has zero length
    vec[0] = 0.0f;
    vec[1] = 0.0f;
    vec[2] = 1.0f;
  }
}

// compute cross product: result = a × b
void crossProduct(const float *a, const float *b, float *result)
{
  result[0] = a[1] * b[2] - a[2] * b[1];
  result[1] = a[2] * b[0] - a[0] * b[2];
  result[2] = a[0] * b[1] - a[1] * b[0];
}

// compute tangent vector at parameter u on the spline
void computeTangent(const float u, const float *M, const float *C, float *T)
{
  // calculate derivative of spline at parameter u (3u², 2u, 1, 0)
  float Uprime[4] = {3 * u * u, 2 * u, 1.0f, 0.0f};

  // multiply by basis matrix to get tangent
  float UM[4];
  MultiplyMatrices(1, 4, 4, Uprime, M, UM);
  MultiplyMatrices(1, 4, 3, UM, C, T);

  // normalize to unit length
  normalize(T);
}

// compute binormal vector from tangent and previous normal
void computeBinormal(const float *tangent, const float *prevNormal, float *binormal)
{
  crossProduct(tangent, prevNormal, binormal);
  normalize(binormal);
}

// compute normal vector to complete the orthonormal frame
void computeNormal(const float *tangent, float *normal, float *binormal)
{
  crossProduct(binormal, tangent, normal);
  normalize(normal);
}

// function to compute normal matrix from modelview matrix
void computeNormalMatrix(const float* modelViewMatrix, float* normalMatrix)
{
  // extract the upper-left 3x3 matrix
  float m[9] = {
    modelViewMatrix[0], modelViewMatrix[1], modelViewMatrix[2],
    modelViewMatrix[4], modelViewMatrix[5], modelViewMatrix[6],
    modelViewMatrix[8], modelViewMatrix[9], modelViewMatrix[10]
  };

  // for orthogonal matrices, the inverse is equal to the transpose
  // here we assume no non-uniform scaling
  memcpy(normalMatrix, m, 9 * sizeof(float));
}

//*****************************
// geometry generation functions
//*****************************

// generate the four corners of a rail segment at a given position
void generateRailCorners(const float *center, const float *normal, const float *binormal,
                         const float *railCenter, float corners[4][3])
{
  // bottom-left corner
  corners[0][0] = railCenter[0] - binormal[0] * railWidth / 2 - normal[0] * railHeight / 2;
  corners[0][1] = railCenter[1] - binormal[1] * railWidth / 2 - normal[1] * railHeight / 2;
  corners[0][2] = railCenter[2] - binormal[2] * railWidth / 2 - normal[2] * railHeight / 2;

  // bottom-right corner
  corners[1][0] = railCenter[0] + binormal[0] * railWidth / 2 - normal[0] * railHeight / 2;
  corners[1][1] = railCenter[1] + binormal[1] * railWidth / 2 - normal[1] * railHeight / 2;
  corners[1][2] = railCenter[2] + binormal[2] * railWidth / 2 - normal[2] * railHeight / 2;

  // top-right corner
  corners[2][0] = railCenter[0] + binormal[0] * railWidth / 2 + normal[0] * railHeight / 2;
  corners[2][1] = railCenter[1] + binormal[1] * railWidth / 2 + normal[1] * railHeight / 2;
  corners[2][2] = railCenter[2] + binormal[2] * railWidth / 2 + normal[2] * railHeight / 2;

  // top-left corner
  corners[3][0] = railCenter[0] - binormal[0] * railWidth / 2 + normal[0] * railHeight / 2;
  corners[3][1] = railCenter[1] - binormal[1] * railWidth / 2 + normal[1] * railHeight / 2;
  corners[3][2] = railCenter[2] - binormal[2] * railWidth / 2 + normal[2] * railHeight / 2;
}

// add vertices for a rail segment to the geometry buffers
void addRailVertices(const float corners[4][3], float r, float g, float b, const float* normal)
{
  for (int j = 0; j < 4; j++)
  {
    // add position coordinates
    pointPos.push_back(corners[j][0]);
    pointPos.push_back(corners[j][1]);
    pointPos.push_back(corners[j][2]);

    // store normal in RGB channels (mapped from [-1,1] to [0,1])
    pointCol.push_back(normal[0] * 0.5f + 0.5f);
    pointCol.push_back(normal[1] * 0.5f + 0.5f);
    pointCol.push_back(normal[2] * 0.5f + 0.5f);
    pointCol.push_back(1.0f);
  }
}

// generate a cylindrical tie between rail segments
void generateCylindricalTie(const float *leftRailCenter, const float *rightRailCenter,
                          const float *normal, const float *binormal)
{
  const int numSides = 12;  // number of sides for the cylinder
  const float tieRadius = 0.055f;  // cylinder radius
  
  // store starting vertex index for this tie
  int baseVertex = tiePositions.size() / 3;
  
  // calculate center point between rails
  float tieCenter[3];
  float tieVector[3]; // vector from left to right rail
  float tieLength = 0.0f;
  
  // compute vector between rails and its length
  for (int i = 0; i < 3; i++) {
    tieVector[i] = rightRailCenter[i] - leftRailCenter[i];
    tieLength += tieVector[i] * tieVector[i];
    tieCenter[i] = (leftRailCenter[i] + rightRailCenter[i]) * 0.5f;
  }
  tieLength = sqrt(tieLength);
  
  // normalize the tie vector
  for (int i = 0; i < 3; i++) {
    tieVector[i] /= tieLength;
  }
  
  // move tie center slightly down
  for (int i = 0; i < 3; i++) {
    tieCenter[i] -= normal[i] * (tieRadius * 0.5f);
  }
  
  // create a new coordinate system for the cylinder
  float xAxis[3], yAxis[3], zAxis[3];
  
  // z axis is along the tie
  zAxis[0] = tieVector[0];
  zAxis[1] = tieVector[1];
  zAxis[2] = tieVector[2];
  
  // use normal as a starting point for y axis
  yAxis[0] = normal[0];
  yAxis[1] = normal[1];
  yAxis[2] = normal[2];
  
  // x axis is orthogonal to both
  crossProduct(yAxis, zAxis, xAxis);
  normalize(xAxis);
  
  // recalculate y to ensure orthogonality
  crossProduct(zAxis, xAxis, yAxis);
  normalize(yAxis);
  
  // create vertices for each end of the cylinder
  // make cylinder slightly shorter than rail separation
  float leftEnd[3], rightEnd[3];
  for (int i = 0; i < 3; i++) {
    leftEnd[i] = tieCenter[i] - zAxis[i] * (tieLength * 0.52f);
    rightEnd[i] = tieCenter[i] + zAxis[i] * (tieLength * 0.52f);
  }
  
  // create vertices for the cylinder sides
  for (int side = 0; side < numSides; side++) {
    float angle = 2.0f * M_PI * side / numSides;
    float cosA = cos(angle);
    float sinA = sin(angle);
    
    // calculate offset from center line
    float offsetX[3], offsetY[3];
    for (int i = 0; i < 3; i++) {
      offsetX[i] = xAxis[i] * cosA * tieRadius;
      offsetY[i] = yAxis[i] * sinA * tieRadius;
    }
    
    // left end vertices
    for (int i = 0; i < 3; i++) {
      tiePositions.push_back(leftEnd[i] + offsetX[i] + offsetY[i]);
    }
    
    // dark wooden color
    tieColors.push_back(0.35f);
    tieColors.push_back(0.18f);
    tieColors.push_back(0.05f);
    tieColors.push_back(1.0f);
    
    // right end vertices
    for (int i = 0; i < 3; i++) {
      tiePositions.push_back(rightEnd[i] + offsetX[i] + offsetY[i]);
    }
    
    // same color for consistency
    tieColors.push_back(0.35f);
    tieColors.push_back(0.18f);
    tieColors.push_back(0.05f);
    tieColors.push_back(1.0f);
  }
  
  // generate indices for the cylinder sides
  for (int i = 0; i < numSides; i++) {
    int i2 = (i * 2);
    int next_i2 = ((i + 1) % numSides) * 2;
    
    // first triangle
    tieIndices.push_back(baseVertex + i2);
    tieIndices.push_back(baseVertex + i2 + 1);
    tieIndices.push_back(baseVertex + next_i2);
    
    // second triangle
    tieIndices.push_back(baseVertex + i2 + 1);
    tieIndices.push_back(baseVertex + next_i2 + 1);
    tieIndices.push_back(baseVertex + next_i2);
  }
   
  // add end caps
  // left cap center
  int leftCenterIdx = tiePositions.size() / 3;
  for (int i = 0; i < 3; i++) {
    tiePositions.push_back(leftEnd[i]);
  }
  tieColors.push_back(0.3f); // darker brown for caps
  tieColors.push_back(0.15f);
  tieColors.push_back(0.0f);
  tieColors.push_back(1.0f);
   
  // right cap center
  int rightCenterIdx = leftCenterIdx + 1;
  for (int i = 0; i < 3; i++) {
    tiePositions.push_back(rightEnd[i]);
  }
  tieColors.push_back(0.3f);
  tieColors.push_back(0.15f);
  tieColors.push_back(0.0f);
  tieColors.push_back(1.0f);
   
  // add triangles for left cap
  for (int i = 0; i < numSides; i++) {
    int next = (i + 1) % numSides;
    tieIndices.push_back(leftCenterIdx);
    tieIndices.push_back(baseVertex + next * 2);
    tieIndices.push_back(baseVertex + i * 2);
  }
   
  // add triangles for right cap
  for (int i = 0; i < numSides; i++) {
    int next = (i + 1) % numSides;
    tieIndices.push_back(rightCenterIdx);
    tieIndices.push_back(baseVertex + i * 2 + 1);
    tieIndices.push_back(baseVertex + next * 2 + 1);
  }
}

// add indices for connecting rail segments
void addRailIndices(int baseVertex, int prevBaseVertex)
{
  // bottom face
  railIndices.push_back(prevBaseVertex);
  railIndices.push_back(baseVertex);
  railIndices.push_back(baseVertex + 1);

  railIndices.push_back(prevBaseVertex);
  railIndices.push_back(baseVertex + 1);
  railIndices.push_back(prevBaseVertex + 1);

  // side faces (4 sides of the rail)
  for (int side = 0; side < 4; side++)
  {
    int current = baseVertex + side;
    int next = baseVertex + ((side + 1) % 4);
    int prevCurrent = prevBaseVertex + side;
    int prevNext = prevBaseVertex + ((side + 1) % 4);

    // first triangle
    railIndices.push_back(prevCurrent);
    railIndices.push_back(current);
    railIndices.push_back(next);

    // second triangle
    railIndices.push_back(prevCurrent);
    railIndices.push_back(next);
    railIndices.push_back(prevNext);
  }
}

// compute the catmull-rom spline and generate the full roller coaster geometry
void computeCatmullSpline()
{
  // clear previous geometry data
  pointPos.clear();
  pointCol.clear();
  railIndices.clear();
  tiePositions.clear();
  tieColors.clear();
  tieIndices.clear();

  // catmull-rom basis matrix
  float basisM[16] = {
      -0.5f, 1.0f, -0.5f, 0.0f,
      1.5f, -2.5f, 0.0f, 1.0f,
      -1.5f, 2.0f, 0.5f, 0.0f,
      0.5f, -0.5f, 0.0f, 0.0f};

  // initialize parameters for geometry generation
  bool isFirstPointOverall = true;
  float prevNormal[3] = {0.0f, 1.0f, 0.0f};  // start with up vector
  float tieDistance = 0.0f;  // accumulated distance for placing ties
  numVertices = 0;

  // for each spline segment defined by 4 control points
  for (int i = 0; i < spline.numControlPoints - 3; i++)
  {
    // get control points for current segment
    Point p1 = spline.points[i];
    Point p2 = spline.points[i + 1];
    Point p3 = spline.points[i + 2];
    Point p4 = spline.points[i + 3];

    // set up control matrix with positions
    float controlM[12] = {
        p1.x, p2.x, p3.x, p4.x,
        p1.y, p2.y, p3.y, p4.y,
        p1.z, p2.z, p3.z, p4.z};

    float u = 0.0f;  // parameter along spline (0 to 1)
    float step = 0.01f;  // step size for sampling the spline
    bool isFirstPointInSegment = true;
    unsigned int lastSegmentVertexCount = numVertices;

    // copy last vertices from previous segment for continuity if not the first segment
    if (i > 0)
    {
      for (int j = 0; j < 8; j++)
      {
        pointPos.push_back(pointPos[lastSegmentVertexCount - 24 + j * 3]);
        pointPos.push_back(pointPos[lastSegmentVertexCount - 24 + j * 3 + 1]);
        pointPos.push_back(pointPos[lastSegmentVertexCount - 24 + j * 3 + 2]);

        pointCol.push_back(pointCol[lastSegmentVertexCount - 32 + j * 4]);
        pointCol.push_back(pointCol[lastSegmentVertexCount - 32 + j * 4 + 1]);
        pointCol.push_back(pointCol[lastSegmentVertexCount - 32 + j * 4 + 2]);
        pointCol.push_back(pointCol[lastSegmentVertexCount - 32 + j * 4 + 3]);
      }
      numVertices += 8;  // 8 vertices added (4 for each rail)
    }

    // generate points along the spline at regular intervals
    while (u <= 1.01f)  // slightly over 1 to ensure we include the endpoint
    {
      // calculate position on spline using cubic polynomial
      float U[4] = {u * u * u, u * u, u, 1.0f};  // u³, u², u, 1
      float UM[4];  // will hold U * basisM
      float center[3];  // will hold the final position
      MultiplyMatrices(1, 4, 4, U, basisM, UM);
      MultiplyMatrices(1, 4, 3, UM, controlM, center);

      // calculate frame vectors (tangent, normal, binormal)
      float tangent[3];
      computeTangent(u, basisM, controlM, tangent);
      normalize(tangent);

      float binormal[3];
      float normal[3];
      computeBinormal(tangent, prevNormal, binormal);
      computeNormal(tangent, normal, binormal);

      // ensure smooth normal transition (flip if necessary)
      if (!isFirstPointOverall)
      {
        float dotProduct = prevNormal[0] * normal[0] +
                           prevNormal[1] * normal[1] +
                           prevNormal[2] * normal[2];
        if (dotProduct < 0)
        {
          // flip normal if it flipped direction
          normal[0] = -normal[0];
          normal[1] = -normal[1];
          normal[2] = -normal[2];
        }
      }

      // store normal for next iteration
      prevNormal[0] = normal[0];
      prevNormal[1] = normal[1];
      prevNormal[2] = normal[2];

      // calculate rail centers (offset from center by half the rail separation)
      float leftRailCenter[3], rightRailCenter[3];
      for (int j = 0; j < 3; j++)
      {
        leftRailCenter[j] = center[j] - binormal[j] * (railSeparation / 2);
        rightRailCenter[j] = center[j] + binormal[j] * (railSeparation / 2);
      }

      // generate rail corners for left and right rails
      float leftCorners[4][3], rightCorners[4][3];
      generateRailCorners(center, normal, binormal, leftRailCenter, leftCorners);
      generateRailCorners(center, normal, binormal, rightRailCenter, rightCorners);

      // add vertices for both rails
      addRailVertices(leftCorners, 0.3f, 0.45f, 0.2f, normal);  // dark greenish-brown color
      addRailVertices(rightCorners, 0.3f, 0.45f, 0.2f, normal);

      // generate indices for rails (connect to previous point)
      if (!isFirstPointOverall)
      {
        if (isFirstPointInSegment)
        {
          // connect to last point of previous segment
          addRailIndices(numVertices, lastSegmentVertexCount - 8);
          addRailIndices(numVertices + 4, lastSegmentVertexCount - 4);
          isFirstPointInSegment = false;
        }
        else
        {
          // connect to previous point in this segment
          addRailIndices(numVertices, numVertices - 8);
          addRailIndices(numVertices + 4, numVertices - 4);
        }
      }

      // add ties at regular intervals
      tieDistance += step;
      if (tieDistance >= tieSpacing)
      {
        tieDistance = 0.0f;

        // only place ties away from segment boundaries
        if (u > 0.05f && u < 0.95f)
        {
          generateCylindricalTie(leftRailCenter, rightRailCenter, normal, binormal);
        }
      }

      isFirstPointOverall = false;
      numVertices += 8;  // 8 vertices added (4 for each rail)
      u += step;
    }
  }
  
  // debug output
  std::cout << "generated " << pointPos.size() / 3 << " vertices for rails" << std::endl;
  std::cout << "generated " << railIndices.size() << " indices for rails" << std::endl;
  std::cout << "generated " << tiePositions.size() / 3 << " vertices for ties" << std::endl;
  std::cout << "generated " << tieIndices.size() << " indices for ties" << std::endl;
}

//*****************************
// camera and view functions
//*****************************

// set up camera for bird's eye view
void setBirdViewCamera()
{
  // set to 45 degree elevation angle
  float elevationAngle = 45.0f * (M_PI / 180.0f);
  float azimuthAngle = -90.0f * (M_PI / 180.0f);

  // calculate camera position in spherical coordinates
  float cameraX = birdEyeHeight * cos(elevationAngle) * sin(azimuthAngle);
  float cameraY = birdEyeHeight * sin(elevationAngle);
  float cameraZ = birdEyeHeight * cos(elevationAngle) * cos(azimuthAngle);

  // set up modelview matrix
  matrix.SetMatrixMode(OpenGLMatrix::ModelView);
  matrix.LoadIdentity();
  matrix.LookAt(
      cameraX, cameraY, cameraZ, // camera position
      0.0f, 0.0f, 0.0f,          // look-at point (origin)
      0.0f, 1.0f, 0.0f           // up vector
  );
}

// helper function to set up camera at a position on the track
void cameraHelper(const float *position, const float *tangent, const float *normal, const float *binormal)
{
  // camera positioning parameters
  float heightOffset = 0.5f; // height above track
  
  // calculate camera position (slightly above track)
  float cameraPosition[3] = {
      position[0] + normal[0] * heightOffset,
      position[1] + normal[1] * heightOffset,
      position[2] + normal[2] * heightOffset};

  // calculate look-at point (ahead of camera position)
  float lookAt[3] = {
      cameraPosition[0] + tangent[0],
      cameraPosition[1] + tangent[1],
      cameraPosition[2] + tangent[2]};

  // set the camera matrix
  matrix.SetMatrixMode(OpenGLMatrix::ModelView);
  matrix.LoadIdentity();
  matrix.LookAt(
      cameraPosition[0], cameraPosition[1], cameraPosition[2], // camera position
      lookAt[0], lookAt[1], lookAt[2],                         // look-at point
      normal[0], normal[1], normal[2]                          // up vector
  );
}

// compute camera position and orientation at parameter u on spline segment
void computeCamera(float u, int currentSpline)
{
  // if in bird's eye view, don't calculate roller coaster camera
  if (birdEyeView)
  {
    return;
  }

  // get control points for current segment
  Point p1 = spline.points[currentSpline];
  Point p2 = spline.points[currentSpline + 1];
  Point p3 = spline.points[currentSpline + 2];
  Point p4 = spline.points[currentSpline + 3];

  // set up control matrix
  float controlM[12] = {
      p1.x, p2.x, p3.x, p4.x,
      p1.y, p2.y, p3.y, p4.y,
      p1.z, p2.z, p3.z, p4.z};

  // catmull-rom basis matrix
  float basisM[16] = {
      -0.5f, 1.0f, -0.5f, 0.0f,
      1.5f, -2.5f, 0.0f, 1.0f,
      -1.5f, 2.0f, 0.5f, 0.0f,
      0.5f, -0.5f, 0.0f, 0.0f};

  // calculate position on spline
  float U[4] = {u * u * u, u * u, u, 1.0f};
  float UM[4];
  float position[3];
  MultiplyMatrices(1, 4, 4, U, basisM, UM);
  MultiplyMatrices(1, 4, 3, UM, controlM, position);

  // calculate tangent, normal, and binormal
  float tangent[3];
  computeTangent(u, basisM, controlM, tangent);
  normalize(tangent);

  // use persistent normal to maintain continuity
  static float prevNormal[3] = {0.0f, 1.0f, 0.0f};
  float binormal[3];
  float normal[3];

  computeBinormal(tangent, prevNormal, binormal);
  computeNormal(tangent, normal, binormal);

  // store normal for next iteration
  prevNormal[0] = normal[0];
  prevNormal[1] = normal[1];
  prevNormal[2] = normal[2];

  // update camera
  cameraHelper(position, tangent, normal, binormal);
}

//*****************************
// texture and material functions
//*****************************

// For debugging the skybox coordinate system
void debugSkyboxCoordinates() {
  std::cout << "Skybox direction vectors:" << std::endl;
  std::cout << "  Right (+X): " << "vec3(1.0, 0.0, 0.0)" << std::endl;
  std::cout << "  Left  (-X): " << "vec3(-1.0, 0.0, 0.0)" << std::endl;
  std::cout << "  Top   (+Y): " << "vec3(0.0, 1.0, 0.0)" << std::endl;
  std::cout << "  Bottom(-Y): " << "vec3(0.0, -1.0, 0.0)" << std::endl;
  std::cout << "  Front (+Z): " << "vec3(0.0, 0.0, 1.0)" << std::endl;
  std::cout << "  Back  (-Z): " << "vec3(0.0, 0.0, -1.0)" << std::endl;
}

// load a cubemap texture from a set of image files
GLuint loadCubemap(std::vector<std::string> faces)
{
    GLuint textureID;
    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_CUBE_MAP, textureID);

    for (unsigned int i = 0; i < faces.size(); i++)
    {
        ImageIO img;
        ImageIO::fileFormatType imgFormat;
        ImageIO::errorType err = img.load(faces[i].c_str(), &imgFormat);

        if (err == ImageIO::OK)
        {
            glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 
                         0, GL_RGB, img.getWidth(), img.getHeight(), 0, 
                         GL_RGB, GL_UNSIGNED_BYTE, img.getPixels());
            printf("loaded skybox texture: %s\n", faces[i].c_str());
        }
        else
        {
            printf("failed to load skybox texture: %s\n", faces[i].c_str());
        }
    }

    // set texture parameters
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    return textureID;
}

// create a directory for storing screenshot frames
void createDirectory(const char *path)
{
  mkdir(path, 0777);
}

// load a texture from a file and set up texture parameters
void loadTexture(const char *filename)
{
  ImageIO *textureImage = new ImageIO();
  if (textureImage->loadJPEG(filename) != ImageIO::OK)
  {
    cout << "error loading texture image." << endl;
    return;
  }

  glGenTextures(1, &textureID);
  glBindTexture(GL_TEXTURE_2D, textureID);

  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB,
               textureImage->getWidth(), textureImage->getHeight(),
               0, GL_RGB, GL_UNSIGNED_BYTE, textureImage->getPixels());

  // set texture parameters
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

  delete textureImage;
}

int initTexture(const char *imageFilename, GLuint textureHandle)
{
  ImageIO img;
  ImageIO::fileFormatType imgFormat;
  ImageIO::errorType err = img.load(imageFilename, &imgFormat);

  if (err != ImageIO::OK)
  {
    printf("loading texture from %s failed.\n", imageFilename);
    return -1;
  }

  if (img.getWidth() * img.getBytesPerPixel() % 4)
  {
    printf("error (%s): the width*numChannels in the loaded image must be a multiple of 4.\n", imageFilename);
    return -1;
  }

  int width = img.getWidth();
  int height = img.getHeight();
  unsigned char *pixelsRGBA = new unsigned char[4 * width * height]; // 4 bytes per pixel (RGBA)

  memset(pixelsRGBA, 0, 4 * width * height); // set all bytes to 0
  for (int h = 0; h < height; h++)
    for (int w = 0; w < width; w++)
    {
      pixelsRGBA[4 * (h * width + w) + 0] = 0;   // red
      pixelsRGBA[4 * (h * width + w) + 1] = 0;   // green
      pixelsRGBA[4 * (h * width + w) + 2] = 0;   // blue
      pixelsRGBA[4 * (h * width + w) + 3] = 255; // alpha channel; fully opaque

      int numChannels = img.getBytesPerPixel();
      for (int c = 0; c < numChannels; c++) // only set as many channels as are available
        pixelsRGBA[4 * (h * width + w) + c] = img.getPixel(w, h, c);
    }

  glBindTexture(GL_TEXTURE_2D, textureHandle);

  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixelsRGBA);


  glGenerateMipmap(GL_TEXTURE_2D);


  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  

  GLfloat fLargest;
  glGetFloatv(GL_MAX_TEXTURE_MAX_ANISOTROPY_EXT, &fLargest);
  printf("max available anisotropic samples: %f\n", fLargest);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, 0.5f * fLargest);


  GLenum errCode = glGetError();
  if (errCode != 0)
  {
    printf("texture initialization error. error code: %d.\n", errCode);
    return -1;
  }


  delete[] pixelsRGBA;

  return 0;
}

// write a screenshot to the specified filename
void saveScreenshot(const char *filename)
{
  unsigned char *screenshotData = new unsigned char[windowWidth * windowHeight * 3];
  glReadPixels(0, 0, windowWidth, windowHeight, GL_RGB, GL_UNSIGNED_BYTE, screenshotData);

  ImageIO screenshotImg(windowWidth, windowHeight, 3, screenshotData);

  if (screenshotImg.save(filename, ImageIO::FORMAT_JPEG) == ImageIO::OK)
    cout << "file " << filename << " saved successfully." << endl;
  else
    cout << "failed to save file " << filename << '.' << endl;

  delete[] screenshotData;
}

//*****************************
// file loading functions
//*****************************


void loadSpline(char *argv)
{
  FILE *fileSpline = fopen(argv, "r");
  if (fileSpline == NULL)
  {
    printf("cannot open file %s.\n", argv);
    exit(1);
  }


  fscanf(fileSpline, "%d\n", &spline.numControlPoints);
  printf("detected %d control points.\n", spline.numControlPoints);


  spline.points = (Point *)malloc(spline.numControlPoints * sizeof(Point));
  

  for (int i = 0; i < spline.numControlPoints; i++)
  {
    if (fscanf(fileSpline, "%f %f %f",
               &spline.points[i].x,
               &spline.points[i].y,
               &spline.points[i].z) != 3)
    {
      printf("error: incorrect number of control points in file %s.\n", argv);
      exit(1);
    }
  }
}

// compute the lowest point on the spline (for ground positioning)
// doesnt really worok
void computeLowestPoint()
{
  float lowestY = std::numeric_limits<float>::max();
  for (int i = 0; i < spline.numControlPoints; i++)
  {
    if (spline.points[i].y < lowestY)
    {
      lowestY = spline.points[i].y;
    }
  }
  std::cout << "lowest spline point y: " << lowestY << std::endl;
}

//*****************************
// event callbacks
//*****************************

// glut idle callback for animation
void idleFunc()
{

  if (birdEyeView)
  {
    // in bird's eye view, we don't move the camera
    setBirdViewCamera();
  }
  
  // move camera along the spline
  cameraP += cameraSpeed;
  if (cameraP > 1.0f)
  {
    cameraP -= 1.0f;
    currentSpline++;
    if (currentSpline >= spline.numControlPoints - 3)
    {
      currentSpline = 0;
    }
  }

  // compute new camera position
  computeCamera(cameraP, currentSpline);

  // Handle recording if enabled
  if (isRecording)
  {
    // Get elapsed time and calculate frame rate control
    static int lastFrameTime = glutGet(GLUT_ELAPSED_TIME);
    int currentTime = glutGet(GLUT_ELAPSED_TIME);
    int elapsedMs = currentTime - lastFrameTime;
    
    // Calculate milliseconds per frame for target frame rate
    int msPerFrame = 1000 / FRAMES_PER_SECOND;
    
    // Only capture frame when enough time has passed
    if (elapsedMs >= msPerFrame)
    {
      // Update the last frame time
      lastFrameTime = currentTime;
      
      // Create the output directory if it doesn't exist
      static bool dirCreated = false;
      if (!dirCreated)
      {
        createDirectory(outputDirectory);
        dirCreated = true;
      }
      
      // Generate filename
      char filename[1024];
      sprintf(filename, "%sframe%04d.jpg", outputDirectory, frameCount);
      
      // Save the screenshot
      saveScreenshot(filename);
      std::cout << "Saved frame " << frameCount + 1 << " of " << MAX_FRAMES << std::endl;
      
      // Increment frame counter and stop recording if maximum is reached
      frameCount++;
      if (frameCount >= MAX_FRAMES)
      {
        std::cout << "Recording complete: " << MAX_FRAMES << " frames saved to " 
                  << outputDirectory << std::endl;
        isRecording = false;
        frameCount = 0;
      }
    }
  }
  else
  {
    // Reset the frame counter if recording is off
    frameCount = 0;
  }

  // request display update
  glutPostRedisplay();
}

// glut reshape callback - called when window is resized
void reshapeFunc(int w, int h)
{
  glViewport(0, 0, w, h);

  // when window has been resized, we need to re-set the projection matrix
  matrix.SetMatrixMode(OpenGLMatrix::Projection);
  matrix.LoadIdentity();
  
  // set up perspective projection
  const float zNear = 0.1f;
  const float zFar = 10000.0f;
  const float fieldOfView = birdEyeView ? BIRD_EYE_FOV : NORMAL_FOV;
  matrix.Perspective(fieldOfView, 1.0f * w / h, zNear, zFar);
}

// glut mouse motion callback (with button pressed)
void mouseMotionDragFunc(int x, int y)
{
  // mouse has moved, and one of the mouse buttons is pressed
  cout << "dragging. x: " << x << " y:" << y << endl;

  // the change in mouse position since the last invocation of this function
  int mousePosDelta[2] = {x - mousePos[0], y - mousePos[1]};

  switch (controlState)
  {
  // translate the terrain
  case TRANSLATE:
    if (leftMouseButton)
    {
      // control x,y translation via the left mouse button
      terrainTranslate[0] += mousePosDelta[0] * 0.01f;
      terrainTranslate[1] -= mousePosDelta[1] * 0.01f;
    }
    if (middleMouseButton)
    {
      // control z translation via the middle mouse button
      terrainTranslate[2] += mousePosDelta[1] * 0.01f;
    }
    break;

  // rotate the terrain
  case ROTATE:
    if (leftMouseButton)
    {
      // control x,y rotation via the left mouse button
      terrainRotate[0] += mousePosDelta[1];
      terrainRotate[1] += mousePosDelta[0];
    }
    if (middleMouseButton)
    {
      // control z rotation via the middle mouse button
      terrainRotate[2] += mousePosDelta[1];
    }
    break;

  // scale the terrain
  case SCALE:
    if (leftMouseButton)
    {
      // control x,y scaling via the left mouse button
      terrainScale[0] *= 1.0f + mousePosDelta[0] * 0.01f;
      terrainScale[1] *= 1.0f - mousePosDelta[1] * 0.01f;
    }
    if (middleMouseButton)
    {
      // control z scaling via the middle mouse button
      terrainScale[2] *= 1.0f - mousePosDelta[1] * 0.01f;
    }
    break;
  }

  // store the new mouse position
  mousePos[0] = x;
  mousePos[1] = y;
}

// glut mouse motion callback (no button pressed)
void mouseMotionFunc(int x, int y)
{
  // mouse has moved, store the new position
  mousePos[0] = x;
  mousePos[1] = y;
}

// glut mouse button callback
void mouseButtonFunc(int button, int state, int x, int y)
{
  // a mouse button has been pressed or released
  
  // keep track of the mouse button state
  switch (button)
  {
  case GLUT_LEFT_BUTTON:
    leftMouseButton = (state == GLUT_DOWN);
    break;

  case GLUT_MIDDLE_BUTTON:
    middleMouseButton = (state == GLUT_DOWN);
    break;

  case GLUT_RIGHT_BUTTON:
    rightMouseButton = (state == GLUT_DOWN);
    break;
  }

  // determine control state based on modifier keys
  // shifted to use tab instead of ctrl for macbook compatibility
  if (glutGetModifiers() & GLUT_ACTIVE_SHIFT)
  {
    controlState = SCALE;
  }
  else
  {
    controlState = ROTATE;
  }

  // store the new mouse position
  mousePos[0] = x;
  mousePos[1] = y;
}

// glut keyboard callback
void keyboardFunc(unsigned char key, int x, int y)
{
  switch (key)
  {
  case 27:   // ESC key
    exit(0); // exit the program
    break;

  case 'b': // toggle bird's eye view
  case 'B':
    birdEyeView = !birdEyeView;
    std::cout << "bird's eye view: " << (birdEyeView ? "ON" : "OFF") << std::endl;
    break;
    
  case 's': // toggle skybox
  case 'S':
    showSkybox = !showSkybox;
    std::cout << "skybox: " << (showSkybox ? "ON" : "OFF") << std::endl;
    break;
    
  case 'g': // toggle ground plane
  case 'G':
    showGround = !showGround;
    std::cout << "ground plane: " << (showGround ? "ON" : "OFF") << std::endl;
    break;

  case ' ':
    cout << "you pressed the spacebar." << endl;
    break;

  case 'x':
    // take a screenshot
    saveScreenshot("screenshot.jpg");
    break;

  case 't':
    // toggle recording
    isRecording = !isRecording;
    std::cout << "recording: " << (isRecording ? "ON" : "OFF") << std::endl;
    break;
  }
}

//*****************************
// rendering functions
//*****************************

// glut display callback - render the scene
void displayFunc()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // get current matrices
  float modelViewMatrix[16];
  matrix.SetMatrixMode(OpenGLMatrix::ModelView);
  matrix.GetMatrix(modelViewMatrix);

  float projectionMatrix[16];
  matrix.SetMatrixMode(OpenGLMatrix::Projection);
  matrix.GetMatrix(projectionMatrix);

  // set up lighting direction
  float worldLightDir[3] = {1.0f, 1.0f, 0.2f};  // light coming from upper right
  
  // normalize the world light direction
  float length = sqrt(worldLightDir[0]*worldLightDir[0] + 
                     worldLightDir[1]*worldLightDir[1] + 
                     worldLightDir[2]*worldLightDir[2]);
  worldLightDir[0] /= length;
  worldLightDir[1] /= length;
  worldLightDir[2] /= length;

  // transform light direction to view space using the normal matrix
  float normalMatrix[9];
  computeNormalMatrix(modelViewMatrix, normalMatrix);
  
  // draw skybox if enabled
  if (showSkybox) {
    glDepthFunc(GL_LEQUAL);  // change depth function to ensure skybox is drawn at the far plane
    skyboxPipeline->Bind();
    skyboxPipeline->SetUniformVariableMatrix4fv("projectionMatrix", GL_FALSE, projectionMatrix);
    skyboxPipeline->SetUniformVariableMatrix4fv("modelViewMatrix", GL_FALSE, modelViewMatrix);
    skyboxPipeline->SetUniformVariablei("skyboxTexture", 0); // set texture unit 0

    // bind skybox texture
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_CUBE_MAP, skyboxTexture);

    // draw skybox
    skyboxVAO->Bind();
    glDrawElements(GL_TRIANGLES, skyboxIndices.size(), GL_UNSIGNED_INT, 0);
    
    // reset depth function
    glDepthFunc(GL_LESS);
  }

  // draw ground plane with texture if enabled
  if (showGround) {
    texturePipelineProgram->Bind();
    texturePipelineProgram->SetUniformVariableMatrix4fv("modelViewMatrix", GL_FALSE, modelViewMatrix);
    texturePipelineProgram->SetUniformVariableMatrix4fv("projectionMatrix", GL_FALSE, projectionMatrix);

    // bind texture and draw ground
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, groundTextureHandle);
    groundVAO->Bind();
    groundEBO->Bind();
    glDrawElements(GL_TRIANGLES, groundIndices.size(), GL_UNSIGNED_INT, 0);
  }

  // draw spline (rails)
  pipelineProgram->Bind();
  pipelineProgram->SetUniformVariableMatrix4fv("modelViewMatrix", GL_FALSE, modelViewMatrix);
  pipelineProgram->SetUniformVariableMatrix4fv("projectionMatrix", GL_FALSE, projectionMatrix);
  pipelineProgram->SetUniformVariableMatrix3fv("normalMatrix", GL_FALSE, normalMatrix);
  
  // set lighting parameters
  pipelineProgram->SetUniformVariable4f("La", 0.5f, 0.5f, 0.5f, 1.0f);  // ambient light
  pipelineProgram->SetUniformVariable4f("Ld", 0.7f, 0.7f, 0.7f, 1.0f);  // diffuse light
  pipelineProgram->SetUniformVariable4f("Ls", 0.9f, 0.9f, 0.9f, 1.0f);  // specular light
  pipelineProgram->SetUniformVariable3f("viewLightDirection", 0.0f, 0.0f, -1.0f);  // light direction

  // set material properties for rails
  pipelineProgram->SetUniformVariable4f("ka", 0.15f, 0.25f, 0.1f, 1.0f);  // dark greenish-brownish ambient
  pipelineProgram->SetUniformVariable4f("kd", 0.3f, 0.45f, 0.2f, 1.0f);   // dark greenish-brownish diffuse 
  pipelineProgram->SetUniformVariable4f("ks", 0.2f, 0.2f, 0.1f, 1.0f);    // subtle specular
  pipelineProgram->SetUniformVariablef("alpha", 1.0f);                   // sharper highlights

  // draw rails
  splineVAO->Bind();
  railEBO->Bind();
  glDrawElements(GL_TRIANGLES, railIndices.size(), GL_UNSIGNED_INT, 0);

  // draw ties if they exist
  if (!tiePositions.empty() && !tieIndices.empty()) {
    tieVAO->Bind();
    tieEBO->Bind();
    glDrawElements(GL_TRIANGLES, tieIndices.size(), GL_UNSIGNED_INT, 0);
  }

  glutSwapBuffers();
}

//*****************************
// scene setup and initialization
//*****************************

// initialize the scene, including OpenGL state
void initScene(int argc, char *argv[])
{
  // enable depth testing
  glEnable(GL_DEPTH_TEST);
  
  // set white background
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

  // initialize camera parameters
  cameraP = 0.0f;
  cameraSpeed = 0.01f;
  currentSpline = 0;

  // compute spline points and generate roller coaster geometry
  computeCatmullSpline();

  // create and set up basic shader program for rail and tie rendering
  pipelineProgram = new PipelineProgram();
  if (pipelineProgram->BuildShadersFromFiles(shaderBasePath,
                                             "vertexShader.glsl", "fragmentShader.glsl") != 0)
  {
    std::cout << "failed to build pipeline program" << std::endl;
    exit(1);
  }

  // initialize spline VAO/VBOs
  splineVAO = new VAO();
  splineVBO = new VBO(pointPos.size() / 3, 3, pointPos.data(), GL_STATIC_DRAW);
  colorVBO = new VBO(pointCol.size() / 4, 4, pointCol.data(), GL_STATIC_DRAW);
  railEBO = new EBO(railIndices.size(), 1, railIndices.data(), GL_STATIC_DRAW);

  splineVAO->Bind();
  splineVAO->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, splineVBO, "position");
  splineVAO->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, colorVBO, "color");
  railEBO->Bind();

  // initialize tie VAO
  tieVAO = new VAO();
  tieVAO->Bind();
  tieVBOPosition = new VBO(tiePositions.size() / 3, 3, tiePositions.data(), GL_STATIC_DRAW);
  tieVBOColor = new VBO(tieColors.size() / 4, 4, tieColors.data(), GL_STATIC_DRAW);
  tieEBO = new EBO(tieIndices.size(), 1, tieIndices.data(), GL_STATIC_DRAW);
  
  tieVAO->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, tieVBOPosition, "position");
  tieVAO->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, tieVBOColor, "color");
  tieEBO->Bind();

  // initialize ground plane
  groundVAO = new VAO();
  groundVAO->Bind();

  // create VBO for positions
  groundVBOPosition = new VBO(positions.size() / 3, // number of vertices
                              3, // 3 floats per vertex
                              positions.data(),
                              GL_STATIC_DRAW);

  // create VBO for texture coordinates
  groundVBOTexCoord = new VBO(texCoords.size() / 2, // number of texture coordinates
                              2, // 2 floats per vertex
                              texCoords.data(),
                              GL_STATIC_DRAW);

  // create EBO with the correct size calculation
  groundEBO = new EBO(groundIndices.size(), // number of indices
                      1, // 1 uint per index
                      groundIndices.data(),
                      GL_STATIC_DRAW);

  // set up texture pipeline program
  texturePipelineProgram = new PipelineProgram();
  if (texturePipelineProgram->BuildShadersFromFiles(shaderBasePath,
                                                   "textureVertexShader.glsl", "textureFragmentShader.glsl") != 0)
  {
    std::cout << "failed to build texture pipeline program" << std::endl;
    exit(1);
  }

  // connect position attribute
  groundVAO->ConnectPipelineProgramAndVBOAndShaderVariable(
      texturePipelineProgram,
      groundVBOPosition,
      "position");

  // connect texture coordinate attribute
  groundVAO->ConnectPipelineProgramAndVBOAndShaderVariable(
      texturePipelineProgram,
      groundVBOTexCoord,
      "texCoord");

  // make sure EBO is bound
  groundEBO->Bind();

  // load and initialize texture
  glGenTextures(1, &groundTextureHandle);
  
  // debug path information
  char currentDir[1024];
  if (getcwd(currentDir, sizeof(currentDir)) != NULL) {
      printf("current working directory: %s\n", currentDir);
  } else {
      printf("could not get current working directory\n");
  }
  
  // try to load texture
  printf("attempting to load texture from: %s\n", "textures/water_texture.jpg");
  initTexture("textures/water_texture.jpg", groundTextureHandle);

  // initialize skybox
  skyboxPipeline = new PipelineProgram();
  if (skyboxPipeline->BuildShadersFromFiles(shaderBasePath,
                                           "skyboxVertexShader.glsl", 
                                           "skyboxFragmentShader.glsl") != 0)
  {
      std::cout << "failed to build skybox pipeline program" << std::endl;
      exit(1);
  }

  // Debug vertex attributes
  std::cout << "Skybox vertex attribute setup:" << std::endl;
  std::cout << "  - Number of vertices: " << skyboxVertices.size() / 5 << std::endl;
  std::cout << "  - Vertex stride: " << 5 * sizeof(float) << " bytes" << std::endl;
  std::cout << "  - Position attribute: 3 floats at offset 0" << std::endl;
  std::cout << "  - Texture coord attribute: 2 floats at offset " << 3 * sizeof(float) << std::endl;

  // Setup skybox VAO and VBO
  skyboxVAO = new VAO();
  skyboxVAO->Bind();
  skyboxVBO = new VBO(skyboxVertices.size() / 5, 5, skyboxVertices.data(), GL_STATIC_DRAW);

  // Connect the position attribute (first 3 components)
  skyboxVAO->ConnectPipelineProgramAndVBOAndShaderVariable(skyboxPipeline, skyboxVBO, "position");

  // Need to fix the attribute pointers to account for the stride
  GLuint posAttrib = glGetAttribLocation(skyboxPipeline->GetProgramHandle(), "position");
  glVertexAttribPointer(posAttrib, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);

  // Manually setup the texture coordinate attribute (next 2 components)
  GLuint texAttrib = glGetAttribLocation(skyboxPipeline->GetProgramHandle(), "texCoord");
  glEnableVertexAttribArray(texAttrib);
  glVertexAttribPointer(texAttrib, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));

  // Create and bind skybox EBO
  glGenBuffers(1, &skyboxEBO);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, skyboxEBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, skyboxIndices.size() * sizeof(unsigned int), skyboxIndices.data(), GL_STATIC_DRAW);

  // load the cubemap texture
  skyboxTexture = loadCubemap(skyboxFaces);

  // Debug skybox coordinate system
  debugSkyboxCoordinates();
}

// free all allocated resources
void cleanup()
{
  // delete shader programs
  delete pipelineProgram;
  delete texturePipelineProgram;
  delete skyboxPipeline;
  
  // delete VAOs and VBOs for rails
  delete splineVAO;
  delete splineVBO;
  delete colorVBO;
  delete railEBO;
  
  // delete VAOs and VBOs for ground
  delete groundVAO;
  delete groundVBOPosition;
  delete groundVBOTexCoord;
  delete groundEBO;
  
  // delete VAOs and VBOs for ties
  delete tieVAO;
  delete tieVBOPosition;
  delete tieVBOColor;
  delete tieEBO;

  // delete textures
  if (textureID)
  {
    glDeleteTextures(1, &textureID);
  }
  
  if (groundTextureHandle)
  {
    glDeleteTextures(1, &groundTextureHandle);
  }

  // delete skybox
  delete skyboxVAO;
  delete skyboxVBO;
  glDeleteBuffers(1, &skyboxEBO);
  glDeleteTextures(1, &skyboxTexture);
}

//*****************************
// main function
//*****************************

int main(int argc, char *argv[])
{
  // check command line arguments
  if (argc != 2)
  {
    cout << "the arguments are incorrect." << endl;
    cout << "usage: ./hw1 <spline file>" << endl;
    exit(EXIT_FAILURE);
  }

  // load spline data
  loadSpline(argv[1]);
  cout << "loaded spline with " << spline.numControlPoints << " control point(s)." << endl;

  // calculate scene dimensions
  computeLowestPoint();

  // initialize OpenGL and GLUT
  cout << "initializing GLUT..." << endl;
  glutInit(&argc, argv);

  cout << "initializing OpenGL..." << endl;

  // set up display mode based on platform
#ifdef __APPLE__
  glutInitDisplayMode(GLUT_3_2_CORE_PROFILE | GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL);
#else
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL);
#endif

  // create window
  glutInitWindowSize(windowWidth, windowHeight);
  glutInitWindowPosition(0, 0);
  glutCreateWindow(windowTitle);

  // print OpenGL info
  cout << "OpenGL Version: " << glGetString(GL_VERSION) << endl;
  cout << "OpenGL Renderer: " << glGetString(GL_RENDERER) << endl;
  cout << "Shading Language Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;

#ifdef __APPLE__
  // needed on recent Mac OS X versions to correctly display the window
  glutReshapeWindow(windowWidth - 1, windowHeight - 1);
#endif

  // register GLUT callback functions
  glutDisplayFunc(displayFunc);        // rendering
  glutIdleFunc(idleFunc);              // animation
  glutMotionFunc(mouseMotionDragFunc); // mouse movement with button pressed
  glutPassiveMotionFunc(mouseMotionFunc); // mouse movement without button pressed
  glutMouseFunc(mouseButtonFunc);      // mouse button events
  glutReshapeFunc(reshapeFunc);        // window resize
  glutKeyboardFunc(keyboardFunc);      // keyboard events

  // initialize GLEW
#ifdef __APPLE__
  // not needed on Apple
#else
  // Windows, Linux
  GLint result = glewInit();
  if (result != GLEW_OK)
  {
    cout << "error: " << glewGetErrorString(result) << endl;
    exit(EXIT_FAILURE);
  }
#endif

  // perform initialization
  initScene(argc, argv);

  // start main loop
  glutMainLoop();
  
  return 0;
}
