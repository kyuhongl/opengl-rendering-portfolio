/*
  CSCI 420 Computer Graphics, Computer Science, USC
  Assignment 1: Height Fields with Shaders.
  C/C++ starter code

  Student username: kyuhongl
*/

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

void cleanup();

int mousePos[2]; // x,y screen coordinates of the current mouse position

int leftMouseButton = 0; // 1 if pressed, 0 if not 
int middleMouseButton = 0; // 1 if pressed, 0 if not
int rightMouseButton = 0; // 1 if pressed, 0 if not

typedef enum { ROTATE, TRANSLATE, SCALE } CONTROL_STATE;
CONTROL_STATE controlState = ROTATE;

// Transformations of the terrain.
float terrainRotate[3] = { 0.0f, 0.0f, 0.0f }; 
// terrainRotate[0] gives the rotation around x-axis (in degrees)
// terrainRotate[1] gives the rotation around y-axis (in degrees)
// terrainRotate[2] gives the rotation around z-axis (in degrees)
float terrainTranslate[3] = { 0.0f, 0.0f, 0.0f };
float terrainScale[3] = { 1.0f, 1.0f, 1.0f };

// extra credit
float colorCycle = 0.0f;
bool rainbowMode = false;

// Width and height of the OpenGL window, in pixels.
int windowWidth = 1280;
int windowHeight = 720;
char windowTitle[512] = "CSCI 420 Homework 1";

// Stores the image loaded from disk.
ImageIO * heightmapImage;

// Number of vertices in the single triangle (starter code).
int numVertices;

// textures
GLuint textureID;
2
// recording
int frameCount = 0;
const int MAX_FRAMES = 150;
bool isRecording = false;
const char* outputDirectory = "frames/";

// CSCI 420 helper classes.
OpenGLMatrix matrix;

PipelineProgram * pipelineProgram = nullptr;

VBO * vboVertices = nullptr;
VBO * vboColors = nullptr;
VBO * vboCenterVertices = nullptr;
VBO * vboLeftVertices = nullptr;
VBO * vboRightVertices = nullptr;
VBO * vboUpVertices = nullptr;
VBO * vboDownVertices = nullptr;

EBO * eboWireframe = nullptr;
EBO * eboTriangles = nullptr;

VAO * vao = nullptr;

VBO * vboTexCoords = nullptr;

float scale = 1.0f;
float exponent = 1.0f;
int mode = 0;
int renderType = 0;
int numWireframeIndices;
int numTriangleIndices;

bool isTabActive = false;
// Write a screenshot to the specified filename.
void saveScreenshot(const char * filename)
{
  unsigned char * screenshotData = new unsigned char[windowWidth * windowHeight * 3];
  glReadPixels(0, 0, windowWidth, windowHeight, GL_RGB, GL_UNSIGNED_BYTE, screenshotData);

  ImageIO screenshotImg(windowWidth, windowHeight, 3, screenshotData);

  if (screenshotImg.save(filename, ImageIO::FORMAT_JPEG) == ImageIO::OK)
    cout << "File " << filename << " saved successfully." << endl;
  else cout << "Failed to save file " << filename << '.' << endl;

  delete [] screenshotData;
}

void createDirectory(const char* path) {
    mkdir(path, 0777);
}

void loadTexture(const char* filename) {
    ImageIO* textureImage = new ImageIO();
    if (textureImage->loadJPEG(filename) != ImageIO::OK) {
        cout << "Error loading texture image." << endl;
        return;
    }

    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_2D, textureID);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 
                 textureImage->getWidth(), textureImage->getHeight(), 
                 0, GL_RGB, GL_UNSIGNED_BYTE, textureImage->getPixels());

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    delete textureImage;
}

void idleFunc()
{
  // Do some stuff... 
  // For example, here, you can save the screenshots to disk (to make the animation).

  // EXTRA CREIT TEXTURE MODE
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, textureID);
  GLint useTexLoc = glGetUniformLocation(pipelineProgram->GetProgramHandle(), "useTexture");
  glUniform1i(useTexLoc, renderType == 5); 
 
  if (isRecording && frameCount < MAX_FRAMES) {
      // Create filename with leading zeros
      char filename[256];
      createDirectory(outputDirectory);
      sprintf(filename, "%s%03d.jpg", outputDirectory, frameCount);
      

      // ROTATION 
      terrainRotate[1] += 2.0f; // Rotate around Y-axis
      
      

      saveScreenshot(filename);
      
      frameCount++;
      
      if (frameCount >= MAX_FRAMES) {
          isRecording = false;
          cout << "Animation recording complete!" << endl;
      }

  }
  // RAINBOW MODE
  if (rainbowMode) {
    
    colorCycle += 0.01f;
    if (colorCycle > 1.0f)
    {
      colorCycle -= 1.0f;
    }
    pipelineProgram->Bind();
    pipelineProgram->SetUniformVariablef("colorCycle", colorCycle);
  }

  // Notify GLUT that it should call displayFunc.
  glutPostRedisplay();
}

void reshapeFunc(int w, int h)
{
  glViewport(0, 0, w, h);

  // When the window has been resized, we need to re-set our projection matrix.
  matrix.SetMatrixMode(OpenGLMatrix::Projection);
  matrix.LoadIdentity();
  // You need to be careful about setting the zNear and zFar. 
  // Anything closer than zNear, or further than zFar, will be culled.
  const float zNear = 0.1f;
  const float zFar = 10000.0f;
  const float humanFieldOfView = 60.0f;
  matrix.Perspective(humanFieldOfView, 1.0f * w / h, zNear, zFar);
}

void mouseMotionDragFunc(int x, int y)
{
  // Mouse has moved, and one of the mouse buttons is pressed (dragging).
  cout << "Dragging. x: " << x << "y:" << y << endl;

  // the change in mouse position since the last invocation of this function
  int mousePosDelta[2] = { x - mousePos[0], y - mousePos[1] };

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

void mouseMotionFunc(int x, int y)
{
  // Mouse has moved.
  // Store the new mouse position.
  mousePos[0] = x;
  mousePos[1] = y;
}

void mouseButtonFunc(int button, int state, int x, int y)
{
  // A mouse button has has been pressed or depressed.

  // Keep track of the mouse button state, in leftMouseButton, middleMouseButton, rightMouseButton variables.
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

  // Keep track of whether TAB and SHIFT keys are pressed.
  // CTRl changed to TAB because im on a macbook
  if (isTabActive)
  {
    controlState = TRANSLATE;
  }
  else if (glutGetModifiers() & GLUT_ACTIVE_SHIFT)
  {
    controlState = SCALE;
  }
  else
  {
    controlState = ROTATE;
  }

  // Store the new mouse position.
  mousePos[0] = x;
  mousePos[1] = y;
}

void keyboardFunc(unsigned char key, int x, int y)
{
    switch (key)
    {
        case 27: // ESC 
            cleanup(); 
            exit(0);    
        break;

    case 'a':
    case 'A':
        if (!isRecording) {
        cout << "Starting Recording." << endl;
        frameCount = 0;
        isRecording = true;
        }
    break;
    
    // POINT MODE
    case '1':
        renderType = 1;
        mode = 0;
        pipelineProgram->SetUniformVariablei("mode", mode);
    break;

    // WIREFRAME MODE
    case '2':
        renderType = 2;
        mode = 0;
        pipelineProgram->SetUniformVariablei("mode", mode);
    break;

    // TRIANGLE MODE
    case '3':
        renderType = 3;
        mode = 0;
        pipelineProgram->SetUniformVariablei("mode", mode);
    break;

    // SMOOTHING MODE
    case '4':
        renderType = 3;
        mode = 1;
        pipelineProgram->SetUniformVariablei("mode", mode);
    break;

    // TEXTURES 5-7
    case '5': 
        loadTexture("texture1.jpeg");
        renderType = 5;
        mode = 0;
        pipelineProgram->SetUniformVariablei("mode", mode);
        pipelineProgram->SetUniformVariablei("useTexture", 1);  
    break;

    case '6': 
        loadTexture("texture2.jpeg");
        renderType = 5;
        mode = 0;
        pipelineProgram->SetUniformVariablei("mode", mode);
        pipelineProgram->SetUniformVariablei("useTexture", 1);  
    break;

    case '7': 
        loadTexture("texture3.jpeg");
        renderType = 5;
        mode = 0;
        pipelineProgram->SetUniformVariablei("mode", mode);
        pipelineProgram->SetUniformVariablei("useTexture", 1);  

    break;

    // SCALE INCREASE
    case ('='):
        pipelineProgram->Bind();
        scale *= 2.0f;
        cout << " PRESSED" << scale << endl;
        pipelineProgram->SetUniformVariablef("scale", scale);
    break;

    // SCALE DECREASE
    case '-':
        scale /= 2.0f;
        pipelineProgram->SetUniformVariablef("scale", scale);
    break;

    // EXPONENT INCREASE
    case '9':
        exponent *= 2.0f;
        pipelineProgram->SetUniformVariablef("exponent", exponent);
    break;

    // EXPONENT DECREASE
    case '0':
        exponent /= 2.0f;
        pipelineProgram->SetUniformVariablef("exponent", exponent);
    break;
    
    // TAB
    case '\t':
      isTabActive = true;
    break;

    // RAINBOW MODE
    case 'r':
    case 'R':
      rainbowMode = !rainbowMode;
      pipelineProgram->SetUniformVariablei("rainbowMode", rainbowMode ? 1 : 0);
    break;

    // SCREENSHOT
    case 'x':
      // Take a screenshot.
      saveScreenshot("screenshot.jpg");
    break;
  }
}
// key being released function for TAB
void keyboardUpFunc(unsigned char key, int x, int y)
{
  switch (key)
  {
    case '\t': // Tab key
      isTabActive = false;
    break;
  }
}

void displayFunc()
{
  // This function performs the actual rendering.

  // First, clear the screen.
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Set up the camera position, focus point, and the up vector.
  matrix.SetMatrixMode(OpenGLMatrix::ModelView);
  matrix.LoadIdentity();
  matrix.LookAt(0.0, 1.0, 2.0,
                0.0, 0.0, 0.0,
                0.0, 1.0, 0.0);

  // In here, you can do additional modeling on the object, such as performing translations, rotations and scales.
  // ...
  matrix.Translate(terrainTranslate[0], terrainTranslate[1], terrainTranslate[2]);
      matrix.Rotate(terrainRotate[0], 1.0f, 0.0f, 0.0f);
      matrix.Rotate(terrainRotate[1], 0.0f, 1.0f, 0.0f); 
      matrix.Rotate(terrainRotate[2], 0.0f, 0.0f, 1.0f); 
  matrix.Scale(terrainScale[0], terrainScale[1], terrainScale[2]);

  // Read the current modelview and projection matrices from our helper class.
  // The matrices are only read here; nothing is actually communicated to OpenGL yet.
  float modelViewMatrix[16];
  matrix.SetMatrixMode(OpenGLMatrix::ModelView);
  matrix.GetMatrix(modelViewMatrix);

  float projectionMatrix[16];
  matrix.SetMatrixMode(OpenGLMatrix::Projection);
  matrix.GetMatrix(projectionMatrix);
 

  // Upload the modelview and projection matrices to the GPU. Note that these are "uniform" variables.
  // Important: these matrices must be uploaded to *all* pipeline programs used.
  // In hw1, there is only one pipeline program, but in hw2 there will be several of them.
  // In such a case, you must separately upload to *each* pipeline program.
  // Important: do not make a typo in the variable name below; otherwise, the program will malfunction.
    pipelineProgram->Bind();
    pipelineProgram->SetUniformVariableMatrix4fv("modelViewMatrix", GL_FALSE, modelViewMatrix);
    pipelineProgram->SetUniformVariableMatrix4fv("projectionMatrix", GL_FALSE, projectionMatrix);
    pipelineProgram->SetUniformVariablei("mode", mode);
    pipelineProgram->SetUniformVariablef("scale", scale);
    pipelineProgram->SetUniformVariablef("exponent", exponent);

  vao->Bind();
  // Execute the rendering.
  // Bind the VAO that we want to render. Remember, one object = one VAO. 
  if (renderType == 1) {
    glPointSize(1.0f);
    glDrawArrays(GL_POINTS, 0, numVertices); 
  } else if (renderType == 2) {
    glLineWidth(1.0f);
    eboWireframe->Bind();
    glDrawElements(GL_LINES, numWireframeIndices, GL_UNSIGNED_INT, 0);
     
  } else if (renderType == 3) {
    eboTriangles->Bind();
     glDrawElements(GL_TRIANGLES, numTriangleIndices, GL_UNSIGNED_INT, 0); 
     glGetError();
  } else if (renderType == 5) {
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, textureID);
    eboTriangles->Bind();
    glDrawElements(GL_TRIANGLES, numTriangleIndices, GL_UNSIGNED_INT, 0);
  } else {
     glDrawArrays(GL_POINTS, 0, numVertices); 
  }
  // Swap the double-buffers.
  glutSwapBuffers();
}

void initScene(int argc, char *argv[])
{
  // Load the image from a jpeg disk file into main memory.
  heightmapImage = new ImageIO();
  if (heightmapImage->loadJPEG(argv[1]) != ImageIO::OK)
  {
    cout << "Error reading image " << argv[1] << "." << endl;
    exit(EXIT_FAILURE);
  }

  // Set the background color.
  glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black color.

  // Enable z-buffering (i.e., hidden surface removal using the z-buffer algorithm).
  glEnable(GL_DEPTH_TEST);

  // Create a pipeline program. This operation must be performed BEFORE we initialize any VAOs.
  // A pipeline program contains our shaders. Different pipeline programs may contain different shaders.
  // In this homework, we only have one set of shaders, and therefore, there is only one pipeline program.
  // In hw2, we will need to shade different objects with different shaders, and therefore, we will have
  // several pipeline programs (e.g., one for the rails, one for the ground/sky, etc.).
  pipelineProgram = new PipelineProgram(); // Load and set up the pipeline program, including its shaders.
  // Load and set up the pipeline program, including its shaders.
  if (pipelineProgram->BuildShadersFromFiles(shaderBasePath, "vertexShader.glsl", "fragmentShader.glsl") != 0)
  {
    cout << "Failed to build the pipeline program." << endl;
    throw 1;
  } 
  cout << "Successfully built the pipeline program." << endl;
    
  // Bind the pipeline program that we just created. 
  // The purpose of binding a pipeline program is to activate the shaders that it contains, i.e.,
  // any object rendered from that point on, will use those shaders.
  // When the application starts, no pipeline program is bound, which means that rendering is not set up.
  // So, at some point (such as below), we need to bind a pipeline program.
  // From that point on, exactly one pipeline program is bound at any moment of time.
  pipelineProgram->Bind();
  pipelineProgram->SetUniformVariablef("colorCycle", 0.0f);
  pipelineProgram->SetUniformVariablei("rainbowMode", false);

  // Prepare the triangle position and color data for the VBO. 
  // The code below sets up a single triangle (3 vertices).
  // The triangle will be rendered using GL_TRIANGLES (in displayFunc()).

  int img_width = heightmapImage->getWidth();
  int img_height = heightmapImage->getHeight();
  
  numVertices = (img_width) * (img_height); // This must be a global variable, so that we know how many vertices to render in glDrawArrays.

  // Vertex positions.
  float * positions = (float*) malloc (numVertices * 3 * sizeof(float));
  float * p_center = (float*) malloc (numVertices * 3 * sizeof(float));
  float * p_left = (float*) malloc (numVertices * 3 * sizeof(float));
  float * p_right = (float*) malloc (numVertices * 3 * sizeof(float));
  float * p_up = (float*) malloc (numVertices * 3 * sizeof(float));
  float * p_down = (float*) malloc (numVertices * 3 * sizeof(float));
  float * colors = (float*) malloc (numVertices * 4 * sizeof(float));

  float heightScale = 0.5f;
  
  // NORMAL POINT/POSITION LOGIC
  int loc = 0;
  for (int height = 0; height < img_height; height++) {
      for (int width = 0; width < img_width; width++) {
          float heightValue = heightmapImage->getPixel(width, height, 0) / 255.0f;
          
          float nx = (float)width / (img_width - 1) * 2.0f - 1.0f;
          float nz = (float)height / (img_height - 1) * 2.0f - 1.0f;

          positions[loc * 3 + 0] = nx;
          positions[loc * 3 + 1] = heightValue * heightScale;
          positions[loc * 3 + 2] = nz;
          
          colors[loc * 4 + 0] = heightValue;
          colors[loc * 4 + 1] = heightValue;
          colors[loc * 4 + 2] = heightValue;
          colors[loc * 4 + 3] = 1.0f;
          
          loc++;
      }
  }

  // TEXTURE VERTEX LOGIC: Same as before, but with 2 indices each 
  float* texCoords = (float*)malloc(numVertices * 2 * sizeof(float));
  int texLoc = 0;
  for (int height = 0; height < img_height; height++) {
      for (int width = 0; width < img_width; width++) {
          texCoords[texLoc * 2 + 0] = (float)width / (img_width - 1);  
          texCoords[texLoc * 2 + 1] = (float)height / (img_height - 1); 
          texLoc++;
      }
  }


  // LEFT RIGHT UP DOWN NEIGHBOR POINTS AND AVERAGE POINT LOGIC (FOR SMOOTHING)
  loc = 0;
  for (int height = 0; height < img_height; height++) {
    for (int width = 0; width < img_width; width++) {
        float heightValue = heightmapImage->getPixel(width, height, 0) / 255.0f;
        
        float nx = (float)width / (img_width - 1) * 2.0f - 1.0f;
        float nz = (float)height / (img_height - 1) * 2.0f - 1.0f;

        // set positions
        positions[loc * 3 + 0] = nx;
        positions[loc * 3 + 1] = heightValue * heightScale;
        positions[loc * 3 + 2] = nz;

        // set centers (same as position)
        p_center[loc * 3 + 0] = nx;
        p_center[loc * 3 + 1] = heightValue;
        p_center[loc * 3 + 2] = nz;

        // set left neighbors
        float leftHeight = (width > 0) ? heightmapImage->getPixel(width-1, height, 0) / 255.0f : heightValue;
        p_left[loc * 3 + 0] = (width > 0) ? (float)(width-1) / (img_width - 1) * 2.0f - 1.0f : nx;
        p_left[loc * 3 + 1] = leftHeight;
        p_left[loc * 3 + 2] = nz;

        // set right neighbors
        float rightHeight = (width < img_width-1) ? heightmapImage->getPixel(width+1, height, 0) / 255.0f : heightValue;
        p_right[loc * 3 + 0] = (width < img_width-1) ? (float)(width+1) / (img_width - 1) * 2.0f - 1.0f : nx;
        p_right[loc * 3 + 1] = rightHeight;
        p_right[loc * 3 + 2] = nz;

        // set up neighbors
        float upHeight = (height > 0) ? heightmapImage->getPixel(width, height-1, 0) / 255.0f : heightValue;
        p_up[loc * 3 + 0] = nx;
        p_up[loc * 3 + 1] = upHeight;
        p_up[loc * 3 + 2] = (height > 0) ? (float)(height-1) / (img_height - 1) * 2.0f - 1.0f : nz;

        // set down neighbors
        float downHeight = (height < img_height-1) ? heightmapImage->getPixel(width, height+1, 0) / 255.0f : heightValue;
        p_down[loc * 3 + 0] = nx;
        p_down[loc * 3 + 1] = downHeight;
        p_down[loc * 3 + 2] = (height < img_height-1) ? (float)(height+1) / (img_height - 1) * 2.0f - 1.0f : nz;

        // set colors
        colors[loc * 4 + 0] = heightValue;
        colors[loc * 4 + 1] = heightValue;
        colors[loc * 4 + 2] = heightValue;
        colors[loc * 4 + 3] = 1.0f;
        
        loc++;
    }
  }

  // WIREFRAME MODE 

  numWireframeIndices = 2 * ((img_width - 1) * img_height + (img_height - 1) * img_width);
  unsigned int* wireframeIndices = (unsigned int*) malloc(numWireframeIndices * sizeof(unsigned int));

  loc = 0;
  
  // horizontal lines
  for (int h = 0; h < img_height; h++) {
      for (int w = 0; w < img_width - 1; w++) {
          wireframeIndices[loc++] = h * img_width + w;
          wireframeIndices[loc++] = h * img_width + (w + 1);
      }
  }

  // vertical lines
  for (int w = 0; w < img_width; w++) {
      for (int h = 0; h < img_height - 1; h++) {
          wireframeIndices[loc++] = h * img_width + w;
          wireframeIndices[loc++] = (h + 1) * img_width + w;
      }
  }

// TRIANGLE MODE
  numTriangleIndices = (img_width-1) * (img_height-1) * 6;
  unsigned int* triangleIndices = (unsigned int*) malloc(numTriangleIndices * sizeof(unsigned int));

  int index = 0;
  for (int z = 0; z < img_height-1; z++) {
      for (int x = 0; x < img_width-1; x++) {
          // calculate vertex indices for this grid cell, which should just be 4 points of a square.
          unsigned int topLeft = z * img_width + x;
          unsigned int topRight = topLeft + 1;
          unsigned int bottomLeft = (z + 1) * img_width + x;
          unsigned int bottomRight = bottomLeft + 1;

          // first triangle of square
          triangleIndices[index++] = topLeft;
          triangleIndices[index++] = bottomLeft;
          triangleIndices[index++] = topRight;

          // second triangle of square
          triangleIndices[index++] = topRight;
          triangleIndices[index++] = bottomLeft;
          triangleIndices[index++] = bottomRight;

      }
  }
 



  // Create the VBOs. 
  // We make a separate VBO for vertices and colors. 
  // This operation must be performed BEFORE we initialize any VAOs.
  vboVertices = new VBO(numVertices, 3, positions, GL_STATIC_DRAW); // 3 values per position
  vboCenterVertices = new VBO(numVertices, 3, p_center, GL_STATIC_DRAW); 
  vboLeftVertices = new VBO(numVertices, 3, p_left, GL_STATIC_DRAW); 
  vboRightVertices = new VBO(numVertices, 3, p_right, GL_STATIC_DRAW); 
  vboUpVertices = new VBO(numVertices, 3, p_up, GL_STATIC_DRAW); 
  vboDownVertices = new VBO(numVertices, 3, p_down, GL_STATIC_DRAW); 
  vboColors = new VBO(numVertices, 4, colors, GL_STATIC_DRAW); // 4 values per color

  eboWireframe = new EBO(numWireframeIndices, 1, wireframeIndices, GL_STATIC_DRAW);
  eboTriangles = new EBO(numTriangleIndices, 1, triangleIndices, GL_STATIC_DRAW);



  // Create the VAOs. There is a single VAO in this example.
  // Important: this code must be executed AFTER we created our pipeline program, and AFTER we set up our VBOs.
  // A VAO contains the geometry for a single object. There should be one VAO per object.
  // In this homework, "geometry" means vertex positions and colors. In homework 2, it will also include
  // vertex normal and vertex texture coordinates for texture mapping.
  vao = new VAO();

  // Set up the relationship between the "position" shader variable and the VAO.
  // Important: any typo in the shader variable name will lead to malfunction.
  vao->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboVertices, "position");
  vao->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboCenterVertices, "p_center");
  vao->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboLeftVertices, "p_left");
  vao->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboRightVertices, "p_right");
  vao->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboUpVertices, "p_up");
  vao->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboDownVertices, "p_down");
  


  eboWireframe->Bind();
  eboTriangles->Bind();

  // Set up the relationship between the "color" shader variable and the VAO.
  // Important: any typo in the shader variable name will lead to malfunction.
  vao->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboColors, "color");

  // texture
  vboTexCoords = new VBO(numVertices, 2, texCoords, GL_STATIC_DRAW);
  vao->ConnectPipelineProgramAndVBOAndShaderVariable(pipelineProgram, vboTexCoords, "texCoord");

  pipelineProgram->Bind();

  texLoc = glGetUniformLocation(pipelineProgram->GetProgramHandle(), "texture");
  glUniform1i(texLoc, 0); 
  
  GLint useTexLoc = glGetUniformLocation(pipelineProgram->GetProgramHandle(), "textureSampler");
  glUniform1i(useTexLoc, 0); 
  
  // We don't need this data any more, as we have already uploaded it to the VBO. And so we can destroy it, to avoid a memory leak.
  free(positions);
  free(p_center);
  free(p_left);
  free(p_right);
  free(p_up);
  free(p_down);
  free(colors);
  free(wireframeIndices);
  free(triangleIndices);
  free(texCoords);

}

// Just in case ?
void cleanup() {
    delete pipelineProgram;
    delete vboVertices;
    delete vboColors;
    delete vboCenterVertices;
    delete vboLeftVertices;
    delete vboRightVertices;
    delete vboUpVertices;
    delete vboDownVertices;
    delete eboWireframe;
    delete eboTriangles;
    delete vao;
    delete heightmapImage;
    
    if (textureID) {
        glDeleteTextures(1, &textureID);
    }
}

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    cout << "The arguments are incorrect." << endl;
    cout << "usage: ./hw1 <heightmap file>" << endl;
    exit(EXIT_FAILURE);
  }

  cout << "Initializing GLUT..." << endl;
  glutInit(&argc,argv);

  cout << "Initializing OpenGL..." << endl;

  #ifdef __APPLE__
    glutInitDisplayMode(GLUT_3_2_CORE_PROFILE | GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL);
  #else
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL);
  #endif

  glutInitWindowSize(windowWidth, windowHeight);
  glutInitWindowPosition(0, 0);  
  glutCreateWindow(windowTitle);


  cout << "OpenGL Version: " << glGetString(GL_VERSION) << endl;
  cout << "OpenGL Renderer: " << glGetString(GL_RENDERER) << endl;
  cout << "Shading Language Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;

  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(windowWidth - 1, windowHeight - 1);
  #endif

  // Tells GLUT to use a particular display function to redraw.
  glutDisplayFunc(displayFunc);
  // Perform animation inside idleFunc.
  glutIdleFunc(idleFunc);
  // callback for mouse drags
  glutMotionFunc(mouseMotionDragFunc);
  // callback for idle mouse movement
  glutPassiveMotionFunc(mouseMotionFunc);
  // callback for mouse button changes
  glutMouseFunc(mouseButtonFunc);
  // callback for resizing the window
  glutReshapeFunc(reshapeFunc);
  // callback for pressing the keys on the keyboard
  glutKeyboardFunc(keyboardFunc);
  // new callback for releasing a key on the keyboard
  glutKeyboardUpFunc(keyboardUpFunc);

  // init glew
  #ifdef __APPLE__
    // nothing is needed on Apple
  #else
    // Windows, Linux
    GLint result = glewInit();
    if (result != GLEW_OK)
    {
      cout << "error: " << glewGetErrorString(result) << endl;
      exit(EXIT_FAILURE);
    }
  #endif

  // Perform the initialization.
  initScene(argc, argv);

  // Sink forever into the GLUT loop.
  glutMainLoop();
}

