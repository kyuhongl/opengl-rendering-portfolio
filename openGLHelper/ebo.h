#ifndef _EBO_H_
#define _EBO_H_

#include "openGLHeader.h"


class EBO
{
public:

  // Initialize the VBO.
  // "numVertices": the number of vertices
  // "numFloatsPerVertex": the number of floating point values per vertex; e.g. 3 for vertex positions, 4 for colors, and 2 for texture coordinates.
  // "usage" must be either GL_STATIC_DRAW or GL_DYNAMIC_DRAW
  EBO(int numVertices, int numFloatsPerVertex, unsigned int * data, const GLenum usage = GL_STATIC_DRAW);
  virtual ~EBO();

  // Binds (activates) this VBO.
  void Bind();

  // Get handle to this VBO.
  GLuint GetHandle() { return handle; }
  // Get the number of vertices in this VBO.
  int GetNumVertices() { return numVertices; }
  // Get the number of floating point values per vertex in this VBO.
  int GetNumFloatsPerVertex() { return numFloatsPerVertex; }

protected:
  GLuint handle; // the handle to the VBO
  int numVertices;
  int numFloatsPerVertex;
};

#endif

