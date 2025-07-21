#include "ebo.h"

#include <cstring>
#include <cstdio>
#include <iostream>

using namespace std;

EBO::EBO(int numVertices_, int numFloatsPerVertex_, unsigned int * data,  const GLenum usage) : numVertices(numVertices_), numFloatsPerVertex(numFloatsPerVertex_)
{
  // Create the VBO handle 
  glGenBuffers(1, &handle);

  // Initialize the VBO.
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, handle);
  const int numBytes = numVertices * numFloatsPerVertex * sizeof(float);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, numBytes, data, usage);
}

EBO::~EBO()
{
  // Delete the VBO.
  glDeleteBuffers(1, &handle);
}

void EBO::Bind()
{
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, handle);
}