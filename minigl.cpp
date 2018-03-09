/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>

using namespace std;

/**
 * Useful data types
 */
typedef mat<MGLfloat,4,4> mat4; //data structure storing a 4x4 matrix, see mat.h
typedef mat<MGLfloat,3,3> mat3; //data structure storing a 3x3 matrix, see mat.h
typedef vec<MGLfloat,4> vec4;   //data structure storing a 4 dimensional vector, see vec.h
typedef vec<MGLfloat,3> vec3;   //data structure storing a 3 dimensional vector, see vec.h
typedef vec<MGLfloat,2> vec2;   //data structure storing a 2 dimensional vector, see vec.h
struct Vertex{
  vec4 position;
  vec3 color;
};
struct Triangle{
  Vertex A;
  Vertex B;
  Vertex C;
};

/**
 *Global Variables
 */
//mat4 modelview;
//mat4 projection;
vector<mat4> modelview {{0, 0, 0, 0,
                         0, 0, 0, 0,
                         0, 0, 0, 0,
                         0, 0, 0, 0}};

vector<mat4> projection {{0, 0, 0, 0,
                          0, 0, 0, 0,
                          0, 0, 0, 0,
                          0, 0, 0, 0}};


vec3 currColor;
vector<Vertex> currVertices;
vector<Triangle> currTriangles;
MGLpoly_mode currMode;
MGLmatrix_mode currMatrixMode;

MGLfloat nearPlane = -1;
MGLfloat farPlane = 1;

vector<vector<MGLfloat>>zBuffer;

/** 
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}




float Area(vec2 a, vec2 b, vec2 c)
{
  MGLfloat area;
  area = (a[0] * (b[1] - c[1])) + (a[1] * (c[0] - b[0])) + ((b[0] * c[1]) - (b[1] * c[0]));
  return area;
}

void Rasterize_Triangle(const Triangle& tri, int width, int height, MGLpixel* data)
{
  MGLfloat x;
  MGLfloat y;
  MGLfloat i;
  MGLfloat j = 0.0;
  MGLfloat fi;
  MGLfloat fj;
  MGLfloat alpha;
  MGLfloat beta;
  MGLfloat gamma;
  MGLfloat newAlpha;
  MGLfloat newBeta;
  MGLfloat newGamma;
  MGLfloat zValue;
  MGLfloat k;
  int min_x;
  int min_y;
  int max_x;
  int max_y;
  MGLfloat w_a;
  MGLfloat w_b;
  MGLfloat w_c;


  vec3 colors; 




  x = tri.A.position[0];
  y = tri.A.position[1];
  w_a = tri.A.position[3];
  

  i = (x + 1) * 0.5 * width;
  j = (y + 1) * 0.5 * height;
  
  min_x = i;
  min_y = j;
  max_x = i;
  max_y = j;


  fi = i - 0.5;
  fj = j - 0.5;
  
  vec2 A(fi, fj);

  
  x = tri.B.position[0];
  y = tri.B.position[1];
  w_b = tri.B.position[3];

  
  i = (x + 1) * 0.5 * width;
  j = (y + 1) * 0.5 * height;

  if(i < min_x)
    min_x = i;
  if(j < min_y)
    min_y = j;
  if(i > max_x)
    max_x = i;
  if(j > max_y)
    max_y = j;
  
  fi = i - 0.5;
  fj = j - 0.5;
  
  vec2 B(fi, fj);
  
  x = tri.C.position[0];
  y = tri.C.position[1];
  w_c = tri.C.position[3];

  i = (x + 1) * 0.5 * width;
  j = (y + 1) * 0.5 * height;

  if(i < min_x)
    min_x = i;
  if(j < min_y)
    min_y = j;
  if(i > max_x)
    max_x = i;
  if(j > max_y)
    max_y = j;
  

  fi = i - 0.5;
  fj = j - 0.5;
  
  vec2 C(fi, fj);
  



  //cout <<  "Triangle coords: " << tri.A.position << endl << tri.B.position << endl << tri.C.position << endl;
   //cout << "Min x: " << min_x << " Min y:" << min_y << endl;
  
  if(min_x < 0)
    min_x = 0;
  if(min_y < 0)
    min_y = 0;
  if(max_x > width)
    max_x = width;
  if(max_y > height)
    max_y = height;
  for (int i = min_x; i < max_x; i++)
  {
    for(int j = min_y; j < max_y; j++)
    {
      vec2 P(i, j);
      
      //std::cout << Area(A,B,C) << std::endl;

      alpha = Area(P,B,C) / Area(A,B,C);
      beta = Area(A,P,C) / Area(A,B,C);
      gamma = Area(A,B,P) / Area(A,B,C);



      k = (alpha / w_a) + (beta / w_b) + (gamma/ w_c);

      newAlpha = alpha/(w_a * k); 
      newBeta = beta/(w_b * k);
      newGamma = gamma/(w_c * k);
      
       // cout << newAlpha << endl;
       // cout << newBeta << endl;
       // cout << newGamma << endl;
      
      if((newAlpha > 0 && alpha < 1) && 
      (beta > 0 && beta < 1) &&
      (gamma > 0 && gamma < 1))
      {
        //std::cout << "test" << std::endl;
        colors = ((tri.A.color * newAlpha )+ (tri.B.color * newBeta) + (tri.C.color * newGamma)); 
        zValue = ((tri.A.position[2] * alpha )+ (tri.B.position[2] * beta) + (tri.C.position[2] * gamma));

      
        if(zValue <= zBuffer.at(i).at(j) && zValue >= nearPlane && zValue <= farPlane)
        {

          data[i + j * width] = Make_Pixel(255 * colors[0], 255* colors[1], 255* colors[2]);
          zBuffer.at(i).at(j) = zValue;      
        }


      }
      
    }
    
  }

}


/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{

  zBuffer.resize(width);
  for (int i = 0; i < width; i++)
  {
    zBuffer.at(i).resize(height);
  }

  for(int i = 0; i < width; i++)
  {
    for(int j = 0; j < height; j++)
    {
      data[i+j*width] = Make_Pixel(0,0,0);
      zBuffer.at(i).at(j) = 100;
    }
  }
  
  for(int i = 0; i < currTriangles.size(); i++)
  {
    Rasterize_Triangle(currTriangles.at(i), width, height, data);
  }
  
  currTriangles.clear();
}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
  currMode = mode;
}


/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
  if(currMode == MGL_TRIANGLES)
  {
    for (unsigned i = 0; i < currVertices.size(); i+= 3)
    {
      
      Triangle newTriangle;
      newTriangle.A = currVertices.at(i);
      newTriangle.B = currVertices.at(i+1);
      newTriangle.C = currVertices.at(i+2);

      currTriangles.push_back(newTriangle);

      //CASE 1: Triangle inside cube
      // if(abs(newTriangle.A.position.normalized()[0]) < 1 &&
      //    abs(newTriangle.A.position.normalized()[1]) < 1 &&
      //    abs(newTriangle.A.position.normalized()[2]) < 1 &&

      //    abs(newTriangle.B.position.normalized()[0]) < 1 &&
      //    abs(newTriangle.B.position.normalized()[1]) < 1 &&
      //    abs(newTriangle.B.position.normalized()[2]) < 1 &&

      //    abs(newTriangle.C.position.normalized()[0]) < 1 &&
      //    abs(newTriangle.C.position.normalized()[1]) < 1 &&
      //    abs(newTriangle.C.position.normalized()[2]) < 1)
      // {
      //     currTriangles.push_back(newTriangle);
      // }

      // //CASE 2: Triangle outside cube
      // else if(abs(newTriangle.A.position.normalized()[0]) > 1 &&
      //         abs(newTriangle.A.position.normalized()[1]) > 1 &&
      //         abs(newTriangle.A.position.normalized()[2]) > 1 &&

      //         abs(newTriangle.B.position.normalized()[0]) > 1 &&
      //         abs(newTriangle.B.position.normalized()[1]) > 1 &&
      //         abs(newTriangle.B.position.normalized()[2]) > 1 &&

      //         abs(newTriangle.C.position.normalized()[0]) > 1 &&
      //         abs(newTriangle.C.position.normalized()[1]) > 1 &&
      //         abs(newTriangle.C.position.normalized()[2]) > 1)
      // {
      //   //dont draw triangle
      // }

      //CASE 3: One Vertex inside
      
    }
  }
  else if (currMode == MGL_QUADS)
  {
    for (unsigned i = 0; i < currVertices.size(); i+= 4)
    {
  
      Triangle newTriangle1;
      newTriangle1.A = currVertices.at(i);
      newTriangle1.B = currVertices.at(i+1);
      newTriangle1.C = currVertices.at(i+2);
      currTriangles.push_back(newTriangle1);

      
      Triangle newTriangle2;
      newTriangle2.A = currVertices.at(i);
      newTriangle2.B = currVertices.at(i+2);
      newTriangle2.C = currVertices.at(i+3);
      currTriangles.push_back(newTriangle2);
    }
  }
  currVertices.clear();
}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
  mglVertex3(x,y,0);
}

mat4 top_of_active_matrix_stack()
{
    if(currMatrixMode == MGL_MODELVIEW)
    {
      cout << "MODELVIEW" << endl;
        return modelview.back();
    }
    else if(currMatrixMode == MGL_PROJECTION)
    {
        cout << "PROJECTION" << endl;
        return projection.back();
    }
  
}


/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
  vec4 currPos(x, y, z, 1);

  //cout <<  "--------------"  << endl << top_of_active_matrix_stack() << endl << "-------x-------" << endl << projection.back() << endl;

  currPos = projection.back() * modelview.back() * currPos;
  MGLfloat w = currPos[3];
  currPos = currPos / currPos[3];

  currPos[3] = w;

  Vertex currVert;
  currVert.position = currPos;
  currVert.color = currColor;

  //cout << currPos << " " << currColor << endl << 
  //cout << currVert.position << " " << currVert.color << endl;
  
  currVertices.push_back(currVert);
  
  
  //if this doesnt work, make a constructor?
  //model view and projection matrices here
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
  currMatrixMode = mode;
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
  if(currMatrixMode == MGL_MODELVIEW)
  {
      modelview.push_back(modelview.at(modelview.size()-1));
  }
  else if (currMatrixMode == MGL_PROJECTION)
  {
      projection.push_back(projection.at(projection.size()-1));
  }
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
  if(currMatrixMode == MGL_MODELVIEW)
  {
    modelview.pop_back();
  }
  else if(currMatrixMode == MGL_PROJECTION)
  {
    projection.pop_back();
  }
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
  if(currMatrixMode == MGL_MODELVIEW)
  {
    modelview.at(modelview.size()-1) = {{1, 0, 0, 0,
                                         0, 1, 0, 0,
                                         0, 0, 1, 0,
                                         0, 0, 0, 1}};
  }
  else if (currMatrixMode == MGL_PROJECTION)
  {
    projection.at(projection.size()-1) = {{1, 0, 0, 0,
                                           0, 1, 0, 0,
                                           0, 0, 1, 0,
                                           0, 0, 0, 1}};
  }
}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
  mat4 newMatrix = {{matrix[0], matrix[1], matrix[2], matrix[3],
                      matrix[4], matrix[5], matrix[6], matrix[7],
                      matrix[8], matrix[9], matrix[10], matrix[11],
                      matrix[12], matrix[13], matrix[14], matrix[15]}};

  if (currMatrixMode == MGL_MODELVIEW)
  {
      modelview.at(modelview.size()-1) = newMatrix; 
  }
  else if (currMatrixMode == MGL_PROJECTION)
  {
      projection.at(projection.size()-1) = newMatrix;
  }
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{

  mat4 multMatrix = {{matrix[0], matrix[1], matrix[2], matrix[3],
                      matrix[4], matrix[5], matrix[6], matrix[7],
                      matrix[8], matrix[9], matrix[10], matrix[11],
                      matrix[12], matrix[13], matrix[14], matrix[15]}};

  if (currMatrixMode == MGL_MODELVIEW)
  {
      modelview.at(modelview.size()-1) = modelview.back() * multMatrix; 
  }
  else if (currMatrixMode == MGL_PROJECTION)
  {
      projection.at(projection.size()-1) = projection.back() * multMatrix;
  }
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
  // mat4 translate = {{1, 0, 0, x,
  //                    0, 1, 0, y,
  //                    0, 0, 1, z, 
  //                    0, 0, 0, 1}};

  mat4 translate =    {{1, 0, 0, 0,
                        0, 1, 0, 0,
                        0, 0, 1, 0,
                        x, y, z, 1}};

    if(currMatrixMode == MGL_MODELVIEW)
    { 
        modelview.at(modelview.size()-1) = modelview.back() * translate ; 
    }
    else if(currMatrixMode == MGL_PROJECTION)
    {
        projection.at(projection.size()-1) = projection.back() * translate; 
    }
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{

  MGLfloat radianAngle = (angle * 3.14159) / 180.0;
  MGLfloat cosine =  cos(radianAngle);
  MGLfloat sine = sin(radianAngle);

  vec3 unitVec(x,y,z);
  unitVec = unitVec.normalized();

                    //------------column1--------------------------------------
  mat4 rotation = {{cosine + (unitVec[0] * unitVec[0] * (1 - cosine)) , 
                    unitVec[1] * unitVec[0] * (1 - cosine) + unitVec[2] * sine , 
                    unitVec[2] * unitVec[0] * (1 - cosine) - unitVec[1] * sine , 
                    0, 
                    //------------column2---------------------------------------
                    unitVec[0] * unitVec[1] * (1 - cosine) - unitVec[2] * sine ,
                    cosine + (unitVec[1] * unitVec[1] * (1 - cosine)) ,
                    unitVec[2] * unitVec[1] * (1 - cosine) + unitVec[0] * sine , 
                    0, 
                    //------------column3---------------------------------------
                    unitVec[0] * unitVec[2] * (1 - cosine) + unitVec[1] * sine , 
                    unitVec[1] * unitVec[2] * (1 - cosine) - unitVec[0] * sine , 
                    cosine + (unitVec[2] * unitVec[2] * (1 - cosine)) ,
                    0,
                    //-----------column4----------------------------------------
                    0,
                    0,
                    0,
                    1}};

    if(currMatrixMode == MGL_MODELVIEW)
    { 
        modelview.at(modelview.size()-1) = modelview.back() * rotation; 
    }
    else if(currMatrixMode == MGL_PROJECTION)
    {
        projection.at(projection.size()-1) = projection.back() * rotation; 
    }
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
  mat4 scale = {{x, 0, 0, 0,
                 0, y, 0, 0,
                 0, 0, z, 0,
                 0, 0, 0, 1}};

    if(currMatrixMode == MGL_MODELVIEW)
    { 
        modelview.at(modelview.size()-1) = modelview.back() * scale; 
    }
    else if(currMatrixMode == MGL_PROJECTION)
    {
        projection.at(projection.size()-1) = projection.back() * scale; 
    }
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{

    mat4 perspectiveMatrix = {{((2 * near)/(right - left)),0,0,0,
                     0, ((2 * near) / (top - bottom)), 0, 0, 
                     ((right + left) / (right - left)), ((top + bottom)/(top - bottom)), (-(far + near)/(far - near)), -1,
                     0, 0, ((-2*far*near)/(far - near)), 0}};

    if(currMatrixMode == MGL_MODELVIEW)
    {
      modelview.at(modelview.size()-1) = modelview.back() * perspectiveMatrix;
    }
    else if (currMatrixMode == MGL_PROJECTION)
    {
      projection.at(projection.size()-1) = projection.back() * perspectiveMatrix;
    }


  // if(currMatrixMode == MGL_MODELVIEW)
  // {
  //   modelview.push_back({{((2 * near)/(right - left)),0,0,0,
  //                   0, ((2 * near) / (top - bottom)), 0, 0, 
  //                   ((right + left) / (right - left)), ((top + bottom)/(top - bottom)), (-(far + near)/(far - near)), -1,
  //                   0, 0, ((-2*far*near)/(far - near)), 0}});  
  // }
  // else if (currMatrixMode == MGL_PROJECTION)
  // {
  //  projection.push_back({{((2 * near)/(right - left)),0,0,0,
  //                   0, ((2 * near) / (top - bottom)), 0, 0, 
  //                   ((right + left) / (right - left)), ((top + bottom)/(top - bottom)), (-(far + near)/(far - near)), -1,
  //                   0, 0, ((-2*far*near)/(far - near)), 0}}); 
  // }
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{


    mat4 orthoMatrix = {{2/(right-left), 0, 0, 0,
                        0, 2/(top-bottom), 0, 0,
                        0, 0, -2/(far-near), 0, 
                        -(right + left)/(right-left), -(top+bottom)/(top-bottom), -(far + near)/(far-near), 1}};

    if(currMatrixMode == MGL_MODELVIEW)
    {
      modelview.at(modelview.size()-1) = modelview.back() * orthoMatrix;
    }
    else if (currMatrixMode == MGL_PROJECTION)
    {
      //cout << projection.back() << endl;
      projection.at(projection.size()-1) = projection.back() * orthoMatrix;
    }

  // if(currMatrixMode == MGL_MODELVIEW)
  // {
  //   modelview.push_back({{2/(right-left), 0, 0, 0,
  //                0, 2/(top-bottom), 0, 0,
  //                0, 0, -2/(far-near), 0, 
  //               -(right + left)/(right-left), -(top+bottom)/(top-bottom), -(far + near)/(far-near), 1}});  
  // }
  // else if (currMatrixMode == MGL_PROJECTION)
  // {
  // projection.push_back({{2/(right-left), 0, 0, 0,
  //                0, 2/(top-bottom), 0, 0,
  //                0, 0, -2/(far-near), 0, 
  //               -(right + left)/(right-left), -(top+bottom)/(top-bottom), -(far + near)/(far-near), 1}});
  // }

    
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
    
  vec3 color(red, green, blue);
  currColor = color;
  
  //curr_color = glColor3f(red, green, blue);
  
}
