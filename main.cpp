#include <Novice.h>
#include <cassert>
#include <cmath>
#define _USE_MATH_DEFINES
#include <algorithm>
#include <imgui.h>
#include <math.h>

const char kWindowTitle[] = "LE2B_09_コジマユウヤ";

struct Vector3 {
  float x;
  float y;
  float z;
};

struct Matrix4x4 {
  float m[4][4];
};

struct Sphere {
  Vector3 center;
  float radius;
};

struct Plane {
  Vector3 normal;
  float distance;
};

struct Segment {
  Vector3 origin;
  Vector3 diff;
};

struct Triangle {
  Vector3 vertices[3]; // 頂点
};

struct AABB {
  Vector3 min; // 最小点
  Vector3 max; // 最大点
};

// 行列の積
Matrix4x4 Multiply(Matrix4x4 matrix1, Matrix4x4 matrix2) {

  Matrix4x4 result;

  result.m[0][0] =
      matrix1.m[0][0] * matrix2.m[0][0] + matrix1.m[0][1] * matrix2.m[1][0] +
      matrix1.m[0][2] * matrix2.m[2][0] + matrix1.m[0][3] * matrix2.m[3][0];

  result.m[0][1] =
      matrix1.m[0][0] * matrix2.m[0][1] + matrix1.m[0][1] * matrix2.m[1][1] +
      matrix1.m[0][2] * matrix2.m[2][1] + matrix1.m[0][3] * matrix2.m[3][1];

  result.m[0][2] =
      matrix1.m[0][0] * matrix2.m[0][2] + matrix1.m[0][1] * matrix2.m[1][2] +
      matrix1.m[0][2] * matrix2.m[2][2] + matrix1.m[0][3] * matrix2.m[3][2];

  result.m[0][3] =
      matrix1.m[0][0] * matrix2.m[0][3] + matrix1.m[0][1] * matrix2.m[1][3] +
      matrix1.m[0][2] * matrix2.m[2][3] + matrix1.m[0][3] * matrix2.m[3][3];

  result.m[1][0] =
      matrix1.m[1][0] * matrix2.m[0][0] + matrix1.m[1][1] * matrix2.m[1][0] +
      matrix1.m[1][2] * matrix2.m[2][0] + matrix1.m[1][3] * matrix2.m[3][0];

  result.m[1][1] =
      matrix1.m[1][0] * matrix2.m[0][1] + matrix1.m[1][1] * matrix2.m[1][1] +
      matrix1.m[1][2] * matrix2.m[2][1] + matrix1.m[1][3] * matrix2.m[3][1];

  result.m[1][2] =
      matrix1.m[1][0] * matrix2.m[0][2] + matrix1.m[1][1] * matrix2.m[1][2] +
      matrix1.m[1][2] * matrix2.m[2][2] + matrix1.m[1][3] * matrix2.m[3][2];

  result.m[1][3] =
      matrix1.m[1][0] * matrix2.m[0][3] + matrix1.m[1][1] * matrix2.m[1][3] +
      matrix1.m[1][2] * matrix2.m[2][3] + matrix1.m[1][3] * matrix2.m[3][3];

  result.m[2][0] =
      matrix1.m[2][0] * matrix2.m[0][0] + matrix1.m[2][1] * matrix2.m[1][0] +
      matrix1.m[2][2] * matrix2.m[2][0] + matrix1.m[2][3] * matrix2.m[3][0];

  result.m[2][1] =
      matrix1.m[2][0] * matrix2.m[0][1] + matrix1.m[2][1] * matrix2.m[1][1] +
      matrix1.m[2][2] * matrix2.m[2][1] + matrix1.m[2][3] * matrix2.m[3][1];

  result.m[2][2] =
      matrix1.m[2][0] * matrix2.m[0][2] + matrix1.m[2][1] * matrix2.m[1][2] +
      matrix1.m[2][2] * matrix2.m[2][2] + matrix1.m[2][3] * matrix2.m[3][2];

  result.m[2][3] =
      matrix1.m[2][0] * matrix2.m[0][3] + matrix1.m[2][1] * matrix2.m[1][3] +
      matrix1.m[2][2] * matrix2.m[2][3] + matrix1.m[2][3] * matrix2.m[3][3];

  result.m[3][0] =
      matrix1.m[3][0] * matrix2.m[0][0] + matrix1.m[3][1] * matrix2.m[1][0] +
      matrix1.m[3][2] * matrix2.m[2][0] + matrix1.m[3][3] * matrix2.m[3][0];

  result.m[3][1] =
      matrix1.m[3][0] * matrix2.m[0][1] + matrix1.m[3][1] * matrix2.m[1][1] +
      matrix1.m[3][2] * matrix2.m[2][1] + matrix1.m[3][3] * matrix2.m[3][1];

  result.m[3][2] =
      matrix1.m[3][0] * matrix2.m[0][2] + matrix1.m[3][1] * matrix2.m[1][2] +
      matrix1.m[3][2] * matrix2.m[2][2] + matrix1.m[3][3] * matrix2.m[3][2];

  result.m[3][3] =
      matrix1.m[3][0] * matrix2.m[0][3] + matrix1.m[3][1] * matrix2.m[1][3] +
      matrix1.m[3][2] * matrix2.m[2][3] + matrix1.m[3][3] * matrix2.m[3][3];

  return result;
}

// 座標変換
Vector3 Transform(const Vector3 &vector, const Matrix4x4 &matrix) {

  Vector3 result;

  result.x = vector.x * matrix.m[0][0] + vector.y * matrix.m[1][0] +
             vector.z * matrix.m[2][0] + 1.0f * matrix.m[3][0];

  result.y = vector.x * matrix.m[0][1] + vector.y * matrix.m[1][1] +
             vector.z * matrix.m[2][1] + 1.0f * matrix.m[3][1];

  result.z = vector.x * matrix.m[0][2] + vector.y * matrix.m[1][2] +
             vector.z * matrix.m[2][2] + 1.0f * matrix.m[3][2];

  float w = vector.x * matrix.m[0][3] + vector.y * matrix.m[1][3] +
            vector.z * matrix.m[2][3] + 1.0f * matrix.m[3][3];

  assert(w != 0.0f);
  result.x /= w;
  result.y /= w;
  result.z /= w;

  return result;
}

// アフィン行列作成
Matrix4x4 MakeAffineMatrix(const Vector3 &scale, const Vector3 &rotate,
                           const Vector3 &translate) {

  Matrix4x4 scaleMatrix;

  scaleMatrix.m[0][0] = scale.x;
  scaleMatrix.m[0][1] = 0.0f;
  scaleMatrix.m[0][2] = 0.0f;
  scaleMatrix.m[0][3] = 0.0f;
  scaleMatrix.m[1][0] = 0.0f;
  scaleMatrix.m[1][1] = scale.y;
  scaleMatrix.m[1][2] = 0.0f;
  scaleMatrix.m[1][3] = 0.0f;
  scaleMatrix.m[2][0] = 0.0f;
  scaleMatrix.m[2][1] = 0.0f;
  scaleMatrix.m[2][2] = scale.z;
  scaleMatrix.m[2][3] = 0.0f;
  scaleMatrix.m[3][0] = 0.0f;
  scaleMatrix.m[3][1] = 0.0f;
  scaleMatrix.m[3][2] = 0.0f;
  scaleMatrix.m[3][3] = 1.0f;

  Matrix4x4 rotateX;
  rotateX.m[0][0] = 1.0f;
  rotateX.m[0][1] = 0.0f;
  rotateX.m[0][2] = 0.0f;
  rotateX.m[0][3] = 0.0f;
  rotateX.m[1][0] = 0.0f;
  rotateX.m[1][1] = std::cos(rotate.x);
  rotateX.m[1][2] = std::sin(rotate.x);
  rotateX.m[1][3] = 0.0f;
  rotateX.m[2][0] = 0.0f;
  rotateX.m[2][1] = -std::sin(rotate.x);
  rotateX.m[2][2] = std::cos(rotate.x);
  rotateX.m[2][3] = 0.0f;
  rotateX.m[3][0] = 0.0f;
  rotateX.m[3][1] = 0.0f;
  rotateX.m[3][2] = 0.0f;
  rotateX.m[3][3] = 1.0f;

  Matrix4x4 rotateY;
  rotateY.m[0][0] = std::cos(rotate.y);
  rotateY.m[0][1] = 0.0f;
  rotateY.m[0][2] = -std::sin(rotate.y);
  rotateY.m[0][3] = 0.0f;
  rotateY.m[1][0] = 0.0f;
  rotateY.m[1][1] = 1.0f;
  rotateY.m[1][2] = 0.0f;
  rotateY.m[1][3] = 0.0f;
  rotateY.m[2][0] = std::sin(rotate.y);
  rotateY.m[2][1] = 0.0f;
  rotateY.m[2][2] = std::cos(rotate.y);
  rotateY.m[2][3] = 0.0f;
  rotateY.m[3][0] = 0.0f;
  rotateY.m[3][1] = 0.0f;
  rotateY.m[3][2] = 0.0f;
  rotateY.m[3][3] = 1.0f;

  Matrix4x4 rotateZ;
  rotateZ.m[0][0] = std::cos(rotate.z);
  rotateZ.m[0][1] = std::sin(rotate.z);
  rotateZ.m[0][2] = 0.0f;
  rotateZ.m[0][3] = 0.0f;
  rotateZ.m[1][0] = -std::sin(rotate.z);
  rotateZ.m[1][1] = std::cos(rotate.z);
  rotateZ.m[1][2] = 0.0f;
  rotateZ.m[1][3] = 0.0f;
  rotateZ.m[2][0] = 0.0f;
  rotateZ.m[2][1] = 0.0f;
  rotateZ.m[2][2] = 1.0f;
  rotateZ.m[2][3] = 0.0f;
  rotateZ.m[3][0] = 0.0f;
  rotateZ.m[3][1] = 0.0f;
  rotateZ.m[3][2] = 0.0f;
  rotateZ.m[3][3] = 1.0f;

  Matrix4x4 rotateMatrix = Multiply(rotateX, Multiply(rotateY, rotateZ));

  Matrix4x4 translateMatrix;
  translateMatrix.m[0][0] = 1.0f;
  translateMatrix.m[0][1] = 0.0f;
  translateMatrix.m[0][2] = 0.0f;
  translateMatrix.m[0][3] = 0.0f;
  translateMatrix.m[1][0] = 0.0f;
  translateMatrix.m[1][1] = 1.0f;
  translateMatrix.m[1][2] = 0.0f;
  translateMatrix.m[1][3] = 0.0f;
  translateMatrix.m[2][0] = 0.0f;
  translateMatrix.m[2][1] = 0.0f;
  translateMatrix.m[2][2] = 1.0f;
  translateMatrix.m[2][3] = 0.0f;
  translateMatrix.m[3][0] = translate.x;
  translateMatrix.m[3][1] = translate.y;
  translateMatrix.m[3][2] = translate.z;
  translateMatrix.m[3][3] = 1.0f;

  Matrix4x4 worldMatrix;
  worldMatrix = Multiply(scaleMatrix, Multiply(rotateMatrix, translateMatrix));

  return worldMatrix;
}

// 逆行列
Matrix4x4 Inverse(const Matrix4x4 &m) {
  Matrix4x4 inverseMatrix;

  // |A|
  float determinant = m.m[0][0] * m.m[1][1] * m.m[2][2] * m.m[3][3] +
                      m.m[0][0] * m.m[1][2] * m.m[2][3] * m.m[3][1] +
                      m.m[0][0] * m.m[1][3] * m.m[2][1] * m.m[3][2] -
                      m.m[0][0] * m.m[1][3] * m.m[2][2] * m.m[3][1] -
                      m.m[0][0] * m.m[1][2] * m.m[2][1] * m.m[3][3] -
                      m.m[0][0] * m.m[1][1] * m.m[2][3] * m.m[3][2] -
                      m.m[0][1] * m.m[1][0] * m.m[2][2] * m.m[3][3] -
                      m.m[0][2] * m.m[1][0] * m.m[2][3] * m.m[3][1] -
                      m.m[0][3] * m.m[1][0] * m.m[2][1] * m.m[3][2] +
                      m.m[0][3] * m.m[1][0] * m.m[2][2] * m.m[3][1] +
                      m.m[0][2] * m.m[1][0] * m.m[2][1] * m.m[3][3] +
                      m.m[0][1] * m.m[1][0] * m.m[2][3] * m.m[3][2] +
                      m.m[0][1] * m.m[1][2] * m.m[2][0] * m.m[3][3] +
                      m.m[0][2] * m.m[1][3] * m.m[2][0] * m.m[3][1] +
                      m.m[0][3] * m.m[1][1] * m.m[2][0] * m.m[3][2] -
                      m.m[0][3] * m.m[1][2] * m.m[2][0] * m.m[3][1] -
                      m.m[0][2] * m.m[1][1] * m.m[2][0] * m.m[3][3] -
                      m.m[0][1] * m.m[1][3] * m.m[2][0] * m.m[3][2] -
                      m.m[0][1] * m.m[1][2] * m.m[2][3] * m.m[3][0] -
                      m.m[0][2] * m.m[1][3] * m.m[2][1] * m.m[3][0] -
                      m.m[0][3] * m.m[1][1] * m.m[2][2] * m.m[3][0] +
                      m.m[0][3] * m.m[1][2] * m.m[2][1] * m.m[3][0] +
                      m.m[0][2] * m.m[1][1] * m.m[2][3] * m.m[3][0] +
                      m.m[0][1] * m.m[1][3] * m.m[2][2] * m.m[3][0];

  // 1/|A|
  float inverseDeterminant = 1 / determinant;

  inverseMatrix.m[0][0] =
      (m.m[1][1] * m.m[2][2] * m.m[3][3] - m.m[1][1] * m.m[2][3] * m.m[3][2] -
       m.m[1][2] * m.m[2][1] * m.m[3][3] + m.m[1][2] * m.m[2][3] * m.m[3][1] +
       m.m[1][3] * m.m[2][1] * m.m[3][2] - m.m[1][3] * m.m[2][2] * m.m[3][1]) *
      inverseDeterminant;

  inverseMatrix.m[0][1] =
      (-(m.m[0][1] * m.m[2][2] * m.m[3][3]) +
       m.m[0][1] * m.m[2][3] * m.m[3][2] + m.m[0][2] * m.m[2][1] * m.m[3][3] -
       m.m[0][2] * m.m[2][3] * m.m[3][1] - m.m[0][3] * m.m[2][1] * m.m[3][2] +
       m.m[0][3] * m.m[2][2] * m.m[3][1]) *
      inverseDeterminant;

  inverseMatrix.m[0][2] =
      (m.m[0][1] * m.m[1][2] * m.m[3][3] - m.m[0][1] * m.m[1][3] * m.m[3][2] -
       m.m[0][2] * m.m[1][1] * m.m[3][3] + m.m[0][2] * m.m[1][3] * m.m[3][1] +
       m.m[0][3] * m.m[1][1] * m.m[3][2] - m.m[0][3] * m.m[1][2] * m.m[3][1]) *
      inverseDeterminant;

  inverseMatrix.m[0][3] =
      (-(m.m[0][1] * m.m[1][2] * m.m[2][3]) +
       m.m[0][1] * m.m[1][3] * m.m[2][2] + m.m[0][2] * m.m[1][1] * m.m[2][3] -
       m.m[0][2] * m.m[1][3] * m.m[2][1] - m.m[0][3] * m.m[1][1] * m.m[2][2] +
       m.m[0][3] * m.m[1][2] * m.m[2][1]) *
      inverseDeterminant;
  inverseMatrix.m[1][0] =
      (-(m.m[1][0] * m.m[2][2] * m.m[3][3]) +
       m.m[1][0] * m.m[2][3] * m.m[3][2] + m.m[1][2] * m.m[2][0] * m.m[3][3] -
       m.m[1][2] * m.m[2][3] * m.m[3][0] - m.m[1][3] * m.m[2][0] * m.m[3][2] +
       m.m[1][3] * m.m[2][2] * m.m[3][0]) *
      inverseDeterminant;

  inverseMatrix.m[1][1] =
      (m.m[0][0] * m.m[2][2] * m.m[3][3] - m.m[0][0] * m.m[2][3] * m.m[3][2] -
       m.m[0][2] * m.m[2][0] * m.m[3][3] + m.m[0][2] * m.m[2][3] * m.m[3][0] +
       m.m[0][3] * m.m[2][0] * m.m[3][2] - m.m[0][3] * m.m[2][2] * m.m[3][0]) *
      inverseDeterminant;

  inverseMatrix.m[1][2] =
      (-(m.m[0][0] * m.m[1][2] * m.m[3][3]) +
       m.m[0][0] * m.m[1][3] * m.m[3][2] + m.m[0][2] * m.m[1][0] * m.m[3][3] -
       m.m[0][2] * m.m[1][3] * m.m[3][0] - m.m[0][3] * m.m[1][0] * m.m[3][2] +
       m.m[0][3] * m.m[1][2] * m.m[3][0]) *
      inverseDeterminant;

  inverseMatrix.m[1][3] =
      (m.m[0][0] * m.m[1][2] * m.m[2][3] - m.m[0][0] * m.m[1][3] * m.m[2][2] -
       m.m[0][2] * m.m[1][0] * m.m[2][3] + m.m[0][2] * m.m[1][3] * m.m[2][0] +
       m.m[0][3] * m.m[1][0] * m.m[2][2] - m.m[0][3] * m.m[1][2] * m.m[2][0]) *
      inverseDeterminant;

  inverseMatrix.m[2][0] =
      (m.m[1][0] * m.m[2][1] * m.m[3][3] - m.m[1][0] * m.m[2][3] * m.m[3][1] -
       m.m[1][1] * m.m[2][0] * m.m[3][3] + m.m[1][1] * m.m[2][3] * m.m[3][0] +
       m.m[1][3] * m.m[2][0] * m.m[3][1] - m.m[1][3] * m.m[2][1] * m.m[3][0]) *
      inverseDeterminant;

  inverseMatrix.m[2][1] =
      (-(m.m[0][0] * m.m[2][1] * m.m[3][3]) +
       m.m[0][0] * m.m[2][3] * m.m[3][1] + m.m[0][1] * m.m[2][0] * m.m[3][3] -
       m.m[0][1] * m.m[2][3] * m.m[3][0] - m.m[0][3] * m.m[2][0] * m.m[3][1] +
       m.m[0][3] * m.m[2][1] * m.m[3][0]) *
      inverseDeterminant;

  inverseMatrix.m[2][2] =
      (m.m[0][0] * m.m[1][1] * m.m[3][3] - m.m[0][0] * m.m[1][3] * m.m[3][1] -
       m.m[0][1] * m.m[1][0] * m.m[3][3] + m.m[0][1] * m.m[1][3] * m.m[3][0] +
       m.m[0][3] * m.m[1][0] * m.m[3][1] - m.m[0][3] * m.m[1][1] * m.m[3][0]) *
      inverseDeterminant;

  inverseMatrix.m[2][3] =
      (-(m.m[0][0] * m.m[1][1] * m.m[2][3]) +
       m.m[0][0] * m.m[1][3] * m.m[2][1] + m.m[0][1] * m.m[1][0] * m.m[2][3] -
       m.m[0][1] * m.m[1][3] * m.m[2][0] - m.m[0][3] * m.m[1][0] * m.m[2][1] +
       m.m[0][3] * m.m[1][1] * m.m[2][0]) *
      inverseDeterminant;

  inverseMatrix.m[3][0] =
      (-(m.m[1][0] * m.m[2][1] * m.m[3][2]) +
       m.m[1][0] * m.m[2][2] * m.m[3][1] + m.m[1][1] * m.m[2][0] * m.m[3][2] -
       m.m[1][1] * m.m[2][2] * m.m[3][0] - m.m[1][2] * m.m[2][0] * m.m[3][1] +
       m.m[1][2] * m.m[2][1] * m.m[3][0]) *
      inverseDeterminant;

  inverseMatrix.m[3][1] =
      (m.m[0][0] * m.m[2][1] * m.m[3][2] - m.m[0][0] * m.m[2][2] * m.m[3][1] -
       m.m[0][1] * m.m[2][0] * m.m[3][2] + m.m[0][1] * m.m[2][2] * m.m[3][0] +
       m.m[0][2] * m.m[2][0] * m.m[3][1] - m.m[0][2] * m.m[2][1] * m.m[3][0]) *
      inverseDeterminant;

  inverseMatrix.m[3][2] =
      (-(m.m[0][0] * m.m[1][1] * m.m[3][2]) +
       m.m[0][0] * m.m[1][2] * m.m[3][1] + m.m[0][1] * m.m[1][0] * m.m[3][2] -
       m.m[0][1] * m.m[1][2] * m.m[3][0] - m.m[0][2] * m.m[1][0] * m.m[3][1] +
       m.m[0][2] * m.m[1][1] * m.m[3][0]) *
      inverseDeterminant;

  inverseMatrix.m[3][3] =
      (m.m[0][0] * m.m[1][1] * m.m[2][2] - m.m[0][0] * m.m[1][2] * m.m[2][1] -
       m.m[0][1] * m.m[1][0] * m.m[2][2] + m.m[0][1] * m.m[1][2] * m.m[2][0] +
       m.m[0][2] * m.m[1][0] * m.m[2][1] - m.m[0][2] * m.m[1][1] * m.m[2][0]) *
      inverseDeterminant;

  return inverseMatrix;
}

// 透視投影行列
Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio,
                                   float nearClip, float farClip) {
  Matrix4x4 matrix;

  matrix.m[0][0] =
      1.0f / aspectRatio * (std::cos(fovY / 2.0f) / std::sin(fovY / 2.0f));
  matrix.m[0][1] = 0.0f;
  matrix.m[0][2] = 0.0f;
  matrix.m[0][3] = 0.0f;
  matrix.m[1][0] = 0.0f;
  matrix.m[1][1] = std::cos(fovY / 2.0f) / std::sin(fovY / 2.0f);
  matrix.m[1][2] = 0.0f;
  matrix.m[1][3] = 0.0f;
  matrix.m[2][0] = 0.0f;
  matrix.m[2][1] = 0.0f;
  matrix.m[2][2] = farClip / (farClip - nearClip);
  matrix.m[2][3] = 1.0f;
  matrix.m[3][0] = 0.0f;
  matrix.m[3][1] = 0.0f;
  matrix.m[3][2] = (-nearClip * farClip) / (farClip - nearClip);
  matrix.m[3][3] = 0.0f;

  return matrix;
}

// ビューポート変換行列
Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height,
                             float minDepth, float maxDepth) {

  Matrix4x4 matrix;

  matrix.m[0][0] = width / 2.0f;
  matrix.m[0][1] = 0.0f;
  matrix.m[0][2] = 0.0f;
  matrix.m[0][3] = 0.0f;
  matrix.m[1][0] = 0.0f;
  matrix.m[1][1] = -(height / 2.0f);
  matrix.m[1][2] = 0.0f;
  matrix.m[1][3] = 0.0f;
  matrix.m[2][0] = 0.0f;
  matrix.m[2][1] = 0.0f;
  matrix.m[2][2] = maxDepth - minDepth;
  matrix.m[2][3] = 0.0f;
  matrix.m[3][0] = left + (width / 2.0f);
  matrix.m[3][1] = top + (height / 2.0f);
  matrix.m[3][2] = minDepth;
  matrix.m[3][3] = 1.0f;

  return matrix;
}

/// <summary>
/// ベクトルの足し算
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
Vector3 Add(const Vector3 &v1, const Vector3 &v2) {

  Vector3 result;

  result.x = v1.x + v2.x;
  result.y = v1.y + v2.y;
  result.z = v1.z + v2.z;

  return result;
}

/// <summary>
/// 内積を求める
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
float Dot(const Vector3 &v1, const Vector3 &v2) {
  float result;

  result = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;

  return result;
}

/// <summary>
/// 長さを求める
/// </summary>
/// <param name="vector"></param>
/// <returns></returns>
float Length(const Vector3 &vector) {
  float result;

  result = sqrtf(Dot(vector, vector));

  return result;
}

/// <summary>
/// 外積
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
Vector3 Cross(const Vector3 &v1, const Vector3 &v2) {

  Vector3 result;

  result.x = v1.y * v2.z - v1.z * v2.y;
  result.y = v1.z * v2.x - v1.x * v2.z;
  result.z = v1.x * v2.y - v1.y * v2.x;

  return result;
}

/// <summary>
/// Vector3とfloatの積を求める
/// </summary>
/// <param name="f"></param>
/// <param name="vector"></param>
/// <returns></returns>
Vector3 Multiply(const float &f, const Vector3 vector) {

  Vector3 result;

  result = {vector.x * f, vector.y * f, vector.z * f};

  return result;
}

Vector3 Perpendicular(const Vector3 &vector) {

  if (vector.x != 0.0f || vector.y != 0.0f) {
    return {-vector.y, vector.x, 0.0f};
  }

  return {0.0f, -vector.z, vector.y};
}

Vector3 Normalize(const Vector3 &vector) {

  float length = Length(vector);

  Vector3 result;

  result = {vector.x / length, vector.y / length, vector.z / length};

  return result;
}

/// <summary>
/// Gridを描画する
/// </summary>
/// <param name="viewProjectionMatrix">ビュー射影行列</param>
/// <param name="viewPortMatrix">ビューポート行列</param>
void DrawGlid(const Matrix4x4 &viewProjectionMatrix,
              const Matrix4x4 &viewPortMatrix) {

  const float kGridHalfWidth = 2.0f; // Gridの半分の幅
  const uint32_t kSubdivision = 10;  // 分割数
  const float kGridEvery =
      (kGridHalfWidth * 2.0f) / float(kSubdivision); // 一つ分の長さ

  // 奥から手前への線を引いていく
  for (int xIndex = 0; xIndex <= kSubdivision; ++xIndex) {

    // 始点と終点を求める
    Vector3 worldGridStart{xIndex * kGridEvery - kGridHalfWidth, 0.0f,
                           -kGridHalfWidth};
    Vector3 worldGridEnd{xIndex * kGridEvery - kGridHalfWidth, 0.0f,
                         kGridHalfWidth};

    // 座標変換
    Matrix4x4 worldMatrix = MakeAffineMatrix(
        {1.0f, 1.0f, 1.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f});
    Matrix4x4 worldViewProjectionMatrix =
        Multiply(worldMatrix, viewProjectionMatrix);

    Vector3 screenGridStart = Transform(
        Transform(worldGridStart, worldViewProjectionMatrix), viewPortMatrix);
    Vector3 screenGridEnd = Transform(
        Transform(worldGridEnd, worldViewProjectionMatrix), viewPortMatrix);

    // 描画
    Novice::DrawLine(int(screenGridStart.x), int(screenGridStart.y),
                     int(screenGridEnd.x), int(screenGridEnd.y), 0xaaaaaaff);

    // 原点は黒で描画
    if (xIndex == (kSubdivision + 1) / 2) {
      Novice::DrawLine(int(screenGridStart.x), int(screenGridStart.y),
                       int(screenGridEnd.x), int(screenGridEnd.y), 0x000000ff);
    }
  }

  // 左から右
  for (int zIndex = 0; zIndex <= kSubdivision; ++zIndex) {

    // 始点と終点を求める
    Vector3 worldGridStart{-kGridHalfWidth, 0.0f,
                           zIndex * kGridEvery - kGridHalfWidth};
    Vector3 worldGridEnd{kGridHalfWidth, 0.0f,
                         zIndex * kGridEvery - kGridHalfWidth};

    // 座標変換
    Matrix4x4 worldMatrix = MakeAffineMatrix(
        {1.0f, 1.0f, 1.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f});
    Matrix4x4 worldViewProjectionMatrix =
        Multiply(worldMatrix, viewProjectionMatrix);

    Vector3 screenGridStart = Transform(
        Transform(worldGridStart, worldViewProjectionMatrix), viewPortMatrix);
    Vector3 screenGridEnd = Transform(
        Transform(worldGridEnd, worldViewProjectionMatrix), viewPortMatrix);

    // 描画
    Novice::DrawLine(int(screenGridStart.x), int(screenGridStart.y),
                     int(screenGridEnd.x), int(screenGridEnd.y), 0xaaaaaaff);
    // 原点は黒で描画
    if (zIndex == (kSubdivision + 1) / 2) {
      Novice::DrawLine(int(screenGridStart.x), int(screenGridStart.y),
                       int(screenGridEnd.x), int(screenGridEnd.y), 0x000000ff);
    }
  }
}

/// <summary>
/// 球の描画
/// </summary>
/// <param name="sphere">球</param>
/// <param name="viewProjectionMatrix">ビュー射影行列</param>
/// <param name="viewportMatrix">ビューポート行列</param>
/// <param name="color">色</param>
void DrawSphere(const Sphere &sphere, const Matrix4x4 &viewProjectionMatrix,
                const Matrix4x4 viewportMatrix, uint32_t color) {

  const uint32_t kSubdivision = 15; // 分割数
  const float kLonEvery =
      static_cast<float>(M_PI) * 2.0f /
      static_cast<float>(kSubdivision); // 経度分割一つ分の角度
  const float kLatEvery =
      static_cast<float>(M_PI) /
      static_cast<float>(kSubdivision); // 緯度分割一つ分の角度

  // 緯度の方向に分割 -π/2 ~ π/2
  for (uint32_t latIndex = 0; latIndex <= kSubdivision; ++latIndex) {

    float lat =
        static_cast<float>(M_PI) / 2.0f + kLatEvery * latIndex; // 現在の緯度

    // 経度の方向に分割 0~2π
    for (uint32_t lonIndex = 0; lonIndex <= kSubdivision; ++lonIndex) {

      float lon = lonIndex * kLonEvery; // 現在の経度

      // world座標系でのabcを求める
      Vector3 a, b, c;

      a = {sphere.center.x + sphere.radius * std::cosf(lat) * std::cosf(lon),
           sphere.center.y + sphere.radius * std::sinf(lat),
           sphere.center.z + sphere.radius * std::cosf(lat) * std::sinf(lon)};

      b = {sphere.center.x +
               sphere.radius * std::cosf(lat + kLatEvery) * std::cosf(lon),
           sphere.center.y + sphere.radius * std::sinf(lat + kLatEvery),
           sphere.center.z +
               sphere.radius * std::cosf(lat + kLatEvery) * std::sinf(lon)};

      c = {sphere.center.x +
               sphere.radius * std::cosf(lat) * std::cosf(lon + kLonEvery),
           sphere.center.y + sphere.radius * std::sinf(lat),
           sphere.center.z +
               sphere.radius * std::cosf(lat) * std::sinf(lon + kLonEvery)};

      // a,b,cをScreen座標系まで変換

      Matrix4x4 worldMatrix = MakeAffineMatrix(
          {1.0f, 1.0f, 1.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f});
      Matrix4x4 worldViewProjectionMatrix =
          Multiply(worldMatrix, viewProjectionMatrix);

      Vector3 screenA, screenB, screenC;

      screenA =
          Transform(Transform(a, worldViewProjectionMatrix), viewportMatrix);
      screenB =
          Transform(Transform(b, worldViewProjectionMatrix), viewportMatrix);
      screenC =
          Transform(Transform(c, worldViewProjectionMatrix), viewportMatrix);

      // AB
      Novice::DrawLine(static_cast<int>(screenA.x), static_cast<int>(screenA.y),
                       static_cast<int>(screenB.x), static_cast<int>(screenB.y),
                       color);

      // BC
      Novice::DrawLine(static_cast<int>(screenA.x), static_cast<int>(screenA.y),
                       static_cast<int>(screenC.x), static_cast<int>(screenC.y),
                       color);
    }
  }
}

/// <summary>
/// 平面の描画
/// </summary>
/// <param name="plane">平面</param>
/// <param name="viewProjectionMatrix">ビュー射影行列</param>
/// <param name="viewMatrix">ビューポート行列</param>
/// <param name="color">色</param>
void DrawPlane(const Plane &plane, const Matrix4x4 &viewProjectionMatrix,
               const Matrix4x4 &viewportMatrix, uint32_t color) {
  // 中心点を決める
  Vector3 center = Multiply(plane.distance, plane.normal);

  Vector3 perpendiculars[4];

  // 法線と垂直なベクトルを一つ求める
  perpendiculars[0] = Normalize(Perpendicular(plane.normal));

  // ↑の逆ベクトルを求める
  perpendiculars[1] = {
      -perpendiculars[0].x,
      -perpendiculars[0].y,
      -perpendiculars[0].z,
  };

  //[0]の法線とクロス積を求める
  perpendiculars[2] = Cross(plane.normal, perpendiculars[0]);

  // ↑の逆ベクトルを求める
  perpendiculars[3] = {
      -perpendiculars[2].x,
      -perpendiculars[2].y,
      -perpendiculars[2].z,
  };

  // 四頂点を求める
  Vector3 points[4];
  for (int32_t index = 0; index < 4; ++index) {
    Vector3 extend = Multiply(2.0f, perpendiculars[index]);
    Vector3 point = Add(center, extend);
    points[index] =
        Transform(Transform(point, viewProjectionMatrix), viewportMatrix);
  }

  // 矩形を描画

  Novice::DrawLine(static_cast<int>(points[0].x), static_cast<int>(points[0].y),
                   static_cast<int>(points[2].x), static_cast<int>(points[2].y),
                   color);
  Novice::DrawLine(static_cast<int>(points[1].x), static_cast<int>(points[1].y),
                   static_cast<int>(points[3].x), static_cast<int>(points[3].y),
                   color);
  Novice::DrawLine(static_cast<int>(points[2].x), static_cast<int>(points[2].y),
                   static_cast<int>(points[1].x), static_cast<int>(points[1].y),
                   color);
  Novice::DrawLine(static_cast<int>(points[3].x), static_cast<int>(points[3].y),
                   static_cast<int>(points[0].x), static_cast<int>(points[0].y),
                   color);
}

/// <summary>
/// ベクトルの引き算
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
Vector3 Subtract(const Vector3 &v1, const Vector3 &v2) {

  Vector3 result;

  result.x = v1.x - v2.x;
  result.y = v1.y - v2.y;
  result.z = v1.z - v2.z;

  return result;
}

/// <summary>
/// 球の当たり判定
/// </summary>
/// <param name="s1"></param>
/// <param name="s2"></param>
/// <returns></returns>
bool IsCollision(const Sphere &s1, const Sphere &s2) {

  float distance = Length(Subtract(s2.center, s1.center));

  if (distance <= s1.radius + s2.radius) {
    return true;
  }
  return false;
}

/// <summary>
/// 球と平面の当たり判定
/// </summary>
/// <param name="sphere"></param>
/// <param name="plane"></param>
/// <returns></returns>
bool IsCollision(const Sphere &sphere, const Plane &plane) {

  // 平面と球の中心点の距離を求める
  float k = Dot(plane.normal, sphere.center) - plane.distance;

  // 絶対値にする
  if (k < 0) {
    k *= -1.0f;
  }

  if (sphere.radius >= k) {
    return true;
  } else {
    return false;
  }
}

/// <summary>
/// 線分と平面の当たり判定
/// </summary>
/// <param name="segment"></param>
/// <param name="plane"></param>
/// <returns></returns>
bool IsCollision(const Segment &segment, const Plane &plane) {

  float dot = Dot(segment.diff, plane.normal);

  // 平行なので衝突しない
  if (dot == 0.0f) {
    return false;
  }

  float t = (plane.distance - Dot(segment.origin, plane.normal)) / dot;

  // 線分
  if (t > 0.0f && t <= 1.0f) {
    return true;
  } else {
    return false;
  }
}

/// <summary>
/// 線分の描画
/// </summary>
/// <param name="segment"></param>
/// <param name="viewProjectionMatrix"></param>
/// <param name="viewportMatrix"></param>
/// <param name="color"></param>
void DrawSegment(const Segment &segment, const Matrix4x4 &viewProjectionMatrix,
                 const Matrix4x4 &viewportMatrix, uint32_t color) {

  Vector3 start = Transform(Transform(segment.origin, viewProjectionMatrix),
                            viewportMatrix);

  Vector3 end = Transform(
      Transform(Add(segment.origin, segment.diff), viewProjectionMatrix),
      viewportMatrix);

  Novice::DrawLine(static_cast<int>(start.x), static_cast<int>(start.y),
                   static_cast<int>(end.x), static_cast<int>(end.y), color);
}

/// <summary>
/// 三角形と線分の当たり判定
/// </summary>
/// <param name="triangle"></param>
/// <param name="segment"></param>
/// <returns></returns>
bool IsCollision(const Triangle &triangle, const Segment &segment) {

  Vector3 v01 = Subtract(triangle.vertices[1], triangle.vertices[0]);
  Vector3 v12 = Subtract(triangle.vertices[2], triangle.vertices[1]);
  Vector3 v20 = Subtract(triangle.vertices[0], triangle.vertices[2]);

  Vector3 normal = Normalize(Cross(v01, v12));

  float dot = Dot(segment.diff, normal);

  // 平行なので衝突しない
  if (dot == 0.0f) {
    return false;
  }

  float d = Dot(triangle.vertices[0], normal);

  float t = (d - Dot(segment.origin, normal)) / dot;

  // 線分の範囲外は衝突しない
  if (t < 0.0f || t > 1.0f) {
    return false;
  }

  Vector3 p = Add(segment.origin, Multiply(t, segment.diff));

  Vector3 v0p = Subtract(p, triangle.vertices[0]);
  Vector3 v1p = Subtract(p, triangle.vertices[1]);
  Vector3 v2p = Subtract(p, triangle.vertices[2]);

  Vector3 cross01 = Cross(v01, v1p);
  Vector3 cross12 = Cross(v12, v2p);
  Vector3 cross20 = Cross(v20, v0p);

  if (Dot(cross01, normal) >= 0.0f && Dot(cross12, normal) >= 0.0f &&
      Dot(cross20, normal) >= 0.0f) {
    return true;
  } else {
    return false;
  }
}

/// <summary>
/// 三角形の描画
/// </summary>
/// <param name="triangle"></param>
/// <param name="viewProjectionMatrix"></param>
/// <param name="viewportMatrix"></param>
/// <param name="color"></param>
void DrawTriangle(const Triangle &triangle,
                  const Matrix4x4 viewProjectionMatrix,
                  const Matrix4x4 viewportMatrix, uint32_t color) {

  Vector3 v0 = Transform(Transform(triangle.vertices[0], viewProjectionMatrix),
                         viewportMatrix);
  Vector3 v1 = Transform(Transform(triangle.vertices[1], viewProjectionMatrix),
                         viewportMatrix);
  Vector3 v2 = Transform(Transform(triangle.vertices[2], viewProjectionMatrix),
                         viewportMatrix);

  Novice::DrawTriangle(static_cast<int>(v0.x), static_cast<int>(v0.y),
                       static_cast<int>(v1.x), static_cast<int>(v1.y),
                       static_cast<int>(v2.x), static_cast<int>(v2.y), color,
                       kFillModeWireFrame);
}

/// <summary>
/// マウスでカメラ操作
/// </summary>
/// <param name="cameraTranslate"></param>
/// <param name="cameraRotate"></param>
void UpdateCameraByMouse(Vector3 &cameraTranslate, Vector3 &cameraRotate) {
  // マウスでカメラ移動// マウス座標を取得
  // マウス座標取得
  int mouseX, mouseY;
  Novice::GetMousePosition(&mouseX, &mouseY);

  // 状態保持
  static int prevMouseX = mouseX;
  static int prevMouseY = mouseY;
  static bool wasRotating = false;
  static bool wasPanning = false;

  // 左ボタン回転
  bool isRotating = Novice::IsPressMouse(0); // 左ボタン
  bool allowRotate =
      isRotating && !ImGui::IsAnyItemActive() && !ImGui::IsAnyItemHovered();

  // 右ボタン平行移動
  bool isPanning = Novice::IsPressMouse(1); // 右ボタン
  bool allowPan =
      isPanning && !ImGui::IsAnyItemActive() && !ImGui::IsAnyItemHovered();

  // 視点回転（左ドラッグ）
  if (allowRotate) {
    if (!wasRotating) {
      prevMouseX = mouseX;
      prevMouseY = mouseY;
    }

    float dx = static_cast<float>(mouseX - prevMouseX);
    float dy = static_cast<float>(mouseY - prevMouseY);

    cameraRotate.y -= dx * 0.001f;
    cameraRotate.x -= dy * 0.001f;

    const float piOver2 = 3.141592f / 2.0f;
    cameraRotate.x = std::clamp(cameraRotate.x, -piOver2, piOver2);

    prevMouseX = mouseX;
    prevMouseY = mouseY;
  }

  // カメラ平行移動（右ドラッグ）
  if (allowPan) {
    if (!wasPanning) {
      prevMouseX = mouseX;
      prevMouseY = mouseY;
    }

    float dx = static_cast<float>(mouseX - prevMouseX);
    float dy = static_cast<float>(mouseY - prevMouseY);

    // カメラの向きに応じてX方向とY方向に動くように補正
    float speed = 0.01f;
    cameraTranslate.x -= dx * speed;
    cameraTranslate.y += dy * speed;

    prevMouseX = mouseX;
    prevMouseY = mouseY;
  }

  wasRotating = isRotating;
  wasPanning = isPanning;

  // ホイールズーム
  ImGuiIO &io = ImGui::GetIO();
  cameraTranslate.z += io.MouseWheel * 0.5f;
}

/// <summary>
/// AABBの当たり判定
/// </summary>
/// <param name="a"></param>
/// <param name="b"></param>
/// <returns></returns>
bool IsCollision(const AABB &a, const AABB &b) {

  if ((a.min.x <= b.max.x && a.max.x >= b.min.x) && // x軸
      (a.min.y <= b.max.y && a.max.y >= b.min.y) && // y軸
      (a.min.z <= b.max.z && a.max.z >= b.min.z)) { // z軸

    return true;

  } else {
    return false;
  }
}

/// <summary>
/// AABB描画
/// </summary>
/// <param name="aabb"></param>
/// <param name="viewProjectionMatrix"></param>
/// <param name="viewportMatrix"></param>
/// <param name="color"></param>
void DrawAABB(const AABB &aabb, const Matrix4x4 &viewProjectionMatrix,
              const Matrix4x4 &viewportMatrix, uint32_t color) {

  Vector3 vertex[8];

  vertex[0] = {aabb.min.x, aabb.min.y, aabb.min.z};
  vertex[1] = {aabb.max.x, aabb.min.y, aabb.min.z};
  vertex[2] = {aabb.min.x, aabb.max.y, aabb.min.z};
  vertex[3] = {aabb.max.x, aabb.max.y, aabb.min.z};
  vertex[4] = {aabb.min.x, aabb.min.y, aabb.max.z};
  vertex[5] = {aabb.max.x, aabb.min.y, aabb.max.z};
  vertex[6] = {aabb.min.x, aabb.max.y, aabb.max.z};
  vertex[7] = {aabb.max.x, aabb.max.y, aabb.max.z};

  Vector3 screenVertex[8];

  for (int i = 0; i < 8; i++) {
    screenVertex[i] =
        Transform(Transform(vertex[i], viewProjectionMatrix), viewportMatrix);
  }

  Novice::DrawLine(static_cast<int>(screenVertex[0].x),
                   static_cast<int>(screenVertex[0].y),
                   static_cast<int>(screenVertex[1].x),
                   static_cast<int>(screenVertex[1].y), color);

  Novice::DrawLine(static_cast<int>(screenVertex[0].x),
                   static_cast<int>(screenVertex[0].y),
                   static_cast<int>(screenVertex[2].x),
                   static_cast<int>(screenVertex[2].y), color);

  Novice::DrawLine(static_cast<int>(screenVertex[0].x),
                   static_cast<int>(screenVertex[0].y),
                   static_cast<int>(screenVertex[4].x),
                   static_cast<int>(screenVertex[4].y), color);

  Novice::DrawLine(static_cast<int>(screenVertex[1].x),
                   static_cast<int>(screenVertex[1].y),
                   static_cast<int>(screenVertex[3].x),
                   static_cast<int>(screenVertex[3].y), color);

  Novice::DrawLine(static_cast<int>(screenVertex[1].x),
                   static_cast<int>(screenVertex[1].y),
                   static_cast<int>(screenVertex[5].x),
                   static_cast<int>(screenVertex[5].y), color);

  Novice::DrawLine(static_cast<int>(screenVertex[2].x),
                   static_cast<int>(screenVertex[2].y),
                   static_cast<int>(screenVertex[3].x),
                   static_cast<int>(screenVertex[3].y), color);

  Novice::DrawLine(static_cast<int>(screenVertex[2].x),
                   static_cast<int>(screenVertex[2].y),
                   static_cast<int>(screenVertex[6].x),
                   static_cast<int>(screenVertex[6].y), color);

  Novice::DrawLine(static_cast<int>(screenVertex[3].x),
                   static_cast<int>(screenVertex[3].y),
                   static_cast<int>(screenVertex[7].x),
                   static_cast<int>(screenVertex[7].y), color);

  Novice::DrawLine(static_cast<int>(screenVertex[4].x),
                   static_cast<int>(screenVertex[4].y),
                   static_cast<int>(screenVertex[5].x),
                   static_cast<int>(screenVertex[5].y), color);

  Novice::DrawLine(static_cast<int>(screenVertex[4].x),
                   static_cast<int>(screenVertex[4].y),
                   static_cast<int>(screenVertex[6].x),
                   static_cast<int>(screenVertex[6].y), color);

  Novice::DrawLine(static_cast<int>(screenVertex[5].x),
                   static_cast<int>(screenVertex[5].y),
                   static_cast<int>(screenVertex[7].x),
                   static_cast<int>(screenVertex[7].y), color);

  Novice::DrawLine(static_cast<int>(screenVertex[6].x),
                   static_cast<int>(screenVertex[6].y),
                   static_cast<int>(screenVertex[7].x),
                   static_cast<int>(screenVertex[7].y), color);
}

const int kWindowWidth = 1280;
const int kWindowHeight = 720;

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

  // ライブラリの初期化
  Novice::Initialize(kWindowTitle, kWindowWidth, kWindowHeight);

  // キー入力結果を受け取る箱
  char keys[256] = {0};
  char preKeys[256] = {0};

  Vector3 cameraTranslate{0.0f, 1.9f, -6.49f};
  Vector3 cameraRotate{0.26f, 0.0f, 0.0f};

  // 色
  uint32_t color = WHITE;

  AABB aabb1{
      .min{-0.5f, -0.5f, -0.5f},
      .max{0.0f, 0.0f, 0.0f},
  };

  AABB aabb2{
      .min{0.2f, 0.2f, 0.2f},
      .max{1.0f, 1.0f, 1.0f},
  };

  // ウィンドウの×ボタンが押されるまでループ
  while (Novice::ProcessMessage() == 0) {
    // フレームの開始
    Novice::BeginFrame();

    // キー入力を受け取る
    memcpy(preKeys, keys, 256);
    Novice::GetHitKeyStateAll(keys);

    ///
    /// ↓更新処理ここから
    ///

    Matrix4x4 cameraMatrix =
        MakeAffineMatrix({1.0f, 1.0f, 1.0f}, cameraRotate, cameraTranslate);
    Matrix4x4 viewMatrix = Inverse(cameraMatrix);
    Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(
        0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
    Matrix4x4 viewProjectionMatrix = Multiply(viewMatrix, projectionMatrix);
    Matrix4x4 viewportMatrix = MakeViewportMatrix(
        0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);

    // minとmaxが入れ替わらないように
    aabb1.min.x = (std::min)(aabb1.min.x, aabb1.max.x);
    aabb1.max.x = (std::max)(aabb1.min.x, aabb1.max.x);
    aabb1.min.y = (std::min)(aabb1.min.y, aabb1.max.y);
    aabb1.max.y = (std::max)(aabb1.min.y, aabb1.max.y);
    aabb1.min.z = (std::min)(aabb1.min.z, aabb1.max.z);
    aabb1.max.z = (std::max)(aabb1.min.z, aabb1.max.z);

    // minとmaxが入れ替わらないように
    aabb2.min.x = (std::min)(aabb2.min.x, aabb2.max.x);
    aabb2.max.x = (std::max)(aabb2.min.x, aabb2.max.x);
    aabb2.min.y = (std::min)(aabb2.min.y, aabb2.max.y);
    aabb2.max.y = (std::max)(aabb2.min.y, aabb2.max.y);
    aabb2.min.z = (std::min)(aabb2.min.z, aabb2.max.z);
    aabb2.max.z = (std::max)(aabb2.min.z, aabb2.max.z);

    // 当たり判定を取る
    if (IsCollision(aabb1, aabb2)) {
      color = RED;
    } else {
      color = WHITE;
    }

    UpdateCameraByMouse(cameraTranslate, cameraRotate);

    ///
    /// ↑更新処理ここまで
    ///

    ///
    /// ↓描画処理ここから
    ///

    // グリッド
    DrawGlid(viewProjectionMatrix, viewportMatrix);

    // AABB
    DrawAABB(aabb1, viewProjectionMatrix, viewportMatrix, color);

    DrawAABB(aabb2, viewProjectionMatrix, viewportMatrix, WHITE);

    // UI
    ImGui::Begin("Window");
    ImGui::DragFloat3("CameraTranslate", &cameraTranslate.x, 0.01f);
    ImGui::DragFloat3("CameraRotate", &cameraRotate.x, 0.001f);
    ImGui::DragFloat3("aabb1.min", &aabb1.min.x, 0.01f);
    ImGui::DragFloat3("aabb1.max", &aabb1.max.x, 0.01f);
    ImGui::DragFloat3("aabb2.min", &aabb2.min.x, 0.01f);
    ImGui::DragFloat3("aabb2.max", &aabb2.max.x, 0.01f);
    ImGui::End();

    ///
    /// ↑描画処理ここまで
    ///

    // フレームの終了
    Novice::EndFrame();

    // ESCキーが押されたらループを抜ける
    if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
      break;
    }
  }

  // ライブラリの終了
  Novice::Finalize();
  return 0;
}
