#include <Novice.h>
#include <cassert>
#include <cmath>
#define _USE_MATH_DEFINES
#include <imgui.h>
#include <math.h>

const char kWindowTitle[] = "LC2B_09_コジマユウヤ";

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

struct Segment {
  Vector3 origin;
  Vector3 diff;
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
      matrix1.m[0][2] * matrix2.m[2][3] + matrix1.m[3][0] * matrix2.m[3][3];

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
/// 内積求める
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
/// 正射影ベクトルを求める
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
Vector3 Project(const Vector3 &v1, const Vector3 &v2) {

  Vector3 result;

  result.y = (Dot(v1, v2) /
              Dot(v2, v2) /*std::sqrt(Dot(v2, v2)) * std::sqrt(Dot(v2, v2))*/) *
             v2.y;
  result.z = (Dot(v1, v2) /
              Dot(v2, v2) /*std::sqrt(Dot(v2, v2)) * std::sqrt(Dot(v2, v2))*/) *
             v2.z;
  result.x = (Dot(v1, v2) /
              Dot(v2, v2) /*std::sqrt(Dot(v2, v2)) * std::sqrt(Dot(v2, v2))*/) *
             v2.x;

  return result;
}

/// <summary>
/// 最近接点を求める
/// </summary>
/// <param name="point"></param>
/// <param name="segment"></param>
/// <returns></returns>
Vector3 ClosestPoint(const Vector3 &point, const Segment &segment) {

  Vector3 result;

  result = Add(segment.origin, point);

  return result;
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

  Segment segment{
      {-2.0f, -1.0f, 0.0f},
      {3.0f, 2.0f, 2.0f},
  };

  Vector3 point{-1.5f, 0.6f, 0.6f};

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

    // 正射影ベクトル
    Vector3 project = Project(Subtract(point, segment.origin), segment.diff);

    // 最近接点
    Vector3 closestPoint = ClosestPoint(project, segment);

    // 描画用の球
    Sphere pointSphere{point, 0.01f};
    Sphere closestPointSphere{closestPoint, 0.01f};

    // 線分
    Vector3 start = Transform(Transform(segment.origin, viewProjectionMatrix),
                              viewportMatrix);
    Vector3 end = Transform(
        Transform(Add(segment.origin, segment.diff), viewProjectionMatrix),
        viewportMatrix);

    ///
    /// ↑更新処理ここまで
    ///

    ///
    /// ↓描画処理ここから
    ///

    // グリッド描画
    DrawGlid(viewProjectionMatrix, viewportMatrix);

    // 線分描画
    Novice::DrawLine(static_cast<int>(start.x), static_cast<int>(start.y),
                     static_cast<int>(end.x), static_cast<int>(end.y), WHITE);

    // 元の点
    DrawSphere(pointSphere, viewProjectionMatrix, viewportMatrix, RED);

    // 最近接点
    DrawSphere(closestPointSphere, viewProjectionMatrix, viewportMatrix, BLACK);

    // UI
    ImGui::Begin("Window");

    ImGui::InputFloat3("point", &point.x, "%.3f", ImGuiInputTextFlags_ReadOnly);

    ImGui::InputFloat3("segment origin", &segment.origin.x, "%.3f",
                       ImGuiInputTextFlags_ReadOnly);

    ImGui::InputFloat3("segment diff", &segment.diff.x, "%.3f",
                       ImGuiInputTextFlags_ReadOnly);

    ImGui::InputFloat3("Project", &project.x, "%.3f",
                       ImGuiInputTextFlags_ReadOnly);

    ImGui::DragFloat3("camera translate", &cameraTranslate.x, 0.01f);
    ImGui::DragFloat3("camera rotate", &cameraRotate.x, 0.01f);

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
