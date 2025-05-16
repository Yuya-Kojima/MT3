#include <Novice.h>
#include <cassert>
#include <cmath>
#define _USE_MATH_DEFINES
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

  // 球の色
  uint32_t color = 0x000000ff;

  Sphere sphere[2] = {

      {
          // 一つ目(動かす方)
          {0.0f, 0.0f, 0.0f},
          0.5f,
      },

      {
          // 二つ目
          {1.0f, 0.0f, 0.0f},
          {0.5f},
      },
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

    // 当たり判定を取る
    // trueなら赤色にする
    // falseなら白色
    if (IsCollision(sphere[0], sphere[1])) {
      color = 0xff0000ff;
    } else {
      color = 0xffffffff;
    }

    ///
    /// ↑更新処理ここまで
    ///

    ///
    /// ↓描画処理ここから
    ///

    // グリッド
    DrawGlid(viewProjectionMatrix, viewportMatrix);

    // 動かす方の球
    DrawSphere(sphere[0], viewProjectionMatrix, viewportMatrix, color);

    // 固定されてる球
    DrawSphere(sphere[1], viewProjectionMatrix, viewportMatrix, 0xffffffff);

    // UI
    ImGui::Begin("Window");
    ImGui::DragFloat3("CameraTranslate", &cameraTranslate.x, 0.01f);
    ImGui::DragFloat3("CameraRotate", &cameraRotate.x, 0.01f);
    ImGui::DragFloat3("SphereTranslate", &sphere[0].center.x, 0.01f);
    ImGui::DragFloat("SphereRadius", &sphere[0].radius, 0.01f);
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
