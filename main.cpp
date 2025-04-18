#include <Novice.h>
#include <assert.h>
#include <cmath>

const char kWindowTitle[] = "LC2B_09_コジマユウヤ";

struct Vector3 {
  float x;
  float y;
  float z;
};

struct Matrix4x4 {
  float m[4][4];
};

Vector3 Cross(const Vector3 &v1, const Vector3 &v2) {

  Vector3 result;

  result.x = v1.y * v2.z - v1.z * v2.y;
  result.y = v1.z * v2.x - v1.x * v2.z;
  result.z = v1.x * v2.y - v1.y * v2.x;

  return result;
}

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

  assert(w != 0);
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

// 正射影行列
Matrix4x4 MakeOrthographicMatrix(float left, float top, float right,
                                 float bottom, float nearClip, float farClip) {
  Matrix4x4 matrix;

  matrix.m[0][0] = 2.0f / (right - left);
  matrix.m[0][1] = 0.0f;
  matrix.m[0][2] = 0.0f;
  matrix.m[0][3] = 0.0f;
  matrix.m[1][0] = 0.0f;
  matrix.m[1][1] = 2.0f / (top - bottom);
  matrix.m[1][2] = 0.0f;
  matrix.m[1][3] = 0.0f;
  matrix.m[2][0] = 0.0f;
  matrix.m[2][1] = 0.0f;
  matrix.m[2][2] = 1.0f / (farClip / nearClip);
  matrix.m[2][3] = 0.0f;
  matrix.m[3][0] = (left + right) / (left - right);
  matrix.m[3][1] = (top + bottom) / (bottom - top);
  matrix.m[3][2] = nearClip / (nearClip - farClip);
  matrix.m[3][3] = 1.0f;

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

// 描画
static const int kRowHeight = 20;
static const int kColumnWidth = 60;

void MatrixScreenPrintf(int x, int y, const Matrix4x4 &matrix) {

  for (int row = 0; row < 4; ++row) {
    for (int column = 0; column < 4; ++column) {
      Novice::ScreenPrintf(x + column * kColumnWidth, y + row * kRowHeight,
                           "%6.02f", matrix.m[row][column]);
    }
  }
}

void VectorScreenPrintf(int x, int y, const Vector3 &vector,
                        const char *label) {
  Novice::ScreenPrintf(x, y, "%.02f", vector.x);
  Novice::ScreenPrintf(x + kColumnWidth, y, "%.02f", vector.y);
  Novice::ScreenPrintf(x + kColumnWidth * 2, y, "%.02f", vector.z);
  Novice::ScreenPrintf(x + kColumnWidth * 3, y, "%s", label);
};

const int kWindowWidth = 1280;
const int kWindowHeight = 720;

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

  // ライブラリの初期化
  Novice::Initialize(kWindowTitle, kWindowWidth, kWindowHeight);

  // キー入力結果を受け取る箱
  char keys[256] = {0};
  char preKeys[256] = {0};

  // Cross確認用
  Vector3 v1{1.2f, -3.9f, 2.5f};
  Vector3 v2{2.8f, 0.4f, -1.3f};
  Vector3 cross = Cross(v1, v2);

  // 三角形用
  Vector3 translate{};
  Vector3 rotate{};
  Vector3 cameraPosition{0.0f, 0.0f, 0.0f};
  Vector3 kLocalVertices[3];

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

    // 移動処理
    if (keys[DIK_W]) {
      translate.z++;
    } else if (keys[DIK_S]) {
      translate.z--;
    }

    if (keys[DIK_A]) {
      translate.x--;
    } else if (keys[DIK_D]) {
      translate.x++;
    }

    Matrix4x4 worldMatrix =
        MakeAffineMatrix({1.0f, 1.0f, 1.0f}, rotate, translate);
    Matrix4x4 cameraMatrix = MakeAffineMatrix(
        {1.0f, 1.0f, 1.0f}, {0.0f, 0.0f, 0.0f}, cameraPosition);
    Matrix4x4 viewMatrix = Inverse(cameraMatrix);
    Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(
        0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
    Matrix4x4 worldViewProjectionMatrix =
        Multiply(worldMatrix, Multiply(viewMatrix, projectionMatrix));
    Matrix4x4 viewportMatrix = MakeViewportMatrix(
        0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);

    Vector3 screenVertices[3];

    for (uint32_t i = 0; i < 3; ++i) {
      Vector3 ndcVertex =
          Transform(kLocalVertices[i], worldViewProjectionMatrix);
      screenVertices[i] = Transform(ndcVertex, viewMatrix);
    }

    ///
    /// ↑更新処理ここまで
    ///

    ///
    /// ↓描画処理ここから
    ///

    VectorScreenPrintf(0, 0, cross, "Cross");

    Novice::DrawTriangle(int(screenVertices[0].x), int(screenVertices[0].y),
                         int(screenVertices[1].x), int(screenVertices[1].y),
                         int(screenVertices[2].x), int(screenVertices[2].y),
                         RED, kFillModeSolid);

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
