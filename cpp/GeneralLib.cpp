//----------------------------------------------------------------------------//
//             *** КВАТЕРНИОННО-ВЕКТОРНО-МАТРИЧНАЯ БИБЛИОТЕКА ***             //
//                                                                            //
// Файл GeneralLib.cpp                                                        //
//                                                                            //
// Общематематические сущности                                                //
//                                                                            //
//----------------------------------------------------------------------------//


#include "GeneralLib.h"


double TVector::getModul(void)
{
  try
  {
    return static_cast<double>(sqrt(x*x + y*y + z*z));
  }
  catch(...)
  {
    return 0.0;
  }
}


TVector TVector::getBasis(void)
{
double d;
TVector R;
  try
  {
    d = getModul();

    if (d != 0.0)
    {
      R.x = x / d;
      R.y = y / d;
      R.z = z / d;
    }
    else
    {
      R.x = x;
      R.y = y;
      R.z = z;
    }
  }
  catch(...)
  {
    R.Empty();
  }
  return R;
}


void TVector::Empty(void)
{
  x = 0.0;
  y = 0.0;
  z = 0.0;
}


void TVector::Data(double ValueX, double ValueY, double ValueZ)
{
  x = ValueX;
  y = ValueY;
  z = ValueZ;
}


TVector operator+(TVector A, TVector B)
{
TVector R;
  try
  {
    R.x = A.x + B.x;
    R.y = A.y + B.y;
    R.z = A.z + B.z;
  }
  catch(...)
  {
    R.Empty();
  }
  return R;
}


TVector operator-(TVector A, TVector B)
{
TVector R;
  try
  {
    R.x = A.x - B.x;
    R.y = A.y - B.y;
    R.z = A.z - B.z;
  }
  catch(...)
  {
    R.Empty();
  }
  return R;
}


double operator*(TVector A, TVector B)
{
  try
  {
    return static_cast<double>(A.x*B.x + A.y*B.y + A.z*B.z);
  }
  catch(...)
  {
    return 0.0;
  }
}


TVector operator^(TVector A, TVector B)
{
TVector R;
  try
  {
    R.x = A.y*B.z - A.z*B.y;
    R.y = A.z*B.x - A.x*B.z;
    R.z = A.x*B.y - A.y*B.x;
  }
  catch(...)
  {
    R.Empty();
  }
  return R;
}


TVector operator*(TVector V, double S)
{
TVector R;
  try
  {
    R.x = V.x * S;
    R.y = V.y * S;
    R.z = V.z * S;
  }
  catch(...)
  {
    R.Empty();
  }
  return R;
}


TMatrix TMatrix::getTrans(void)
{
TVector FirstColumn, SecondColumn, ThirdColumn;
TMatrix R;
  try
  {
    FirstColumn.x = FirstString.x;
    FirstColumn.y = SecondString.x;
    FirstColumn.z = ThirdString.x;

    SecondColumn.x = FirstString.y;
    SecondColumn.y = SecondString.y;
    SecondColumn.z = ThirdString.y;

    ThirdColumn.x = FirstString.z;
    ThirdColumn.y = SecondString.z;
    ThirdColumn.z = ThirdString.z;

    R.FirstString = FirstColumn;
    R.SecondString = SecondColumn;
    R.ThirdString = ThirdColumn;
  }
  catch(...)
  {
    R.Empty();
  }
  return R;
}


void TMatrix::Empty(void)
{
  FirstString.x = 0.0;
  FirstString.y = 0.0;
  FirstString.z = 0.0;
  SecondString.x = 0.0;
  SecondString.y = 0.0;
  SecondString.z = 0.0;
  ThirdString.x = 0.0;
  ThirdString.y = 0.0;
  ThirdString.z = 0.0;
}


void TMatrix::Data(double a, double b, double c, double d, double e, double f, double g, double h, double k)
{
  FirstString.x = a;
  FirstString.y = b;
  FirstString.z = c;
  SecondString.x = d;
  SecondString.y = e;
  SecondString.z = f;
  ThirdString.x = g;
  ThirdString.y = h;
  ThirdString.z = k;
}


void TMatrix::DataVect(TVector FirstStr, TVector SecondStr, TVector ThirdStr)
{
  FirstString = FirstStr;
  SecondString = SecondStr;
  ThirdString = ThirdStr;
}


void TMatrix::OXMatrix(double Angle)
{
  try
  {
    FirstString.x =  1.0;  FirstString.y =  0.0;         FirstString.z =  0.0;
    SecondString.x = 0.0;  SecondString.y = cos(Angle);  SecondString.z = -sin(Angle);
    ThirdString.x =  0.0;  ThirdString.y =  sin(Angle);  ThirdString.z =  cos(Angle);
  }
  catch(...)
  {
    FirstString.x = 0.0;
    FirstString.y = 0.0;
    FirstString.z = 0.0;
    SecondString.x = 0.0;
    SecondString.y = 0.0;
    SecondString.z = 0.0;
    ThirdString.x = 0.0;
    ThirdString.y = 0.0;
    ThirdString.z = 0.0;
  }
}


void TMatrix::OYMatrix(double Angle)
{
  try
  {
    FirstString.x =  cos(Angle);  FirstString.y =  0.0;  FirstString.z =  sin(Angle);
    SecondString.x = 0.0;         SecondString.y = 1.0;  SecondString.z = 0.0;
    ThirdString.x =  -sin(Angle); ThirdString.y =  0.0;  ThirdString.z =  cos(Angle);
  }
  catch(...)
  {
    FirstString.x = 0.0;
    FirstString.y = 0.0;
    FirstString.z = 0.0;
    SecondString.x = 0.0;
    SecondString.y = 0.0;
    SecondString.z = 0.0;
    ThirdString.x = 0.0;
    ThirdString.y = 0.0;
    ThirdString.z = 0.0;
  }
}


void TMatrix::OZMatrix(double Angle)
{
  try
  {
    FirstString.x =  cos(Angle);  FirstString.y =  -sin(Angle);  FirstString.z =  0.0;
    SecondString.x = sin(Angle);  SecondString.y = cos(Angle);   SecondString.z = 0.0;
    ThirdString.x =  0.0;         ThirdString.y =  0.0;          ThirdString.z =  1.0;
  }
  catch(...)
  {
    FirstString.x = 0.0;
    FirstString.y = 0.0;
    FirstString.z = 0.0;
    SecondString.x = 0.0;
    SecondString.y = 0.0;
    SecondString.z = 0.0;
    ThirdString.x = 0.0;
    ThirdString.y = 0.0;
    ThirdString.z = 0.0;
  }
}


void TMatrix::ThreeAnglesMatrix(double OXAngle, double OYAngle, double OZAngle, int Order)
{
TMatrix M, M1, M2, M3, M4;
  try
  {
    M.FirstString.x = 0.0;  M.FirstString.y = 0.0;  M.FirstString.z = 0.0;
    M.SecondString.x = 0.0; M.SecondString.y = 0.0; M.SecondString.z = 0.0;
    M.ThirdString.x = 0.0;  M.ThirdString.y = 0.0;  M.ThirdString.z = 0.0;

    OXMatrix(OXAngle);

    M1.FirstString.x  = FirstString.x;
    M1.FirstString.y  = FirstString.y;
    M1.FirstString.z  = FirstString.z;
    M1.SecondString.x = SecondString.x;
    M1.SecondString.y = SecondString.y;
    M1.SecondString.z = SecondString.z;
    M1.ThirdString.x  = ThirdString.x;
    M1.ThirdString.y  = ThirdString.y;
    M1.ThirdString.z  = ThirdString.z;

    OYMatrix(OYAngle);

    M2.FirstString.x  = FirstString.x;
    M2.FirstString.y  = FirstString.y;
    M2.FirstString.z  = FirstString.z;
    M2.SecondString.x = SecondString.x;
    M2.SecondString.y = SecondString.y;
    M2.SecondString.z = SecondString.z;
    M2.ThirdString.x  = ThirdString.x;
    M2.ThirdString.y  = ThirdString.y;
    M2.ThirdString.z  = ThirdString.z;

    OZMatrix(OZAngle);

    M3.FirstString.x  = FirstString.x;
    M3.FirstString.y  = FirstString.y;
    M3.FirstString.z  = FirstString.z;
    M3.SecondString.x = SecondString.x;
    M3.SecondString.y = SecondString.y;
    M3.SecondString.z = SecondString.z;
    M3.ThirdString.x  = ThirdString.x;
    M3.ThirdString.y  = ThirdString.y;
    M3.ThirdString.z  = ThirdString.z;

    if (Order == 1)
    {
      M4 = M3*M2;
      M = M4*M1;
    }

    if (Order == 2)
    {
      M4 = M2*M3;
      M = M4*M1;
    }

    if (Order == 3)
    {
      M4 = M3*M1;
      M = M4*M2;
    }

    if (Order == 4)
    {
      M4 = M1*M3;
      M = M4*M2;
    }

    if (Order == 5)
    {
      M4 = M1*M2;
      M = M4*M3;
    }

    if (Order == 6)
    {
      M4 = M2*M1;
      M = M4*M3;
    }

    FirstString.x =  M.FirstString.x;
    FirstString.y =  M.FirstString.y;
    FirstString.z =  M.FirstString.z;
    SecondString.x = M.SecondString.x;
    SecondString.y = M.SecondString.y;
    SecondString.z = M.SecondString.z;
    ThirdString.x =  M.ThirdString.x;
    ThirdString.y =  M.ThirdString.y;
    ThirdString.z =  M.ThirdString.z;
  }
  catch(...)
  {
    FirstString.x = 0.0;
    FirstString.y = 0.0;
    FirstString.z = 0.0;
    SecondString.x = 0.0;
    SecondString.y = 0.0;
    SecondString.z = 0.0;
    ThirdString.x = 0.0;
    ThirdString.y = 0.0;
    ThirdString.z = 0.0;
  }
}


TVector operator*(TMatrix M, TVector V)
{
TVector R;
  try
  {
    R.x = M.FirstString.x*V.x + M.FirstString.y*V.y + M.FirstString.z*V.z;
    R.y = M.SecondString.x*V.x + M.SecondString.y*V.y + M.SecondString.z*V.z;
    R.z = M.ThirdString.x*V.x + M.ThirdString.y*V.y + M.ThirdString.z*V.z;
  }
  catch(...)
  {
    R.Empty();
  }
  return R;
}


TMatrix operator+(TMatrix A, TMatrix B)
{
TMatrix R;
  try
  {
    R.FirstString  = A.FirstString + B.FirstString;
    R.SecondString = A.SecondString + B.SecondString;
    R.ThirdString  = A.ThirdString + B.ThirdString;
  }
  catch(...)
  {
    R.Empty();
  }
  return R;
}


TMatrix operator-(TMatrix A, TMatrix B)
{
TMatrix R;
  try
  {
    R.FirstString  = A.FirstString - B.FirstString;
    R.SecondString = A.SecondString - B.SecondString;
    R.ThirdString  = A.ThirdString - B.ThirdString;
  }
  catch(...)
  {
    R.Empty();
  }
  return R;
}


TMatrix operator*(TMatrix A, TMatrix B)
{
TVector FirstColumn, SecondColumn, ThirdColumn;
TMatrix R;
  try
  {
    FirstColumn.x = B.FirstString.x;
    FirstColumn.y = B.SecondString.x;
    FirstColumn.z = B.ThirdString.x;

    SecondColumn.x = B.FirstString.y;
    SecondColumn.y = B.SecondString.y;
    SecondColumn.z = B.ThirdString.y;

    ThirdColumn.x = B.FirstString.z;
    ThirdColumn.y = B.SecondString.z;
    ThirdColumn.z = B.ThirdString.z;

    R.FirstString.x =  A.FirstString * FirstColumn;
    R.SecondString.x = A.SecondString * FirstColumn;
    R.ThirdString.x =  A.ThirdString * FirstColumn;

    R.FirstString.y =  A.FirstString * SecondColumn;
    R.SecondString.y = A.SecondString * SecondColumn;
    R.ThirdString.y =  A.ThirdString * SecondColumn;

    R.FirstString.z =  A.FirstString * ThirdColumn;
    R.SecondString.z = A.SecondString * ThirdColumn;
    R.ThirdString.z =  A.ThirdString * ThirdColumn;
  }
  catch(...)
  {
    R.Empty();
  }
  return R;
}


double TQuaternion::getModulAngle(void)
{
double Angle;
TQuaternion P;
  try
  {
    P = getBasis();
    Angle = acos(P.q0)*2.0;
  }
  catch(...)
  {
    Angle = 0.0;
  }
  return Angle;
}


double TQuaternion::getNorm(void)
{
double N;
  try
  {
    N = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
  }
  catch(...)
  {
    N = 0.0;
  }
  return N;
}


TQuaternion TQuaternion::getBasis(void)
{
TQuaternion R;
double N;
  try
  {
    N = getNorm();

    R.q0 = q0 / N;
    R.q1 = q1 / N;
    R.q2 = q2 / N;
    R.q3 = q3 / N;
  }
  catch(...)
  {
    R.Empty();
  }
  return R;
}


TQuaternion TQuaternion::getOpposite(void)
{
TQuaternion P;
double N;
  try
  {
    N = getNorm();
    P.q0 = q0 / N;
    P.q1 = -q1 / N;
    P.q2 = -q2 / N;
    P.q3 = -q3 / N;
  }
  catch(...)
  {
    P.Empty();
  }
  return P;
}


TMatrix TQuaternion::getMatrix(void)
{
double wx, wy, wz, xx, yy, yz, xy, xz, zz, x2, y2, z2;
TMatrix M;
  try
  {
    x2 = q1 + q1;
    y2 = q2 + q2;
    z2 = q3 + q3;

    xx = q1 * x2;
    yy = q2 * y2;
    wx = q0 * x2;

    xy = q1 * y2;
    yz = q2 * z2;
    wy = q0 * y2;

    xz = q1 * z2;
    zz = q3 * z2;
    wz = q0 * z2;

    M.FirstString.x  = 1.0-(yy+zz);
    M.SecondString.x = xy+wz;
    M.ThirdString.x  = xz-wy;

    M.FirstString.y  = xy-wz;
    M.SecondString.y = 1.0-(xx+zz);
    M.ThirdString.y  = yz+wx;

    M.FirstString.z  = xz+wy;
    M.SecondString.z = yz-wx;
    M.ThirdString.z  = 1.0-(xx+yy);
  }
  catch(...)
  {
    M.Empty();
  }
  return M;
}


TVector operator*(TQuaternion Q, TVector V)
{
TVector R;
double A, B, C, D, E, F, G, H;
TQuaternion P;
  try
  {
    P = Q.getOpposite();

    A = (Q.q0 + Q.q1) * (V.x);
    B = (Q.q3 - Q.q2) * (V.y - V.z);
    C = (Q.q1 - Q.q0) * (V.y + V.z);
    D = (Q.q2 + Q.q3) * (V.x);
    E = (Q.q1 + Q.q3) * (V.x + V.y);
    F = (Q.q1 - Q.q3) * (V.x - V.y);
    G = (Q.q0 + Q.q2) * (-V.z);
    H = (Q.q0 - Q.q2) * (V.z);

    Q.q0 =  B + (-E - F + G + H) * 0.5;
    Q.q1 =  A - ( E + F + G + H) * 0.5;
    Q.q2 = -C + ( E - F + G - H) * 0.5;
    Q.q3 = -D + ( E - F - G + H) * 0.5;

    A = (Q.q0 + Q.q1) * (P.q0 + P.q1);
    C = (Q.q1 - Q.q0) * (P.q2 + P.q3);
    D = (Q.q2 + Q.q3) * (P.q1 - P.q0);
    E = (Q.q1 + Q.q3) * (P.q1 + P.q2);
    F = (Q.q1 - Q.q3) * (P.q1 - P.q2);
    G = (Q.q0 + Q.q2) * (P.q0 - P.q3);
    H = (Q.q0 - Q.q2) * (P.q0 + P.q3);

    R.x =  A - ( E + F + G + H) * 0.5;
    R.y = -C + ( E - F + G - H) * 0.5;
    R.z = -D + ( E - F - G + H) * 0.5;
  }
  catch(...)
  {
    R.Empty();
  }
  return R;
}


TQuaternion operator+(TQuaternion A, TQuaternion B)
{
TQuaternion R;
  try
  {
    R.q0 = A.q0 + B.q0;
    R.q1 = A.q1 + B.q1;
    R.q2 = A.q2 + B.q2;
    R.q3 = A.q3 + B.q3;
  }
  catch(...)
  {
    R.Empty();
  }
  return R;
}


TQuaternion operator-(TQuaternion A, TQuaternion B)
{
TQuaternion R;
  try
  {
    R.q0 = A.q0 - B.q0;
    R.q1 = A.q1 - B.q1;
    R.q2 = A.q2 - B.q2;
    R.q3 = A.q3 - B.q3;
  }
  catch(...)
  {
    R.Empty();
  }
  return R;
}


TQuaternion operator*(TQuaternion Q1, TQuaternion Q2)
{
TQuaternion R;
double A, B, C, D, E, F, G, H;
  try
  {
    A = (Q1.q0 + Q1.q1) * (Q2.q0 + Q2.q1);
    B = (Q1.q3 - Q1.q2) * (Q2.q2 - Q2.q3);
    C = (Q1.q1 - Q1.q0) * (Q2.q2 + Q2.q3);
    D = (Q1.q2 + Q1.q3) * (Q2.q1 - Q2.q0);
    E = (Q1.q1 + Q1.q3) * (Q2.q1 + Q2.q2);
    F = (Q1.q1 - Q1.q3) * (Q2.q1 - Q2.q2);
    G = (Q1.q0 + Q1.q2) * (Q2.q0 - Q2.q3);
    H = (Q1.q0 - Q1.q2) * (Q2.q0 + Q2.q3);

    R.q0 = B + (-E - F + G + H) * 0.5;
    R.q1 = A - ( E + F + G + H) * 0.5;
    R.q2 = -C + ( E - F + G - H) * 0.5;
    R.q3 = -D + ( E - F - G + H) * 0.5;
  }
  catch(...)
  {
    R.Empty();
  }
  return R;
}


void TQuaternion::Empty(void)
{
  q0 = 0.0;
  q1 = 0.0;
  q2 = 0.0;
  q3 = 0.0;
}


void TQuaternion::Data(double Valueq0, double Valueq1, double Valueq2, double Valueq3)
{
  q0 = Valueq0;
  q1 = Valueq1;
  q2 = Valueq2;
  q3 = Valueq3;
}


void TQuaternion::MakeQuat(double Angle, TVector V)
{
TVector V_;
TQuaternion Q;
  try
  {
    //Вычисляю единичный вектор в направлении заданного вектора
    V_ = V.getBasis();

    //Создаю кватернион (пока он не только вращает, но и растягивает пространство)
    Q.q0 = cos(Angle/2.0);
    Q.q1 = V_.x * sin(Angle/2.0);
    Q.q2 = V_.y * sin(Angle/2.0);
    Q.q3 = V_.z * sin(Angle/2.0);

    //Привожу кватернион к единице (теперь он не растягивает пространство)
    Q = Q.getBasis();

    q0 = Q.q0;
    q1 = Q.q1;
    q2 = Q.q2;
    q3 = Q.q3;
  }
  catch(...)
  {
    q0 = 0.0;
    q1 = 0.0;
    q2 = 0.0;
    q3 = 0.0;
  }
}


void TQuaternion::QuatFromMatrixStanley(TMatrix M)
{
TQuaternion Q;
double T;
double Valueq0, Valueq1, Valueq2, Valueq3;
  try
  {
    Q.q0 = 0.0;
    Q.q1 = 0.0;
    Q.q2 = 0.0;
    Q.q3 = 0.0;

    T = M.FirstString.x + M.SecondString.y + M.ThirdString.z;

    Valueq0 = (1.0 + T)/4.0;
    Valueq1 = (1.0 + 2.0*M.FirstString.x - T)/4.0;
    Valueq2 = (1.0 + 2.0*M.SecondString.y - T)/4.0;
    Valueq3 = (1.0 + 2.0*M.ThirdString.z - T)/4.0;

    if ((Valueq0 >= Valueq1) && (Valueq0 >= Valueq2) && (Valueq0 >= Valueq3))
    {
      Q.q0 = sqrt(Valueq0);
      Q.q1 = (M.SecondString.z - M.ThirdString.y)  / (4.0 * Q.q0);
      Q.q2 = (M.ThirdString.x  - M.FirstString.z)  / (4.0 * Q.q0);
      Q.q3 = (M.FirstString.y  - M.SecondString.x) / (4.0 * Q.q0);
    }

    if ((Valueq1 >= Valueq0) && (Valueq1 >= Valueq2) && (Valueq1 >= Valueq3))
    {
      Q.q1 = sqrt(Valueq1);
      Q.q0 = (M.SecondString.z - M.ThirdString.y)  / (4.0 * Q.q1);
      Q.q2 = (M.FirstString.y  + M.SecondString.x) / (4.0 * Q.q1);
      Q.q3 = (M.ThirdString.x  + M.FirstString.z)  / (4.0 * Q.q1);
    }

    if ((Valueq2 >= Valueq0) && (Valueq2 >= Valueq1) && (Valueq2 >= Valueq3))
    {
      Q.q2 = sqrt(Valueq2);
      Q.q0 = (M.ThirdString.x - M.FirstString.z) / (4.0 * Q.q2);
      Q.q1 = (M.FirstString.y + M.SecondString.x) / (4.0 * Q.q2);
      Q.q3 = (M.SecondString.z + M.ThirdString.y) / (4.0 * Q.q2);
    }

    if ((Valueq3 >= Valueq1) && (Valueq3 >= Valueq2) && (Valueq3 >= Valueq0))
    {
      Q.q3 = sqrt(Valueq3);
      Q.q0 = (M.FirstString.y  - M.SecondString.x) / (4.0 * Q.q3);
      Q.q1 = (M.ThirdString.x  + M.FirstString.z)  / (4.0 * Q.q3);
      Q.q2 = (M.SecondString.z + M.ThirdString.y)  / (4.0 * Q.q3);
    }

    Q = Q.getBasis();

    Q = Q.getOpposite();

    q0 = Q.q0;
    q1 = Q.q1;
    q2 = Q.q2;
    q3 = Q.q3;
  }
  catch(...)
  {
    q0 = 0.0;
    q1 = 0.0;
    q2 = 0.0;
    q3 = 0.0;
  }
}


void TQuaternion::ThreeAnglesQuat(double OXAngle, double OYAngle, double OZAngle, int Order)
{
TQuaternion Q, P1, P2, P3, P4;
TVector V;
  try
  {
    Q.q0 = 0.0;
    Q.q1 = 0.0;
    Q.q2 = 0.0;
    Q.q3 = 0.0;

    V.x = 1.0;
    V.y = 0.0;
    V.z = 0.0;
    P1.MakeQuat(OXAngle,V);

    V.x = 0.0;
    V.y = 1.0;
    V.z = 0.0;
    P2.MakeQuat(OYAngle,V);

    V.x = 0.0;
    V.y = 0.0;
    V.z = 1.0;
    P3.MakeQuat(OZAngle,V);


    if (Order == 1)
    {
      P4 = P3*P2;
      Q = P4*P1;
    }

    if (Order == 2)
    {
      P4 = P2*P3;
      Q = P4*P1;
    }

    if (Order == 3)
    {
      P4 = P3*P1;
      Q = P4*P2;
    }

    if (Order == 4)
    {
      P4 = P1*P3;
      Q = P4*P2;
    }

    if (Order == 5)
    {
      P4 = P1*P2;
      Q = P4*P3;
    }

    if (Order == 6)
    {
      P4 = P2*P1;
      Q = P4*P3;
    }

    q0 = Q.q0;
    q1 = Q.q1;
    q2 = Q.q2;
    q3 = Q.q3;
  }
  catch(...)
  {
    q0 = 0.0;
    q1 = 0.0;
    q2 = 0.0;
    q3 = 0.0;
  }
}
