//----------------------------------------------------------------------------//
//             *** КВАТЕРНИОННО-ВЕКТОРНО-МАТРИЧНАЯ БИБЛИОТЕКА ***             //
//                                                                            //
// Файл GeneralLib.h                                                          //
//                                                                            //
// Общематематические сущности                                                //
//                                                                            //
//----------------------------------------------------------------------------//


#ifndef GeneralLibH
#define GeneralLibH

#include <math.h>


//--- Вектор трёхмерный ---
struct Vector
{
  public:
    Vector() : x(0.0), y(0.0), z(0.0) {}

    //--- Метаданные ---

    double x, y, z;

    //--- Перегруженные операции ---

    friend struct Vector operator+(struct Vector A, struct Vector B);  //Сложение двух векторов [ + ]
    friend struct Vector operator-(struct Vector A, struct Vector B);  //Вычитание двух векторов [ - ]
    friend double operator*(struct Vector A, struct Vector B);         //Скалярное умножение векторов [ * ]
    friend struct Vector operator^(struct Vector A, struct Vector B);  //Векторное умножение векторов [ ^ ]
    friend struct Vector operator*(struct Vector V, double S);         //Умножение вектора на скаляр [ * ]

    //--- Инициализирующие методы ---

    //Инициализирует вектор нулями
    void Empty(void);

    //Инициализирует вектор заданными значениями
    void Data(double ValueX, double ValueY, double ValueZ);

    //--- Свойства только для чтения ---

    //Возвращает модуль вектора
    double getModul(void);

    //Возвращает единичный вектор по направлению заданного
    struct Vector getBasis(void);
};
typedef struct Vector TVector;


//--- Матрица 3х3 ---
struct Matrix
{
  public:
    //--- Метаданные ---

    TVector FirstString, SecondString, ThirdString;

    //--- Перегруженные операции ---

    friend TVector operator*(struct Matrix M, TVector V);              //Произведение матрицы на вектор [ * ]
    friend struct Matrix operator+(struct Matrix A, struct Matrix B);  //Сложение двух матриц [ + ]
    friend struct Matrix operator-(struct Matrix A, struct Matrix B);  //Вычитание двух матриц [ - ]
    friend struct Matrix operator*(struct Matrix A, struct Matrix B);  //Произведение двух матриц [ * ]

    //--- Инициализирующие методы ---

    //Инициализирует матрицу нулями
    void Empty(void);

    //Инициализирует матрицу заданными значениями (запись по строкам)
    void Data(double a, double b, double c, double d, double e, double f, double g, double h, double k);

    //Инициализирует матрицу заданными значениями
    void DataVect(TVector FirstStr, TVector SecondStr, TVector ThirdStr);

    //Матрица вращения вокруг оси OX на заданный угол. Угол вращения считается положительным, если он виден из
    // конца оси, вокруг которой происходит вращение, против часовой стрелки
    void OXMatrix(double Angle);

    //Матрица вращения вокруг оси OY на заданный угол. Угол вращения считается положительным, если он виден из
    // конца оси, вокруг которой происходит вращение, против часовой стрелки
    void OYMatrix(double Angle);

    //Матрица вращения вокруг оси OZ на заданный угол. Угол вращения считается положительным, если он виден из
    // конца оси, вокруг которой происходит вращение, против часовой стрелки
    void OZMatrix(double Angle);

    //Получает матрицу вращения по заданным трем углам поворота вокруг осей OX, OY, OZ. О знаках углов поворота
    // смотри комментарии для функций OXMatrix, OYMatrix, OZMatrix
    // Order задает порядок поворота системы координат вокруг базовых осей
    // 1 - "xyz" - поворот вокруг оси OX, затем вокруг оси OY, затем вокруг оси OZ
    // 2 - "xzy" - поворот вокруг оси OX, затем вокруг оси OZ, затем вокруг оси OY
    // 3 - "yxz" - поворот вокруг оси OY, затем вокруг оси OX, затем вокруг оси OZ
    // 4 - "yzx" - поворот вокруг оси OY, затем вокруг оси OZ, затем вокруг оси OX
    // 5 - "zyx" - поворот вокруг оси OZ, затем вокруг оси OY, затем вокруг оси OX
    // 6 - "zxy" - поворот вокруг оси OZ, затем вокруг оси OX, затем вокруг оси OY
    void ThreeAnglesMatrix(double OXAngle, double OYAngle, double OZAngle, int Order);

    //--- Свойства только для чтения ---

    //Возвращает транспонированную матрицу
    struct Matrix getTrans(void);
};
typedef struct Matrix TMatrix;


//--- Кватернион ---
struct Quaternion
{
  public:
    Quaternion() : q0(0.0), q1(0.0), q2(0.0), q3(0.0) {}

    //--- Метаданные ---

    double q0;  //Скалярная составляющая w
    double q1;  //Векторная составляющая x
    double q2;  //Векторная составляющая y
    double q3;  //Векторная составляющая z

    //--- Перегруженные операции ---

    //Поворот вектора с помощью заданного кватерниона. Положительным считается поворот, выполненный по правилу
    // буравчика: большой палец указывает направление вектора кватерниона, а согнутая на 90 градусов ладонь
    // указывает направление положительного поворота (поворота на положительный угол). Произведение кватерниона на
    // вектор [ * ]
    friend TVector operator*(struct Quaternion Q, TVector V);
    friend struct Quaternion operator+(struct Quaternion A, struct Quaternion B);    //Сложение двух кватернионов [ + ]
    friend struct Quaternion operator-(struct Quaternion A, struct Quaternion B);    //Вычитание двух кватернионов [ - ]
    friend struct Quaternion operator*(struct Quaternion Q1, struct Quaternion Q2);  //Произведение двух кватернионов [ * ]

    //--- Инициализирующие методы ---

    //Инициализирует кватернион нулями
    void Empty(void);

    //Инициализирует кватернион заданными значениями
    void Data(double Valueq0, double Valueq1, double Valueq2, double Valueq3);

    //Создаёт кватернион вращения из вектора, вокруг которого надо выполнить вращение, и угла, на который надо
    // выполнить вращение. Знак величины угла определяется по правилу буравчика
    void MakeQuat(double Angle, TVector V);

    //Получает кватернион вращения из матрицы поворотов по алгоритму Стенли
    void QuatFromMatrixStanley(TMatrix M);

    //Получает кватернион по заданным трем углам поворота вокруг осей OX, OY, OZ. О знаках углов поворота
    // смотри комментарии для функций OXMatrix, OYMatrix, OZMatrix
    // Order задает порядок поворота системы координат вокруг базовых осей
    // 1 - "xyz" - поворот вокруг оси OX, затем вокруг оси OY, затем вокруг оси OZ
    // 2 - "xzy" - поворот вокруг оси OX, затем вокруг оси OZ, затем вокруг оси OY
    // 3 - "yxz" - поворот вокруг оси OY, затем вокруг оси OX, затем вокруг оси OZ
    // 4 - "yzx" - поворот вокруг оси OY, затем вокруг оси OZ, затем вокруг оси OX
    // 5 - "zyx" - поворот вокруг оси OZ, затем вокруг оси OY, затем вокруг оси OX
    // 6 - "zxy" - поворот вокруг оси OZ, затем вокруг оси OX, затем вокруг оси OY
    void ThreeAnglesQuat(double OXAngle, double OYAngle, double OZAngle, int Order);

    //--- Свойства только для чтения ---

    //Возвращает модуль угла поворота
    double getModulAngle(void);

    //Возвращает норму (модуль) кватерниона
    double getNorm(void);

    //Возвращает единичный кватернион по направлению заданного
    struct Quaternion getBasis(void);

    //Возвращает обратный кватернион
    struct Quaternion getOpposite(void);

    //Возвращает матрицу поворотов
    TMatrix getMatrix(void);
};
typedef struct Quaternion TQuaternion;


#endif
