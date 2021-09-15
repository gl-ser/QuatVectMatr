# QuatVectMatr
Кватернионно-векторно-матричная библиотека

---

```
TVector трёхмерный вектор, имеющий координаты x, y, z
Методы TVector:
	- Empty инициализирует вектор нулями;
	- Data инициализирует вектор заданными значениями;
	- getModul возвращает модуль вектора;
	- getBasis возвращает единичный вектор по направлению заданного.
У TVector имеются перегруженные операции:
	- [+] – сложение двух векторов;
	- [–] – вычитание двух векторов;
	- [*] – скалярное умножение векторов;
	- [^] – векторное умножение векторов;
	- [*] – умножение вектора на скаляр.

TMatrix матрица 3x3, состоящая из трёх объектов типа TVector
FirstString, SecondString, ThirdString – построчные элементы матрицы.
Здесь и далее угол вращения считается положительным, если он виден из конца 
оси, вокруг которой происходит вращение, против часовой стрелки.
Методы TMatrix:
	- Empty инициализирует матрицу нулями; 
	- Data инициализирует матрицу заданными значениями (запись по строкам);
	- DataVect инициализирует матрицу заданными значениями (запись векторами);
	- OXMatrix получает матрицу вращения вокруг оси OX на заданный угол;
	- OYMatrix получает матрицу вращения вокруг оси OY на заданный угол;
	- OZMatrix получает матрицу вращения вокруг оси OZ на заданный угол;
	- ThreeAnglesMatrix получает матрицу вращения по заданным трём углам
          поворота вокруг осей OX, OY, OZ и порядку поворота системы координат 
          вокруг базовых осей;
	- getTrans возвращает транспонированную матрицу.
У TMatrix имеются перегруженные операции:
	- [+] – сложение двух матриц;
	- [–] – вычитание двух матриц;
	- [*] – произведение матрицы на вектор;
	- [*] – произведение двух матриц.

TQuaternion кватернион, имеющий атрибуты q0, q1, q2, q3
q0 – скалярная часть, q1, q2, q3 – векторная часть.
Методы TQuaternion:
	- Empty инициализирует кватернион нулями; 
	- Data инициализирует кватернион заданными значениями;
	- MakeQuat создаёт кватернион вращения из вектора, вокруг которого надо 
          выполнить вращение, и угла, на который надо выполнить вращение;
	- QuatFromMatrixStanley получает кватернион вращения из матрицы поворотов 
          по алгоритму Стенли;
	- ThreeAnglesQuat получает кватернион по заданным трём углам поворота вокруг 
          осей OX, OY, OZ и порядку поворота системы координат вокруг базовых осей;
	- getModulAngle возвращает модуль угла поворота;
	- getNorm возвращает норму (модуль) кватерниона;
	- getBasis возвращает единичный кватернион по направлению заданного;
	- getOpposite возвращает обратный кватернион;
	- getMatrix возвращает матрицу поворотов.
У TQuaternion имеются перегруженные операции:
	- [+] – сложение двух кватернионов;
	- [–] – вычитание двух кватернионов;
	- [*] – произведение двух кватернионов;
	- [*] – произведение кватерниона на вектор (поворот вектора с помощью заданного 
          кватерниона по правилу буравчика).
```

---

См. [аналог](https://github.com/gl-ser/QuatVectMatr_Pas) для Object Pascal
