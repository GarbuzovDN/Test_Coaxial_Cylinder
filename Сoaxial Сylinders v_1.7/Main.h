#pragma once

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <windows.h>

using namespace std;

double Pi = 3.14159;

double R0 = 0.2;
double R1 = 1.0;

/* Количество элементов */
int max_str, max_node, max_el;

/* Число Рейнольдса */
double Re = 1.0;

/* Счетчик итераций */
int Iter_Glob;

/* Граничные условия */
double omega_0 = 0.0;
double omega_1 = 1.0;

/* Структура точек и структура элементов*/
struct Point
{

    /* Индивидуальный номер узла */
    int Num_node;

    /*Принадлежность границе*/
    bool Boundary;

    /* Координаты узла */
    double x, y;

    /* КО, которые содержат точку */
    vector<int> Neighb_el;

};
struct Element
{
    /*Индивидуальный номер элемента */
    int Num_el;

    /* Тип элемента */
    int Geom_el;

    /* Тег границы */
    int Num_bound;

    /* Номера вершин элемента */
    int Num_vert[3];

    /* Координаты вершин элемента */
    Point Coord_vert[3];

    /* Длина грани элемента */
    double Length_face_el[3];

    /* Площадь элемента */
    double Area_el;

    /* Координаты центра элемента */
    Point Coord_center_el;

    /* Расстояние от центра до грани элемента */
    double h[3];

    /* Соседний элемент */
    int Neighb_el[3];

    /* Составляющие нормали */
    double Normal[3][2];

    /* Значение градиента U_x, U_y, P */
    double gradU_x[2], gradU_y[2], gradP[2];

    /* Значение температуры в КО на предыдущем шаге */
    double t;

    /* Значение температуры в КО на текущем шаге */
    double T;

    /* Значение ст. отверждения в КО на предыдущем шаге */
    double alfa;

    /* Значение ст. отверждения в КО на текущем шаге */
    double Alfa;

    /* Значение произв. от Alfa */
    double dalfa;

    /* Значение скорости на предыдущем шаге */
    double u_x, u_y;

    /* Значение скорости на текущем шаге */
    double U_x, U_y;

    /* Значение давления на текущем шаге */
    double P;

    /* Коэффициенты А для поправки давления */
    double A_0;

    /* Поправка давления */
    double P_Correction;

    /* Поправка скорости */
    double U_x_Correction;
    double U_y_Correction;

};

enum Topological_Attribute { Input = 2, Output = 4, Wall_1 = 7, Wall_2 = 3, Computational_Cell = 0 };

/* Вектор точек и вектор элементов */
vector<Point> vectorPoint;
vector<Element> vectorElement;

/* Номер элемента для отслеживания */
int num_el_1, num_el_2, num_el_3;
int num_el_1_MUSCL;

/* Координаты контрольного элемента */
double xx_1 = 0.169, yy_1 = 0.129;

/* Параметр установления */
double E_U_x, E_U_x_abs;
double E_U_y, E_U_y_abs;
int E_U_x_Num_el;
int E_U_y_Num_el;

/* Максимальная поправка давления */
double maxP_Corr = 0.0;
int maxP_Cor_num_el = 0;
double maxdivU = 0.0;
int maxdivU_num_el = 0;

/* Директория файла с сеткой и Save */
string File_Mesh_Name =
"Documents/Mesh/Mesh_Coaxial_Cylinders_(El=92).msh";
ifstream File_Mesh(File_Mesh_Name);

bool Read_From_Save = false;
string File_Save_Name =
"Documents/Save/Save_El = 2870/Save_(El = 1515).DAT";

/* Шаг и счетчик времени */
double dt = 0.01;
double _time = 0.0;
double final_time = 1.0;

ofstream Test_n("Documents/Figure/n_Test.DAT", ios_base::trunc);
