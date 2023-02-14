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

/* ���������� ��������� */
int max_str, max_node, max_el;

/* ����� ���������� */
double Re = 1.0;

/* ������� �������� */
int Iter_Glob;

/* ��������� ������� */
double omega_0 = 0.0;
double omega_1 = 1.0;

/* ��������� ����� � ��������� ���������*/
struct Point
{

    /* �������������� ����� ���� */
    int Num_node;

    /*�������������� �������*/
    bool Boundary;

    /* ���������� ���� */
    double x, y;

    /* ��, ������� �������� ����� */
    vector<int> Neighb_el;

};
struct Element
{
    /*�������������� ����� �������� */
    int Num_el;

    /* ��� �������� */
    int Geom_el;

    /* ��� ������� */
    int Num_bound;

    /* ������ ������ �������� */
    int Num_vert[3];

    /* ���������� ������ �������� */
    Point Coord_vert[3];

    /* ����� ����� �������� */
    double Length_face_el[3];

    /* ������� �������� */
    double Area_el;

    /* ���������� ������ �������� */
    Point Coord_center_el;

    /* ���������� �� ������ �� ����� �������� */
    double h[3];

    /* �������� ������� */
    int Neighb_el[3];

    /* ������������ ������� */
    double Normal[3][2];

    /* �������� ��������� U_x, U_y, P */
    double gradU_x[2], gradU_y[2], gradP[2];

    /* �������� ����������� � �� �� ���������� ���� */
    double t;

    /* �������� ����������� � �� �� ������� ���� */
    double T;

    /* �������� ��. ����������� � �� �� ���������� ���� */
    double alfa;

    /* �������� ��. ����������� � �� �� ������� ���� */
    double Alfa;

    /* �������� ������. �� Alfa */
    double dalfa;

    /* �������� �������� �� ���������� ���� */
    double u_x, u_y;

    /* �������� �������� �� ������� ���� */
    double U_x, U_y;

    /* �������� �������� �� ������� ���� */
    double P;

    /* ������������ � ��� �������� �������� */
    double A_0;

    /* �������� �������� */
    double P_Correction;

    /* �������� �������� */
    double U_x_Correction;
    double U_y_Correction;

};

enum Topological_Attribute { Input = 2, Output = 4, Wall_1 = 7, Wall_2 = 3, Computational_Cell = 0 };

/* ������ ����� � ������ ��������� */
vector<Point> vectorPoint;
vector<Element> vectorElement;

/* ����� �������� ��� ������������ */
int num_el_1, num_el_2, num_el_3;
int num_el_1_MUSCL;

/* ���������� ������������ �������� */
double xx_1 = 0.169, yy_1 = 0.129;

/* �������� ������������ */
double E_U_x, E_U_x_abs;
double E_U_y, E_U_y_abs;
int E_U_x_Num_el;
int E_U_y_Num_el;

/* ������������ �������� �������� */
double maxP_Corr = 0.0;
int maxP_Cor_num_el = 0;
double maxdivU = 0.0;
int maxdivU_num_el = 0;

/* ���������� ����� � ������ � Save */
string File_Mesh_Name =
"Documents/Mesh/Mesh_Coaxial_Cylinders_(El=92).msh";
ifstream File_Mesh(File_Mesh_Name);

bool Read_From_Save = false;
string File_Save_Name =
"Documents/Save/Save_El = 2870/Save_(El = 1515).DAT";

/* ��� � ������� ������� */
double dt = 0.01;
double _time = 0.0;
double final_time = 1.0;

ofstream Test_n("Documents/Figure/n_Test.DAT", ios_base::trunc);
