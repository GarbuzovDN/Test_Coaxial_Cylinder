#pragma once
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <windows.h>

using namespace std;

extern double Pi;

extern double R0;
extern double R1;

/* ���������� ��������� */
extern int max_str, max_node, max_el;

/* ����� ���������� */
extern double Re;

/* ������� �������� */
extern int Iter_Glob;

/* ��������� ������� */
extern double omega_0;
extern double omega_1;

/* ����� �������� ��� ������������ */
extern int num_el_1, num_el_2, num_el_3;
extern int num_el_1_MUSCL;

/* ���������� ������������ �������� */
extern double xx_1, yy_1;

/* �������� ������������ */
extern double E_U_x, E_U_x_abs;
extern double E_U_y, E_U_y_abs;
extern int E_U_x_Num_el;
extern int E_U_y_Num_el;

/* ������������ �������� �������� */
extern double maxP_Corr;
extern int maxP_Cor_num_el;
extern double maxdivU;
extern int maxdivU_num_el;

/* ���������� ����� � ������ � Save */
extern string File_Mesh_Name;
extern ifstream File_Mesh;

extern bool Read_From_Save;
extern string File_Save_Name;

/* ��� � ������� ������� */
extern double dt;
extern double _time;
extern double final_time;

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

extern enum Topological_Attribute { Input = 2, Output = 4, Wall_1 = 7, Wall_2 = 3, Computational_Cell = 0 };

extern vector<Point> vectorPoint;
extern vector<Element> vectorElement;

extern ofstream Test_n;