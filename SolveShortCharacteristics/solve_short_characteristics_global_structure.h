#pragma once
#ifndef SHORT_CHARACTERISTICS_GLOBAL_H
#define SHORT_CHARACTERISTICS_GLOBAL_H

#include "solve_short_characteristics_headers.h"

typedef int IntId;
typedef double Type;
typedef std::string Str_Type;
typedef Eigen::Vector3d Vector3;
typedef Eigen::Vector2d Vector2;
typedef Eigen::Matrix3d Matrix3;


#ifdef _MSC_VER
#define fwrite_unlocked _fwrite_nolock
#define fread_unlocked  _fread_nolock
#endif

#define PRINTF(a) {printf(a); return 1;}

using namespace std::chrono;
using namespace std;

struct Normals {
	std::vector<Vector3> n;
	Normals() {
	}

	Normals(const int size) {
		n.resize(size);
	}
};

struct cell {
	//int id;  - ����� � �������
	std::vector<Vector3> nodes_value;
	std::vector<int> neighbours_id_face;

	std::vector<int> in_id;   // �� 9 ����� �� ������ �����������
	std::vector<int> out_id;  // �� 3 ����� �� ������ �����������
	std::vector<Type> s;  // �� 9 ����� �� ������ ����������� // init = -1

	std::vector<Vector2> x0_loc;
	std::vector<Vector3> x;

	std::vector<int> out_count;  // { (0 ��� 1 ��� 2)  {����� �����������}

	int start_dir1;
	int start_dir3;

	cell() {
		//id = -1;
		nodes_value.resize(4, Vector3(-666, -666, -666));	
		start_dir1 = 0;
		start_dir3 = 0;
	}
};

#define PI 3.14159265358979323846
const double eps = 1e-10;

extern Vector3 start_point_plane_coord;   // ������ ��������� ���������
extern Matrix3 transform_matrix;          // ������� �������� �� �������� ��������� � ���������
extern Matrix3 inverse_transform_matrix;  // ������� �������� �� ��������� � ������� ��������

extern Matrix3	straight_face;  // 3 ���� ������������
extern Matrix3 inclined_face;  // 3 ���� ������������ �� ��������� ���������

extern Matrix3	straight_face_inverse;  // 3 ���� ������������
extern Matrix3 inclined_face_inverse;  // 3 ���� ������������ �� ��������� ���������

//std::vector<Type> Illum2;

extern int size_grid;

// ��������� ������ ����� (unstructured_grid)
extern vtkDataArray* density;
extern vtkDataArray* absorp_coef;
extern vtkDataArray* rad_en_loose_rate;

extern Type square_surface;  // ������� ����������� ���������� 

extern Vector3 center_local_sphere;  // ����� ��������� ����� ����� ������������ ���������

extern int num_cur_direction; // ����� �������� �����������

extern int count_negative_interpolation; // ����� ������������� �������� ������������ ������������

extern std::vector<Type> res_inner_bound;  // �������� �� ���������� �������

extern std::vector<int>id_try_surface;

#endif