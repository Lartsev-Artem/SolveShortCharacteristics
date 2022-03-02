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
	//int id;  - номер в массиве
	std::vector<Vector3> nodes_value;
	std::vector<int> neighbours_id_face;

	cell() {
		//id = -1;
		nodes_value.resize(4, Vector3(-666, -666, -666));	
	}
};


#define PI 3.14159265358979323846
const double eps = 1e-10;

extern Vector3 start_point_plane_coord;   // начало координат плоскости
extern Matrix3 transform_matrix;          // матрица перехода из базового тетраэдра в плоскость
extern Matrix3 inverse_transform_matrix;  // матрица перехода из плоскости в базовый тетраэдр

extern Matrix3	straight_face;  // 3 узла интерпол€ции
extern Matrix3 inclined_face;  // 3 узла интерпол€ции на наклонной плоскости

extern Matrix3	straight_face_inverse;  // 3 узла интерпол€ции
extern Matrix3 inclined_face_inverse;  // 3 узла интерпол€ции на наклонной плоскости

//std::vector<Type> Illum2;

extern int size_grid;

// скал€рные данные сетки (unstructured_grid)
extern vtkDataArray* density;
extern vtkDataArray* absorp_coef;
extern vtkDataArray* rad_en_loose_rate;

extern Type square_surface;  // площадь поверхности дискретной 

extern Vector3 center_local_sphere;  // центр описанной сферы около стандартного тетраэдра

extern int num_cur_direction; // номер текущего направлени€

extern int count_negative_interpolation; // число отрицательных значений интерпол€ции интесивности

extern std::vector<Type> res_inner_bound;  // значение на внутренней границе

extern std::vector<int>id_try_surface;


extern int posX;
extern int posX0;
extern int posOutC;
extern int posOut;
extern int posIn;
extern int posS;
extern int pos_in_res;
extern int id_try_pos;


#endif