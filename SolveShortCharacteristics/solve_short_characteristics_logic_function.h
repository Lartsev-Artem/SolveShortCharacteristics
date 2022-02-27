#pragma once
#ifndef SHORT_CHARACTERISTICS_LOGIC_H
#define SHORT_CHARACTERISTICS_LOGIC_H

#include "solve_short_characteristics_global_structure.h"
#include "solve_short_characteristics_headers.h"
#include "solve_short_characteristics_calculations.h"

Type CalculateIllumeOnInnerFace(const int num_cell, const int num_in_face, const int num_in_face_dir, const Vector3& x, const std::vector<cell>& grid);
int GetNodes(const int num_cur_cell, const Vector3& direction, std::vector<cell>& nodes_value, const std::vector<Type>& illum_old,
	const vector<Vector3>& directions, const vector<Type>& squares);

int GetNodes(const int num_cur_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell, const int num_cur_out_face,
	const Eigen::Matrix4d& vertex_tetra, const int* face_state, const Vector3& direction,
	std::vector<cell>& nodes_value, const std::vector<Type>& illum_old,
	const vector<Vector3>& directions, const vector<Type>& squares);

int CalculateNodeValue(const int num_cur_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell,
	const int num_cur_out_face, const int* face_state,
	const Vector3& direction, std::vector<cell>& nodes_value,
	const int num_node, const Eigen::Matrix4d& vertex_tetra, Vector3& x,
	const std::vector<Type>& illum_old, const vector<Vector3>& directions, const vector<Type>& squares);

Type CalculateIllumeOnInnerFace(const int num_cell, const int num_in_face, const Eigen::Matrix4d& vertex_tetra,
	const Vector3& x, const Vector3& x0, const std::vector<cell>& nodes_value);
Type CurGetIllum(const int cur_id, const Vector3 x, const Type s, const Type I_node_prev, const Vector3& cur_direction,
	const std::vector<Type>& illum_old,
	const vector<Vector3>& directions, const vector<Type>& squares);

Type GetValueInCenterCell(const int num_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell, const Vector3 center,
	const Vector3 direction,
	const Eigen::Matrix4d& vertex_tetra, const std::vector<cell>& nodes_value,
	const std::vector<Type>& illum_old, const vector<Vector3>& directions, const vector<Type>& squares);

#endif
