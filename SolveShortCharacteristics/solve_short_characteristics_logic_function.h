#pragma once
#ifndef SHORT_CHARACTERISTICS_LOGIC_H
#define SHORT_CHARACTERISTICS_LOGIC_H

#include "solve_short_characteristics_global_structure.h"
#include "solve_short_characteristics_headers.h"
#include "solve_short_characteristics_calculations.h"

Type CalculateIllumeOnInnerFace(const int num_cell, const int num_in_face, const Vector3& x, const std::vector<Vector2>& X0,
	const std::vector<cell>& grid);

Type CurGetIllum(const int cur_id, const int cur_direction, const Vector3 x, const Type s, const Type I_node_prev,
	const vector<Type>& int_scattering);

Type GetValueInCenterCell(const int num_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell, const Vector3 center,
	const Vector3 direction,
	const Eigen::Matrix4d& vertex_tetra, const std::vector<cell>& nodes_value,
	const std::vector<Type>& illum_old, const vector<Vector3>& directions, const vector<Type>& squares);

#endif
