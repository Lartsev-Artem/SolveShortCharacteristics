#include "solve_short_characteristics_logic_function.h"

Type CalculateIllumeOnInnerFace(const int num_cell, const int num_in_face, const int num_in_face_dir ,const Vector3& x, const std::vector<cell>& grid) {
	Type I_x0 = 0;

	if (grid[num_cell].neighbours_id_face[num_in_face] == -1) {
		/*Граничные условия*/
		//I_x0 = BoundaryFunction(num_cell, x, direction, illum_old, directions, squares);
		return I_x0;
	}
	else if (grid[num_cell].neighbours_id_face[num_in_face] == -2) // внутренняя граница
	{
		static int pos_in_res = 0;
		static int id_try_pos = 0;
		id_try_pos++;

		Type data = res_inner_bound[pos_in_res++]; // защита на выход из диапазона??
		if (data >= 0) //данные от пересечения с диском или шаром
		{
		/*
			 Проверить порядок. Данные в массиве должны лежать в соответствии с упорядоченным графом 
			 от направления к направлению
		*/
			return data;  
		}
		else // определяющимм являются противолежаащие грани (возможен расчет с учетом s=(x-x0).norm())
		{

			// результат не вполне понятен. Пока лучше использовать константу или другие параметры области (шар == граница)

			//+dist_try_surface			
			int id = id_try_surface[id_try_pos - 1];  // будет лежать id грани			
			const int cell = id / 4;
			const int face = id % 4;		
			
			Vector3 coef = grid[cell].nodes_value[face];
			Vector2	x0_local = grid[num_cell].x0_loc[num_in_face_dir];
			
			I_x0 = x0_local[0] * coef[0] + x0_local[1] * coef[1] + coef[2];
			
			//if (I_x0 < 0) I_x0 = 0;
			return 10; // I_x0;

		}
	}else {

		// debug условие
		if (grid[num_cell].nodes_value[num_in_face][0] < -600)
			cout << "Num_dir: " << num_cur_direction << " CalculateIllumeOnInnerFace:  Undefine cell / in face:" << num_cell << " / " << num_in_face << " !!!\n";


		Vector3 coef = grid[num_cell].nodes_value[num_in_face];
		Vector2	x0_local = grid[num_cell].x0_loc[num_in_face_dir];
		I_x0 = x0_local[0] * coef[0] + x0_local[1] * coef[1] + coef[2];
		
		/*
		switch (num_in_face) {
		case 3:
			Vector3 local_plane_x0;
			FromTetraToPlane(transform_matrix, start_point_plane_coord, x0_local, local_plane_x0);

			coef = GetInterpolationCoef(inclined_face, nodes_value[num_cell].nodes_value[num_in_face]);
			I_x0 = local_plane_x0[0] * coef[0] + local_plane_x0[1] * coef[1] + coef[2];

			//I_x0 = x0_local[1] * coef[0] + x0_local[2] * coef[1] + coef[2];
			//Vector3 coef = GetInterpolationCoef(straight_face, nodes_value.find(global_num_in_face)->second); 
			break;
		case 1:
			coef = GetInterpolationCoef(straight_face, nodes_value[num_cell].nodes_value[num_in_face]);
			I_x0 = x0_local[1] * coef[0] + x0_local[2] * coef[1] + coef[2];

			break;
		case 2:
			coef = GetInterpolationCoef(straight_face, nodes_value[num_cell].nodes_value[num_in_face]);
			I_x0 = x0_local[0] * coef[0] + x0_local[2] * coef[1] + coef[2];

			break;
		case 0:
			coef = GetInterpolationCoef(straight_face, nodes_value[num_cell].nodes_value[num_in_face]);
			I_x0 = x0_local[0] * coef[0] + x0_local[1] * coef[1] + coef[2];

			break;
		}
		*/

		if (I_x0 < 0) {
			count_negative_interpolation++;
			return 0;
		}

		return I_x0;
	}
}

int CalculateNodeValue(const int num_cur_cell, const int num_node, const Vector3& direction,
	std::vector<cell>& nodes_value,  Vector3& x, const std::vector<Type>& illum_old, 
	const vector<Vector3>& directions, const vector<Type>& squares) {

	Vector3 x0;
	size_t num_in_face;
	//for (size_t num_in_face = 0; num_in_face < 4; ++num_in_face) {
	//	if (!face_state[num_in_face]) continue;  // обрабатываем только входные грани


	//	IntersectionWithPlane(cur_cell->GetFace(num_in_face), x, direction, x0);

	//	if (InTriangle(num_cur_cell, unstructuredgrid, cur_cell, num_in_face, x0))
		{

			// значение на входящей грани
		Type I_x0 = 0;// CalculateIllumeOnInnerFace(num_cur_cell, num_in_face, x, nodes_value);

			Type s = nodes_value[num_cur_cell].s[num_in_face];

			Type I = CurGetIllum(num_cur_cell, x0, s, I_x0, direction, illum_old, directions, squares);

			//nodes_value[num_cur_cell].nodes_value[num_cur_out_face][num_node] = I;

			//break;
		}

	//}//for num_in_face

	return 0;
}

int GetNodes(const int num_cur_cell, const Vector3& direction, std::vector<cell>& nodes_value, const std::vector<Type>& illum_old,
	const vector<Vector3>& directions, const vector<Type>& squares) {

	int num_cur_out_face = 0;  // 0 -- iый
	for (size_t num_node = 0; num_node < 3; ++num_node) {
		Vector3 x = nodes_value[num_cur_cell].x[num_node];

		CalculateNodeValue(num_cur_cell, num_node, direction, nodes_value, x, illum_old, directions, squares);
	}


	/*
	* Vector3 node;
	switch (num_cur_out_face)
	{

	case 1:// 1->2
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = 0;
			node[1] = straight_face.row(num_node)[0];
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выходящей грани

			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, illum_old, directions, squares);
		}
		//cout << "Number: " << num_cur_out_face << "\n";
		break;
	case 2://2->0
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = 0;
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выходящей грани				
			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, illum_old, directions, squares);
		}
		//cout << "Number: " << num_cur_out_face << "\n";
		break;
	case 0: //0->3
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = straight_face.row(num_node)[1];
			node[2] = 0;
			FromLocalToGlobalTetra(vertex_tetra, node, x);
			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, illum_old, directions, squares);
		}// x->координата узла на выходящей грани		}
	//	cout << "Number: " << num_cur_out_face << "\n";
		break;
	case 3: //3->1
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = inclined_face.row(num_node)[0];
			node[1] = inclined_face.row(num_node)[1];
			node[2] = 0;
			FromPlaneToTetra(inverse_transform_matrix, start_point_plane_coord, node, node);
			FromLocalToGlobalTetra(vertex_tetra, node, x);
			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, illum_old, directions, squares);
		}
		//	cout << "Number: " << num_cur_out_face << "\n";
		break;
	default:
		std::cout << "Number face is not {0,1,2,3}????\n";
		break;
	}
	*/

	int neighbor_id_face = nodes_value[num_cur_cell].neighbours_id_face[num_cur_out_face];
	nodes_value[neighbor_id_face / 4].nodes_value[neighbor_id_face % 4] =
		nodes_value[num_cur_cell].nodes_value[num_cur_out_face];
	return 0;
}

int GetNodes(const int num_cur_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell, const int num_cur_out_face,
	const Eigen::Matrix4d& vertex_tetra, const int* face_state, const Vector3& direction,
	std::vector<cell>& nodes_value, const std::vector<Type>& illum_old,
	const vector<Vector3>& directions, const vector<Type>& squares) {

	Vector3 x;
	Vector3 node;

	switch (num_cur_out_face)
	{

	case 1:// 1->2
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = 0;
			node[1] = straight_face.row(num_node)[0];
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выходящей грани

			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, illum_old, directions, squares);
		}
		//cout << "Number: " << num_cur_out_face << "\n";
		break;
	case 2://2->0
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = 0;
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выходящей грани				
			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, illum_old, directions, squares);
		}
		//cout << "Number: " << num_cur_out_face << "\n";
		break;
	case 0: //0->3
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = straight_face.row(num_node)[1];
			node[2] = 0;
			FromLocalToGlobalTetra(vertex_tetra, node, x);
			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, illum_old, directions, squares);
		}// x->координата узла на выходящей грани		}
	//	cout << "Number: " << num_cur_out_face << "\n";
		break;
	case 3: //3->1
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = inclined_face.row(num_node)[0];
			node[1] = inclined_face.row(num_node)[1];
			node[2] = 0;
			FromPlaneToTetra(inverse_transform_matrix, start_point_plane_coord, node, node);
			FromLocalToGlobalTetra(vertex_tetra, node, x);
			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, illum_old, directions, squares);
		}
		//	cout << "Number: " << num_cur_out_face << "\n";
		break;
	default:
		std::cout << "Number face is not {0,1,2,3}????\n";
		break;
	}

	int neighbor_id_face = nodes_value[num_cur_cell].neighbours_id_face[num_cur_out_face];
	nodes_value[neighbor_id_face / 4].nodes_value[neighbor_id_face % 4] =
		nodes_value[num_cur_cell].nodes_value[num_cur_out_face];
	return 0;
}

int CalculateNodeValue(const int num_cur_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell,
	const int num_cur_out_face, const int* face_state, const Vector3& direction,
	std::vector<cell>& nodes_value, const int num_node, const Eigen::Matrix4d& vertex_tetra, Vector3& x,
	const std::vector<Type>& illum_old, const vector<Vector3>& directions, const vector<Type>& squares) {

	Vector3 x0;

	for (size_t num_in_face = 0; num_in_face < 4; ++num_in_face) {
		if (!face_state[num_in_face]) continue;  // обрабатываем только входные грани


		IntersectionWithPlane(cur_cell->GetFace(num_in_face), x, direction, x0);

		if (InTriangle(num_cur_cell, unstructuredgrid, cur_cell, num_in_face, x0)) {

			Type s = (x - x0).norm();

			// значение на входящей грани
			Type I_x0 = CalculateIllumeOnInnerFace(num_cur_cell, num_in_face, vertex_tetra, x, x0, nodes_value);

			Type I = CurGetIllum(num_cur_cell, x0, s, I_x0, direction, illum_old, directions, squares);

			nodes_value[num_cur_cell].nodes_value[num_cur_out_face][num_node] = I;

			break;
		}

	}//for num_in_face

	return 0;
}

Type CalculateIllumeOnInnerFace(const int num_cell, const int num_in_face, const Eigen::Matrix4d& vertex_tetra,
	const Vector3& x, const Vector3& x0, const std::vector<cell>& nodes_value) {
	Type I_x0 = 0;

	if (nodes_value[num_cell].neighbours_id_face[num_in_face] == -1) {
		/*Граничные условия*/
		//I_x0 = BoundaryFunction(num_cell, x, direction, illum_old, directions, squares);
		return I_x0;
	}
	else {

		if (nodes_value[num_cell].nodes_value[num_in_face][0] < -600)
			cout << "Num_dir: " << num_cur_direction << " CalculateIllumeOnInnerFace:  Undefine cell / in face:" << num_cell << " / " << num_in_face << " !!!\n";

		Vector3 x0_local;

		FromGlobalToLocalTetra(vertex_tetra, x0, x0_local);
		Vector3 coef;// = GetInterpolationCoef(straight_face, nodes_value.find(global_num_in_face)->second);

		switch (num_in_face) {
		case 3:
			Vector3 local_plane_x0;
			FromTetraToPlane(transform_matrix, start_point_plane_coord, x0_local, local_plane_x0);

			coef = GetInterpolationCoef(inclined_face, nodes_value[num_cell].nodes_value[num_in_face]);
			I_x0 = local_plane_x0[0] * coef[0] + local_plane_x0[1] * coef[1] + coef[2];

			/*I_x0 = x0_local[1] * coef[0] + x0_local[2] * coef[1] + coef[2];
			Vector3 coef = GetInterpolationCoef(straight_face, nodes_value.find(global_num_in_face)->second);*/
			break;
		case 1:
			coef = GetInterpolationCoef(straight_face, nodes_value[num_cell].nodes_value[num_in_face]);
			I_x0 = x0_local[1] * coef[0] + x0_local[2] * coef[1] + coef[2];

			break;
		case 2:
			coef = GetInterpolationCoef(straight_face, nodes_value[num_cell].nodes_value[num_in_face]);
			I_x0 = x0_local[0] * coef[0] + x0_local[2] * coef[1] + coef[2];

			break;
		case 0:
			coef = GetInterpolationCoef(straight_face, nodes_value[num_cell].nodes_value[num_in_face]);
			I_x0 = x0_local[0] * coef[0] + x0_local[1] * coef[1] + coef[2];

			break;
		}
		if (I_x0 < 0) {
			count_negative_interpolation++;
			return 0;
		}

		return I_x0;
	}
}

Type CurGetIllum(const int cur_id, const int cur_direction, const Vector3 x, const Type s, const Type I_node_prev,
	   const vector<Type>& int_scattering) {
	// без интеграла рассеивания
		{
	/*		Type Ie = 10;
			Type k = 10;
			if (x.norm() > 0.25) { Ie = 0; k = 1; }

			Type I;
			if (s > 1e-10)
				I = Ie * (1 - exp(-s * k)) + I_node_prev * exp(-s * k);
			else
				I = Ie * (1 + s * k) - I_node_prev * s * k;

			if (I < 0)
				I = 0;
			return I;*/

		}

// test task
		{
			Type S = int_scattering[cur_direction * size_grid + cur_id]; //GetS(cur_id, cur_direction, illum_old, directions, squares);
			Type Ie = rad_en_loose_rate->GetTuple1(cur_id);
			Type alpha = absorp_coef->GetTuple1(cur_id);
			Type betta = 0;
			Type k = alpha + betta;

			Type I = exp(-k * s) * (I_node_prev * k + (exp(k * s) - 1) * (Ie * alpha + S * betta));
			I /= k;

			if (I < 0)
				I = 0;
			return I;
		}

		Type S = int_scattering[cur_direction * size_grid + cur_id]; //GetS(cur_id, cur_direction, illum_old, directions, squares);
		Type Ie = 10.;
		Type alpha = 5.;
		Type betta = 5.;
		Type k = alpha + betta;
		if (x.norm() > 0.3) {
			Ie = 0;
			alpha = 0.5;
			betta = 0.5;
			k = alpha + betta;
		}

		Type I = exp(-k * s) * (I_node_prev * k + (exp(k * s) - 1) * (Ie * alpha + S * betta));
		I /= k;

		if (I < 0)
			I = 0;
		return I;
}

Type GetValueInCenterCell(const int num_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell, const Vector3 center,
	const Vector3 direction,
	const Eigen::Matrix4d& vertex_tetra, const std::vector<cell>& nodes_value,
	const std::vector<Type>& illum_old, const vector<Vector3>& directions, const vector<Type>& squares) {
	/*Все грани должно быть определены*/
	Type value = -666;
	Vector3 x0;

	for (size_t i = 0; i < 4; i++) {

		IntersectionWithPlane(cur_cell->GetFace(i), center, direction, x0);
		if (InTriangle(num_cell, unstructuredgrid, cur_cell, i, x0)) {
			if ((center - x0).dot(direction) <= 0) continue;
			
			Type s = (center - x0).norm();
			Type I_x0 = CalculateIllumeOnInnerFace(num_cell, i, vertex_tetra, center, x0, nodes_value);

			//value = CurGetIllum(num_cell, x0, s, I_x0, direction, illum_old, directions, squares);
			break;
		}
	}
	if (value < 0)
		return 0;
	return value;
}
