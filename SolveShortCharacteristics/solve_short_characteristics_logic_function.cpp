#include "solve_short_characteristics_logic_function.h"

Type CalculateIllumeOnInnerFace(const int num_cell, const int num_in_face ,const Vector3& x, const std::vector<Vector2>&X0, 
	const std::vector<cell>& grid) {
	Type I_x0 = 0;
	
	if (grid[num_cell].neighbours_id_face[num_in_face] == -1) {
		/*√раничные услови€*/
		//I_x0 = BoundaryFunction(num_cell, x, direction, illum_old, directions, squares);
		return I_x0;
	}
	else if (grid[num_cell].neighbours_id_face[num_in_face] == -2) // внутренн€€ граница
	{
		
		id_try_pos++;

		Type data = res_inner_bound[pos_in_res++]; // защита на выход из диапазона??
		if (data >= 0) //данные от пересечени€ с диском или шаром
		{
		/*
			 ѕроверить пор€док. ƒанные в массиве должны лежать в соответствии с упор€доченным графом 
			 от направлени€ к направлению
		*/
			return 10; // I_x0;
			return data;  
		}
		else // определ€ющимм €вл€ютс€ противолежаащие грани (возможен расчет с учетом s=(x-x0).norm())
		{

			// результат не вполне пон€тен. ѕока лучше использовать константу или другие параметры области (шар == граница)

			//+dist_try_surface			
			int id = id_try_surface[id_try_pos - 1];  // будет лежать id грани			
			const int cell = id / 4;
			const int face = id % 4;		
			
			Vector3 coef = grid[cell].nodes_value[face];
			Vector2	x0_local = X0[posX0++];//grid[num_cell].x0_loc[num_in_face_dir];

			I_x0 = x0_local[0] * coef[0] + x0_local[1] * coef[1] + coef[2];
			
			//if (I_x0 < 0) I_x0 = 0;
			return 10; // I_x0;

		}
	}else {

		// debug условие
		if (grid[num_cell].nodes_value[num_in_face][0] < -600)
			cout << "Num_dir: " << num_cur_direction << " CalculateIllumeOnInnerFace:  Undefine cell / in face:" << num_cell << " / " << num_in_face << " !!!\n";


		Vector3 coef = grid[num_cell].nodes_value[num_in_face];
		Vector2	x0_local = X0[posX0++]; // grid[num_cell].x0_loc[num_in_face_dir];
		I_x0 = x0_local[0] * coef[0] + x0_local[1] * coef[1] + coef[2];
		

		if (I_x0 < 0) {
			count_negative_interpolation++;
			return 0;
		}

		return I_x0;
	}
}


Type CurGetIllum(const int cur_id, const int cur_direction, const Vector3 x, const Type s, const Type I_node_prev,
	   const vector<Type>& int_scattering) {
	// без интеграла рассеивани€
		{
			Type Ie = 10;
			Type k = 10;
			if ((x-Vector3(1,0,0)).norm() > 0.09) { Ie = 0; k = 1; }

			Type I;
			if (s > 1e-10)
				I = Ie * (1 - exp(-s * k)) + I_node_prev * exp(-s * k);
			else
				I = Ie * (1 + s * k) - I_node_prev * s * k;

			if (I < 0)
				I = 0;
			return I;

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
	/*¬се грани должно быть определены*/
	Type value = -666;
	Vector3 x0;

	for (size_t i = 0; i < 4; i++) {

		IntersectionWithPlane(cur_cell->GetFace(i), center, direction, x0);
		if (InTriangle(num_cell, unstructuredgrid, cur_cell, i, x0)) {
			if ((center - x0).dot(direction) <= 0) continue;
			
			Type s = (center - x0).norm();
		//	Type I_x0 = CalculateIllumeOnInnerFace(num_cell, i, vertex_tetra, center, x0, nodes_value);

			//value = CurGetIllum(num_cell, x0, s, I_x0, direction, illum_old, directions, squares);
			break;
		}
	}
	if (value < 0)
		return 0;
	return value;
}
