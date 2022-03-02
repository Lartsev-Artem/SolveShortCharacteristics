#include "solve_short_characteristics_main.h"
Vector3 start_point_plane_coord;   // начало координат плоскости
Matrix3 transform_matrix;          // матрица перехода из базового тетраэдра в плоскость
Matrix3 inverse_transform_matrix;  // матрица перехода из плоскости в базовый тетраэдр

Matrix3	straight_face;  // 3 узла интерпол€ции
Matrix3 inclined_face;  // 3 узла интерпол€ции на наклонной плоскости

Matrix3	straight_face_inverse;  // 3 узла интерпол€ции
Matrix3 inclined_face_inverse;  // 3 узла интерпол€ции на наклонной плоскости

//std::vector<Type> Illum2;

// скал€рные данные сетки (unstructured_grid)
vtkDataArray* density;
vtkDataArray* absorp_coef;
vtkDataArray* rad_en_loose_rate;

Type square_surface;  // площадь поверхности дискретной 

Vector3 center_local_sphere;  // центр описанной сферы около стандартного тетраэдра

int num_cur_direction; // номер текущего направлени€=

int count_negative_interpolation; // число отрицательных значений интерпол€ции интесивности

int size_grid;
std::vector<Type> res_inner_bound;  // значение на внутренней границе
std::vector<int> id_try_surface;    // id определ€ющих €чеек внутренней границе


int posX;
int posX0;
int posOutC;
int posOut;
int posIn;
int posS;
int pos_in_res;
int id_try_pos;

int main(int argc, char* argv[])
{
	std::string name_file_settings = "";
	int max_number_of_iter = 1;
	Type accuracy = 1e-5;

	if (argc <= 1)
		name_file_settings = "D:\\Desktop\\FilesCourse\\settings_file_solve.txt";
	else
		name_file_settings = argv[1];
	if (argc > 2)
		max_number_of_iter = std::stoi(argv[2]);

	std::cout << "Max_number_of_iter= " << max_number_of_iter << '\n';

	size_t class_file_vtk;
	std::string name_file_vtk;
	std::string name_file_sphere_direction;
	std::string out_file_grid_vtk;

	std:string main_dir;
	
	/*if (ReadStartSettings(name_file_settings, class_file_vtk, name_file_vtk, name_file_sphere_direction, out_file_grid_vtk,
		name_file_graph, name_file_in_faces, name_file_out_faces, name_file_count_out_faces, name_file_local_x0, name_file_x, name_file_s,
		name_file_id_neighbors, name_file_centers, name_file_dist_try, name_file_id_try, name_file_res, name_file_sizes)) {
		std::cout << "Error reading the start settings\n";
		return 1;
	}*/
	if (ReadStartSettings(name_file_settings, class_file_vtk, name_file_vtk, name_file_sphere_direction, out_file_grid_vtk,
		main_dir)) {
		std::cout << "Error reading the start settings\n";
		return 1;
	}
	std::string name_file_graph = main_dir + "graph";
	std::string name_file_in_faces = main_dir + "InId";
	std::string name_file_out_faces = main_dir + "OutId";
	std::string name_file_count_out_faces = main_dir + "CountOutId";
	std::string name_file_local_x0 = main_dir + "locX0";
	std::string name_file_x = main_dir + "X";
	std::string name_file_s = main_dir + "S";

	std::string name_file_id_neighbors = main_dir + "id_neighbors";
	std::string name_file_centers = main_dir + "centers";

	std::string name_file_dist_try = main_dir + "dist_defining_faces";
	std::string name_file_id_try = main_dir + "id_defining_faces";
	std::string name_file_res = main_dir + "ResBound";
	std::string name_file_sizes = main_dir + "Size";


	Type _clock = -omp_get_wtime();
	if (ReadDataArray(class_file_vtk, name_file_vtk, size_grid, density, absorp_coef, rad_en_loose_rate, true)) {
		std::cout << "Error reading the data array vtk\n";
		return 1;
	}
	_clock += omp_get_wtime();
	std::cout << "\n Reading time of the vtk_grid file: " << _clock << "\n";

	vector<Vector3> directions;
	vector<Type> squares;

	_clock = -omp_get_wtime();
	ReadSphereDirectionDecartToSpherical(name_file_sphere_direction, directions, squares, square_surface);
	_clock += omp_get_wtime();
	std::cout << "\n Reading time of the sphere_direction file: " << _clock << "\n";


	const int count_directions = directions.size();
	const int count_cells = size_grid;

	std::vector<cell> grid;
	std::vector<int> OutC;
	std::vector<int> Out;
	std::vector<int> In;
	std::vector<Type>S;
	std::vector<Vector3> X;
	std::vector<Vector2> X0;

	vector<IntId> sorted_id_cell(count_cells * count_directions); 	// ”пор€доченные индексы €чеек по данному направлению

	_clock = -omp_get_wtime();

	if (ReadCompactFastGridData(count_directions, size_grid, name_file_in_faces, name_file_out_faces, name_file_count_out_faces,
		name_file_local_x0, name_file_x, name_file_s, name_file_id_neighbors, name_file_centers,
		name_file_dist_try, name_file_id_try, name_file_res, name_file_sizes, name_file_graph, grid, OutC, Out, In, S, X, X0,
		res_inner_bound, id_try_surface, sorted_id_cell)) {
		std::cout << "Error fast reading the data grid vtk\n";
		return 1;
	}
	_clock += omp_get_wtime();
	std::cout << "\n Reading time of the data_grid file: " << _clock << "\n";


	InitGlobalValue(start_point_plane_coord, transform_matrix, inverse_transform_matrix, straight_face, inclined_face);

	//Eigen::Matrix4d vertex_tetra;
	/* x1 x2 x3 x4
	*  y1 y2 y3 y4
	*  z1 z2 z3 z4
	*  1  1  1  1
	*/

	std::vector<Type> Illum(count_cells * directions.size(), 0);
	std::vector<Type> Illum2(count_cells * directions.size(), 0);

	std::vector<Type> int_scattering(count_cells * directions.size(), 0);


	int count = 0;
	Type norm = 0;
	Vector3 direction;

	ofstream ofile;
	ofile.open("File_with_Logs.txt");

	do {
		Type _clock = -omp_get_wtime();
		{

			posX = 0;
			posX0 = 0;
			posOutC = 0;
			posOut = 0;
			posIn = 0;
			posS = 0;

			pos_in_res = 0;
			id_try_pos = 0;


			if (count) CalculateInt(count_cells, count_directions, Illum2, directions, squares, int_scattering);

			/*---------------------------------- далее FOR по направлени€м----------------------------------*/
			for (int num_direction = 0; num_direction < count_directions; ++num_direction)
			{
				num_cur_direction = num_direction;
				direction = directions[num_direction];

				int num_cell;

				/*---------------------------------- далее FOR по €чейкам----------------------------------*/
				for (int h = 0; h < count_cells; ++h) {
					num_cell = sorted_id_cell[num_direction * count_cells + h]; //sorted_id_cell[h];

					int n_out = OutC[posOutC++]; //grid[num_cell].out_count[num_direction];


					Vector3 I;

					for (int i = 0; i < n_out; ++i) {
						int num_out_face = Out[posOut++];//grid[num_cell].out_id[start_1 + i];

						//GetNodes

						for (size_t num_node = 0; num_node < 3; ++num_node) {
							Vector3 x = X[posX++];//grid[num_cell].x[start_3 + i + num_node];

							//CalculateNodeValue:
							int num_in_face = In[posIn++]; //grid[num_cell].in_id[start_3 + i + num_node];
							Type I_x0 = CalculateIllumeOnInnerFace(num_cell, num_in_face, x, X0, grid);

							Type s = S[posS++]; //grid[num_cell].s[start_3 + i + num_node];

							I[num_node] = CurGetIllum(num_cell, num_direction, x, s, I_x0, int_scattering);
						}//num_node

						Vector3 coef;
						if (num_out_face == 3)
							coef = GetInterpolationCoefInverse(inclined_face_inverse, I);
						else
							coef = GetInterpolationCoefInverse(straight_face_inverse, I);

						grid[num_cell].nodes_value[num_out_face] = coef;

						int neighbor_id_face = grid[num_cell].neighbours_id_face[num_out_face];
						if (neighbor_id_face >= 0)
							grid[neighbor_id_face / 4].nodes_value[neighbor_id_face % 4] = coef;


					} //num_out_face


					Illum[num_direction * count_cells + num_cell] = (I[0] + I[1] + I[2]) / 3;

				}
				/*---------------------------------- конец FOR по €чейкам----------------------------------*/

				//std::cout << "End direction number: " << num_direction << '\n';
			}
			/*---------------------------------- конец FOR по направлени€м----------------------------------*/
		}
		Illum.swap(Illum2);

		_clock += omp_get_wtime();

		if (count % 2 == 0) norm = NormIllumOmp(Illum, Illum2);

		std::cout << "Error:= " << norm << '\n';
		std::cout << "Time of iter: " << _clock << '\n';
		std::cout << "End iter_count number: " << count << '\n';

		ofile << "Error:= " << norm << '\n';
		ofile << "Time of iter: " << _clock << '\n';
		ofile << "End iter_count number: " << count << '\n';

		count++;
	} while (norm > accuracy && count < max_number_of_iter);


	Illum.swap(Illum2);
	ofile.close();

	vector<Type> energy(count_cells);
	MakeEnergy(Illum, squares, square_surface, energy);

	WriteFileSolution(out_file_grid_vtk, Illum, energy, name_file_vtk);
	return 0;
}


