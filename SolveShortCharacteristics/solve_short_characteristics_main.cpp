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


size_t ReadFastGridData(const Type count_dir, Str_Type& name_file_graph,
	Str_Type& name_file_in_faces, Str_Type& name_file_out_faces, Str_Type& name_file_count_out_faces, Str_Type& name_file_local_x0, Str_Type& name_file_x,
	Str_Type& name_file_s, Str_Type& name_file_id_neighbors, Str_Type& name_file_centers,
	Str_Type& name_file_dist_try, Str_Type& name_file_id_try, Str_Type& name_file_res, std::vector<cell>& grid) {

	int N = 1542;//2233;

	grid.resize(N);

	std::unique_ptr<FILE, int(*)(FILE*)> file_id_neighbors(fopen((name_file_id_neighbors + ".bin").c_str(), "rb"), fclose);
	if (!file_id_neighbors) PRINTF("Error file_id_neighbors\n")
		for (int i = 0; i < N; i++) {
			grid[i].neighbours_id_face.resize(4, -5);
			fread_unlocked(grid[i].neighbours_id_face.data(), sizeof(int), 4, file_id_neighbors.get());
		}
	fclose(file_id_neighbors.get());


	// —формирован иначе. ћен€ть вывод в make struct
	std::unique_ptr<FILE, int(*)(FILE*)> file_count_out_id(fopen((name_file_count_out_faces + ".bin").c_str(), "rb"), fclose);
	if (!file_count_out_id) PRINTF("Error file_count_out_id\n")
		for (int i = 0; i < N; i++) {
			grid[i].out_count.resize(count_dir);
			fread_unlocked(grid[i].out_count.data(), sizeof(int), count_dir, file_count_out_id.get());
		}
	fclose(file_count_out_id.get());


	std::vector<int> size_data_cell(N, 0);
	for (int i = 0; i < N; i++)
		for (int num_dir = 0; num_dir < count_dir; num_dir++)
			size_data_cell[i] += grid[i].out_count[num_dir];


	for (int num_cell = 0; num_cell < N; ++num_cell)
	{
		grid[num_cell].in_id.resize(3 * size_data_cell[num_cell]);
		grid[num_cell].out_id.resize(size_data_cell[num_cell]);
		grid[num_cell].s.resize(3 * size_data_cell[num_cell]);
		grid[num_cell].x.resize(3 * size_data_cell[num_cell], Vector3(0, 0, 0));
		grid[num_cell].x0_loc.resize(3 * size_data_cell[num_cell], Vector2(0, 0));
	}

	std::vector<int> start_dir1(N, 0);
	std::vector<int> start_dir3(N, 0);
	std::vector<int> sorted_id_cell(N);

	for (int num_dir = 0; num_dir < count_dir; num_dir++)
	{

		std::unique_ptr<FILE, int(*)(FILE*)> file_out_id(fopen((name_file_out_faces + to_string(num_dir) + ".bin").c_str(), "rb"), fclose);
		if (!file_out_id) PRINTF("Error file_out_id\n")
			std::unique_ptr<FILE, int(*)(FILE*)> file_in_id(fopen((name_file_in_faces + to_string(num_dir) + ".bin").c_str(), "rb"), fclose);
		if (!file_in_id) PRINTF("Error file_in_id\n")
			std::unique_ptr<FILE, int(*)(FILE*)> file_x(fopen((name_file_x + to_string(num_dir) + ".bin").c_str(), "rb"), fclose);
		if (!file_x) PRINTF("Error file_x\n")
			std::unique_ptr<FILE, int(*)(FILE*)> file_x0_local(fopen((name_file_local_x0 + to_string(num_dir) + ".bin").c_str(), "rb"), fclose);
		if (!file_x0_local) PRINTF("Error file_x0_local\n")
			std::unique_ptr<FILE, int(*)(FILE*)> file_s(fopen((name_file_s + to_string(num_dir) + ".bin").c_str(), "rb"), fclose);
		if (!file_s) PRINTF("Error file_s\n")

			if (ReadGraphBin(name_file_graph + to_string(num_dir) + ".bin", sorted_id_cell)) PRINTF("Graphd%d not open\n", num_dir)

				int num_cell;
		for (int h = 0; h < N; ++h)
		{
			num_cell = sorted_id_cell[h];
			for (int i = 0; i < grid[num_cell].out_count[num_dir]; i++)
			{
				fread_unlocked(&grid[num_cell].out_id[start_dir1[num_cell]], sizeof(int), 1, file_out_id.get());

				fread_unlocked(grid[num_cell].in_id.data() + start_dir3[num_cell], sizeof(int), 3, file_in_id.get());

				fread_unlocked(grid[num_cell].s.data() + start_dir3[num_cell], sizeof(double), 3, file_s.get());

				fread_unlocked(grid[num_cell].x[start_dir3[num_cell]].data(), sizeof(double), 3, file_x.get());
				fread_unlocked(grid[num_cell].x[start_dir3[num_cell] + 1].data(), sizeof(double), 3, file_x.get());
				fread_unlocked(grid[num_cell].x[start_dir3[num_cell] + 2].data(), sizeof(double), 3, file_x.get());

				fread_unlocked(grid[num_cell].x0_loc[start_dir3[num_cell]].data(), sizeof(double), 2, file_x0_local.get());
				fread_unlocked(grid[num_cell].x0_loc[start_dir3[num_cell] + 1].data(), sizeof(double), 2, file_x0_local.get());
				fread_unlocked(grid[num_cell].x0_loc[start_dir3[num_cell] + 2].data(), sizeof(double), 2, file_x0_local.get());

				start_dir3[num_cell] += 3;
				start_dir1[num_cell]++;

			}// for дл€ i_го узла
		}// for по €чейкам

		fclose(file_out_id.get());
		fclose(file_in_id.get());
		fclose(file_x.get());
		fclose(file_x0_local.get());
		fclose(file_s.get());
	} // for по направлени€м


	int n;
	std::unique_ptr<FILE, int(*)(FILE*)> file_res(fopen((name_file_res + ".bin").c_str(), "rb"), fclose);
	if (!file_res) PRINTF("Error file_res\n")

	fread_unlocked(&n, sizeof(int), 1, file_res.get());
	res_inner_bound.resize(n);
	fread_unlocked(res_inner_bound.data(), sizeof(Type), n, file_res.get());

	fclose(file_res.get());


	std::unique_ptr<FILE, int(*)(FILE*)> file_id(fopen((name_file_id_try + ".bin").c_str(), "rb"), fclose);
	if (!file_id) PRINTF("Error file_id\n")

	fread_unlocked(&n, sizeof(int), 1, file_id.get());
	id_try_surface.resize(n);
	fread_unlocked(id_try_surface.data(), sizeof(int), n, file_id.get());

	fclose(file_id.get());

	return 0;
}


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
	std::string out_file_E1d;
	
	std::string name_file_graph;
	std::string name_file_in_faces;
	std::string name_file_out_faces;
	std::string name_file_count_out_faces;
	std::string name_file_local_x0;
	std::string name_file_x;
	std::string name_file_s;

	std::string name_file_id_neighbors;
	std::string name_file_centers;
	
	std::string name_file_dist_try;
	std::string name_file_id_try;
	std::string name_file_res;	

	if (ReadStartSettings(name_file_settings, class_file_vtk, name_file_vtk, name_file_sphere_direction, out_file_grid_vtk, out_file_E1d, 
		name_file_graph, name_file_in_faces, name_file_out_faces, name_file_count_out_faces, name_file_local_x0, name_file_x, name_file_s,
		name_file_id_neighbors, name_file_centers, name_file_dist_try, name_file_id_try, name_file_res)) {
		std::cout << "Error reading the start settings\n";
		return 1;
	}


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

	std::vector<cell> grid;

	std::vector<cell> fast_grid;
	
	const int count_directions =  directions.size();
	const int count_cells = size_grid;

	_clock = -omp_get_wtime();
	if (ReadFastGridData(count_directions, name_file_graph, name_file_in_faces, name_file_out_faces, name_file_count_out_faces,
		name_file_local_x0, name_file_x, name_file_s, name_file_id_neighbors, name_file_centers,
		name_file_dist_try, name_file_id_try, name_file_res, grid)) {
		std::cout << "Error fast reading the data grid vtk\n";
		return 1;
	}
	_clock += omp_get_wtime();
	std::cout << "\n Reading time of the data_grid file: " << _clock << "\n";

		
	InitGlobalValue(start_point_plane_coord, transform_matrix, inverse_transform_matrix, straight_face, inclined_face);

	// ”пор€доченные индексы €чеек по данному направлению
	vector<IntId> sorted_id_cell(count_cells);


	//Eigen::Matrix4d vertex_tetra;
	/* x1 x2 x3 x4
	*  y1 y2 y3 y4
	*  z1 z2 z3 z4
	*  1  1  1  1
	*/

	std::vector<Type> Illum(count_cells * directions.size() , 0);
	std::vector<Type> Illum2(count_cells * directions.size(), 0);


	/*std::vector<Normals> normals;
	ReadNormalFile(name_file_normals, normals);*/


	int count = 0;
	Type norm = 0;
	Vector3 direction;

	ofstream ofile;
	ofile.open("File_with_Logs.txt");

	vector<IntId> full_sorted_id_cell(count_cells*count_directions);

	for (size_t i = 0; i < count_directions; i++)
	{
		ReadGraphBin(name_file_graph + to_string(i) + ".bin", sorted_id_cell);
		for (size_t j = 0; j < count_cells; j++)
		{
			full_sorted_id_cell[i*count_cells + j] = sorted_id_cell[j];
		}	
	}
	
	do {
		Type _clock = -omp_get_wtime();
		{

			/*---------------------------------- далее FOR по направлени€м----------------------------------*/
			for (int num_direction = 0; num_direction < count_directions; ++num_direction)
			{
				num_cur_direction = num_direction;
				direction = directions[num_direction];

				//ReadGraph(name_file_graph + to_string(num_direction) + ".txt", sorted_id_cell);
				//ResetNodesValue(grid);

				int num_cell;

				/*---------------------------------- далее FOR по €чейкам----------------------------------*/
				for (int h = 0; h < count_cells; ++h) {
					num_cell = full_sorted_id_cell[num_direction * count_cells + h]; //sorted_id_cell[h];

					int n_out = grid[num_cell].out_count[num_direction];

					if (num_direction == 0) {
						grid[num_cell].start_dir1 = 0;
						grid[num_cell].start_dir3 = 0;
					}

					int start_1 = grid[num_cell].start_dir1;
					int start_3 = grid[num_cell].start_dir3;
					Vector3 I;

					for (int i = 0; i < n_out; ++i) {
						int num_out_face = grid[num_cell].out_id[start_1 + i];
						
//GetNodes
						
						for (size_t num_node = 0; num_node < 3; ++num_node) {
							Vector3 x = grid[num_cell].x[start_3 + i + num_node];

//CalculateNodeValue:
							int num_in_face = grid[num_cell].in_id[start_3 + i + num_node];
							Type I_x0 = CalculateIllumeOnInnerFace(num_cell, num_in_face, start_3 + i + num_node, x, grid);

							Type s = grid[num_cell].s[start_3 + i+ num_node];

							I[num_node] = CurGetIllum(num_cell, x, s, I_x0, direction, Illum2, directions, squares);
						}//num_node

						Vector3 coef;
						if (num_out_face == 3)
							coef = GetInterpolationCoefInverse(inclined_face_inverse, I);
						else
							coef = GetInterpolationCoefInverse(straight_face_inverse, I);

						grid[num_cell].nodes_value[num_out_face] = coef;

						int neighbor_id_face = grid[num_cell].neighbours_id_face[num_out_face];
						if(neighbor_id_face>=0)
						grid[neighbor_id_face / 4].nodes_value[neighbor_id_face % 4] = coef;


					} //num_out_face
							

					Illum[num_direction * count_cells + num_cell] = (I[0] + I[1] + I[2]) / 3;

					grid[num_cell].start_dir1 = start_1 + n_out;
					grid[num_cell].start_dir3 = start_3 + 3 * n_out;

				}
				/*---------------------------------- конец FOR по €чейкам----------------------------------*/

				//std::cout << "End direction number: " << num_direction << '\n';
			}
			/*---------------------------------- конец FOR по направлени€м----------------------------------*/
		}
		Illum.swap(Illum2);
		
		_clock += omp_get_wtime();
		count++;
		norm = NormIllum(Illum, Illum2);
		std::cout << "Error:= " << norm << '\n';
		std::cout << "Time of iter: " << _clock << '\n';
		std::cout << "End iter_count number: " << count << '\n';

		ofile << "Error:= " << norm << '\n';
		ofile << "Time of iter: " << _clock << '\n';
		ofile << "End iter_count number: " << count << '\n';
		

	} while (norm > accuracy && count < max_number_of_iter);


	Illum.swap(Illum2);
	ofile.close();

	vector<Type> energy(count_cells);
	MakeEnergy(Illum, squares, square_surface, energy);

	WriteFileSolution(out_file_grid_vtk, Illum, energy, name_file_vtk);
	return 0;
}


