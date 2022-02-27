#pragma once
#ifndef SHORT_CHARACTERISTICS_MAIN_H
#define SHORT_CHARACTERISTICS_MAIN_H

#include "solve_short_characteristics_headers.h"
#include "solve_short_characteristics_global_structure.h"
#include "solve_short_characteristics_calculations.h"
#include "solve_short_characteristics_logic_function.h"

template<typename Type>
size_t ReadStartSettings(std::string name_file_settings, Type& class_file_vtk, std::string& name_file_vtk,
	std::string& name_file_sphere_direction, std::string& out_file_grid_vtk,  std::string& out_file_E1d, std::string& name_file_graph,
	std::string& name_file_in_faces, std::string& name_file_out_faces, std::string& name_file_count_out_faces, std::string& name_file_local_x0, std::string& name_file_x,
	std::string& name_file_s, std::string& name_file_id_neighbors, std::string& name_file_centers,
	std::string& name_file_dist_try, std::string& name_file_id_try, std::string& name_file_res) {

	std::ifstream ifile;
	ifile.open(name_file_settings);
	if (!ifile.is_open()) {
		std::cerr << " Error : file settings SolveShortCharacteristics is not open !\n";
		return 1;
	}

	std::string str; // переменна€ дл€ перевода строки при чтении из файла

	ifile >> class_file_vtk;
	getline(ifile, str);
	getline(ifile, name_file_vtk);
	getline(ifile, name_file_sphere_direction);
	getline(ifile, out_file_grid_vtk);
	getline(ifile, out_file_E1d);
	getline(ifile, name_file_graph);
	
	getline(ifile, name_file_in_faces);
	getline(ifile, name_file_out_faces); 
	getline(ifile, name_file_count_out_faces);
	getline(ifile, name_file_local_x0);
	getline(ifile, name_file_x);
	getline(ifile, name_file_s);

	getline(ifile, name_file_id_neighbors);
	getline(ifile, name_file_centers);

	getline(ifile, name_file_dist_try);
	getline(ifile, name_file_id_try);
	getline(ifile, name_file_res);	

	ifile.close();
	return 0;
}


template<typename Str_Type, typename T>
size_t ReadSimpleFile(const Str_Type name_file, std::vector<T>& data_array) {
	// ‘айл должен содержать в первой строке число элементов. ƒалее последовательные данные
	
	std::ifstream ifile;
	ifile.open(name_file);
	if (!ifile.is_open()) {
		std::cerr << " Error : file <<Simple>> is not open !\n";
		return 1;
	}
	
	int N;
	ifile >> N;

	data_array.resize(N, 0);

	for (int i = 0; i < N; i++)
	{
		ifile >> data_array[i];
	}

	ifile.close();
	return 0;
}

template<typename Type, typename Str_Type>
size_t ReadGridData(const Type count_dir, Str_Type& name_file_graph,
	Str_Type& name_file_in_faces, Str_Type& name_file_out_faces, Str_Type& name_file_count_out_faces, Str_Type& name_file_local_x0, Str_Type& name_file_x,
	Str_Type& name_file_s, Str_Type& name_file_id_neighbors, Str_Type& name_file_centers, std::vector<cell>& grid) {

	std::ifstream ifile_graph;
	std::ifstream ifile_out_id;
	std::ifstream ifile_count_out_id;
	std::ifstream ifile_in_id;
	std::ifstream ifile_x;
	std::ifstream ifile_x0_local;
	std::ifstream ifile_s;	

	std::ifstream ifile_id_neighbors(name_file_id_neighbors);
	std::ifstream ifile_centers(name_file_centers);
	if (!ifile_id_neighbors.is_open()) { printf("File id neighbors not open\n"); return 1; }
	if (!ifile_centers.is_open()){ printf("File centers not open\n"); return 1; }

	int N;
	ifile_id_neighbors >> N;
	N /= 4;

	grid.resize(N);

	for (int i = 0; i < N; i++) {
		grid[i].neighbours_id_face.resize(4, -5);
		for (int j = 0; j < 4; j++)
			ifile_id_neighbors >> grid[i].neighbours_id_face[j];
	}

	ifile_id_neighbors.close();
	ifile_centers.close();

	vector<IntId> sorted_id_cell(N);
	

	for (int i = 0; i < N; i++)
		grid[i].out_count.resize(count_dir);

	std::vector<int> size_data_cell(N,0);

	for (int num_dir = 0; num_dir < count_dir; num_dir++)
	{
		ReadGraph(name_file_graph + to_string(num_dir) + ".txt", sorted_id_cell);

		ifile_count_out_id.open(name_file_count_out_faces + std::to_string(num_dir) + ".txt");
		if (!ifile_count_out_id.is_open()) { printf("CountOutid%d not open\n", num_dir); return 1; }

		ifile_count_out_id >> N;
		
		int buf;
		int h;
		for (int i = 0; i < N; i++) {
			ifile_count_out_id >> buf;
			h = sorted_id_cell[i];
			grid[h].out_count[num_dir] = buf;
			size_data_cell[h] += buf;
		}

		ifile_count_out_id.close();
	}

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
	for (int num_dir = 0; num_dir < count_dir; num_dir++)
	{

		std::ifstream ifile_out_id(name_file_out_faces + std::to_string(num_dir) + ".txt");
		std::ifstream ifile_in_id(name_file_in_faces + std::to_string(num_dir) + ".txt");
		std::ifstream ifile_x(name_file_x + std::to_string(num_dir) + ".txt");
		std::ifstream ifile_x0_local(name_file_local_x0 + std::to_string(num_dir) + ".txt");
		std::ifstream ifile_s(name_file_s + std::to_string(num_dir) + ".txt");
		
		if(ReadGraph(name_file_graph + to_string(num_dir) + ".txt", sorted_id_cell)) { printf("Graphd%d not open\n", num_dir); return 1; }
		if (!ifile_out_id.is_open()) { printf("Outid%d not open\n", num_dir); return 1; }
		if (!ifile_in_id.is_open()){ printf("Inid%d not open\n", num_dir); return 1; }
		if (!ifile_x.is_open()){ printf("X%d not open\n", num_dir); return 1;}
		if (!ifile_x0_local.is_open()){ printf("file X%d not open\n", num_dir); return 1; }
		if (!ifile_s.is_open()){ printf("S%d not open\n", num_dir); return 1; }
		

		ifile_out_id >> N;		
		ifile_in_id >> N;
		ifile_x >> N;
		ifile_x0_local >> N;
		ifile_s >> N;

		int num_cell;
		for (int h = 0; h < N; ++h)
		{
			num_cell = sorted_id_cell[h];
			for (int i = 0; i < grid[num_cell].out_count[num_dir]; i++)
			{
				ifile_out_id >> grid[num_cell].out_id[start_dir1[num_cell]];

				ifile_in_id >> grid[num_cell].in_id[start_dir3[num_cell]] 
							>> grid[num_cell].in_id[start_dir3[num_cell] + 1] 
							>> grid[num_cell].in_id[start_dir3[num_cell] + 2];

				ifile_s >> grid[num_cell].s[start_dir3[num_cell]]
						>> grid[num_cell].s[start_dir3[num_cell] + 1]
						>> grid[num_cell].s[start_dir3[num_cell] + 2];

				ifile_x >> grid[num_cell].x[start_dir3[num_cell]](0) >> grid[num_cell].x[start_dir3[num_cell]](1) >> grid[num_cell].x[start_dir3[num_cell]](2)
						>> grid[num_cell].x[start_dir3[num_cell] + 1](0) >> grid[num_cell].x[start_dir3[num_cell] + 1](1) >> grid[num_cell].x[start_dir3[num_cell] + 1](2)
						>> grid[num_cell].x[start_dir3[num_cell] + 2](0) >> grid[num_cell].x[start_dir3[num_cell] + 2](1) >> grid[num_cell].x[start_dir3[num_cell] + 2](2);

				ifile_x0_local >> grid[num_cell].x0_loc[start_dir3[num_cell]](0) >> grid[num_cell].x0_loc[start_dir3[num_cell]](1)
							   >> grid[num_cell].x0_loc[start_dir3[num_cell] + 1](0) >> grid[num_cell].x0_loc[start_dir3[num_cell] + 1](1)
							   >> grid[num_cell].x0_loc[start_dir3[num_cell] + 2](0) >> grid[num_cell].x0_loc[start_dir3[num_cell] + 2](1);

							start_dir3[num_cell] += 3;
							start_dir1[num_cell]++;
					
			}// for дл€ i_го узла
		}// for по €чейкам

		
		ifile_out_id.close();
		ifile_in_id.close();
		ifile_x.close();
		ifile_x0_local.close();
		ifile_s.close();
	} // for по направлени€м


	return 0;
}

#include <chrono>
#include <cstdio>
#include<inttypes.h>
#include <memory>

#ifdef _MSC_VER
#define fwrite_unlocked_fwrite_nolock
#endif

using namespace std::chrono;

template<typename Type, typename Str_Type>
size_t NoReadFastGridData(const Type count_dir, Str_Type& name_file_graph,
	Str_Type& name_file_in_faces, Str_Type& name_file_out_faces, Str_Type& name_file_count_out_faces, Str_Type& name_file_local_x0, Str_Type& name_file_x,
	Str_Type& name_file_s, Str_Type& name_file_id_neighbors, Str_Type& name_file_centers, std::vector<cell>& grid) {

	std::ifstream ifile_graph;
	std::ifstream ifile_out_id;
	std::ifstream ifile_count_out_id;
	std::ifstream ifile_in_id;
	std::ifstream ifile_x;
	std::ifstream ifile_x0_local;
	std::ifstream ifile_s;

	std::ifstream ifile_id_neighbors(name_file_id_neighbors);
	std::ifstream ifile_centers(name_file_centers);
	if (!ifile_id_neighbors.is_open()) { printf("File id neighbors not open\n"); return 1; }
	if (!ifile_centers.is_open()) { printf("File centers not open\n"); return 1; }

	int N;
	ifile_id_neighbors >> N;
	N /= 4;

	grid.resize(N);

	for (int i = 0; i < N; i++) {
		grid[i].neighbours_id_face.resize(4, -5);
		for (int j = 0; j < 4; j++)
			ifile_id_neighbors >> grid[i].neighbours_id_face[j];
	}

	ifile_id_neighbors.close();
	ifile_centers.close();

	vector<IntId> sorted_id_cell(N);


	for (int i = 0; i < N; i++)
		grid[i].out_count.resize(count_dir);

	std::vector<int> size_data_cell(N, 0);

	for (int num_dir = 0; num_dir < count_dir; num_dir++)
	{
		ReadGraph(name_file_graph + to_string(num_dir) + ".txt", sorted_id_cell);

		ifile_count_out_id.open(name_file_count_out_faces + std::to_string(num_dir) + ".txt");
		if (!ifile_count_out_id.is_open()) { printf("CountOutid%d not open\n", num_dir); return 1; }

		ifile_count_out_id >> N;

		int buf;
		int h;
		for (int i = 0; i < N; i++) {
			ifile_count_out_id >> buf;
			h = sorted_id_cell[i];
			grid[h].out_count[num_dir] = buf;
			size_data_cell[h] += buf;
		}

		ifile_count_out_id.close();
	}

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
	for (int num_dir = 0; num_dir < count_dir; num_dir++)
	{

		std::ifstream ifile_out_id(name_file_out_faces + std::to_string(num_dir) + ".txt");
		std::ifstream ifile_in_id(name_file_in_faces + std::to_string(num_dir) + ".txt");
		std::ifstream ifile_x(name_file_x + std::to_string(num_dir) + ".txt");
		std::ifstream ifile_x0_local(name_file_local_x0 + std::to_string(num_dir) + ".txt");
		std::ifstream ifile_s(name_file_s + std::to_string(num_dir) + ".txt");

		if (ReadGraph(name_file_graph + to_string(num_dir) + ".txt", sorted_id_cell)) { printf("Graphd%d not open\n", num_dir); return 1; }
		if (!ifile_out_id.is_open()) { printf("Outid%d not open\n", num_dir); return 1; }
		if (!ifile_in_id.is_open()) { printf("Inid%d not open\n", num_dir); return 1; }
		if (!ifile_x.is_open()) { printf("X%d not open\n", num_dir); return 1; }
		if (!ifile_x0_local.is_open()) { printf("file X%d not open\n", num_dir); return 1; }
		if (!ifile_s.is_open()) { printf("S%d not open\n", num_dir); return 1; }


		ifile_out_id >> N;
		ifile_in_id >> N;
		ifile_x >> N;
		ifile_x0_local >> N;
		ifile_s >> N;

		int num_cell;
		for (int h = 0; h < N; ++h)
		{
			num_cell = sorted_id_cell[h];
			for (int i = 0; i < grid[num_cell].out_count[num_dir]; i++)
			{
				ifile_out_id >> grid[num_cell].out_id[start_dir1[num_cell]];

				ifile_in_id >> grid[num_cell].in_id[start_dir3[num_cell]]
					>> grid[num_cell].in_id[start_dir3[num_cell] + 1]
					>> grid[num_cell].in_id[start_dir3[num_cell] + 2];

				ifile_s >> grid[num_cell].s[start_dir3[num_cell]]
					>> grid[num_cell].s[start_dir3[num_cell] + 1]
					>> grid[num_cell].s[start_dir3[num_cell] + 2];

				ifile_x >> grid[num_cell].x[start_dir3[num_cell]](0) >> grid[num_cell].x[start_dir3[num_cell]](1) >> grid[num_cell].x[start_dir3[num_cell]](2)
					>> grid[num_cell].x[start_dir3[num_cell] + 1](0) >> grid[num_cell].x[start_dir3[num_cell] + 1](1) >> grid[num_cell].x[start_dir3[num_cell] + 1](2)
					>> grid[num_cell].x[start_dir3[num_cell] + 2](0) >> grid[num_cell].x[start_dir3[num_cell] + 2](1) >> grid[num_cell].x[start_dir3[num_cell] + 2](2);

				ifile_x0_local >> grid[num_cell].x0_loc[start_dir3[num_cell]](0) >> grid[num_cell].x0_loc[start_dir3[num_cell]](1)
					>> grid[num_cell].x0_loc[start_dir3[num_cell] + 1](0) >> grid[num_cell].x0_loc[start_dir3[num_cell] + 1](1)
					>> grid[num_cell].x0_loc[start_dir3[num_cell] + 2](0) >> grid[num_cell].x0_loc[start_dir3[num_cell] + 2](1);

				start_dir3[num_cell] += 3;
				start_dir1[num_cell]++;

			}// for дл€ i_го узла
		}// for по €чейкам


		ifile_out_id.close();
		ifile_in_id.close();
		ifile_x.close();
		ifile_x0_local.close();
		ifile_s.close();
	} // for по направлени€м


	return 0;
}
#endif