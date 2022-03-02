#pragma once
#ifndef SHORT_CHARACTERISTICS_MAIN_H
#define SHORT_CHARACTERISTICS_MAIN_H

#include "solve_short_characteristics_headers.h"
#include "solve_short_characteristics_global_structure.h"
#include "solve_short_characteristics_calculations.h"
#include "solve_short_characteristics_logic_function.h"

template<typename Type>
size_t ReadStartSettings(std::string name_file_settings, Type& class_file_vtk, std::string& name_file_vtk,
	std::string& name_file_sphere_direction, std::string& out_file_grid_vtk, std::string& main_dir/* ,
	std::string& name_file_in_faces, std::string& name_file_out_faces, std::string& name_file_count_out_faces, std::string& name_file_local_x0, std::string& name_file_x,
	std::string& name_file_s, std::string& name_file_id_neighbors, std::string& name_file_centers,
	std::string& name_file_dist_try, std::string& name_file_id_try, std::string& name_file_res, std::string& name_file_sizes*/) {

	std::ifstream ifile;
	ifile.open(name_file_settings);
	if (!ifile.is_open()) {
		std::cerr << " Error : file settings SolveShortCharacteristics is not open !\n";
		return 1;
	}

	std::string str; // переменная для перевода строки при чтении из файла

	ifile >> class_file_vtk;
	getline(ifile, str);
	getline(ifile, name_file_vtk);
	getline(ifile, name_file_sphere_direction);
	getline(ifile, out_file_grid_vtk);
	
	getline(ifile, main_dir);
	//getline(ifile, name_file_graph);
	
	/*getline(ifile, name_file_in_faces);
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

	getline(ifile, name_file_sizes);*/

	ifile.close();
	return 0;
}


template<typename Str_Type, typename T>
size_t ReadSimpleFile(const Str_Type name_file, std::vector<T>& data_array) {
	// Файл должен содержать в первой строке число элементов. Далее последовательные данные
	
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


#include <chrono>
#include <cstdio>
#include<inttypes.h>
#include <memory>

#ifdef _MSC_VER
#define fwrite_unlocked_fwrite_nolock
#endif

using namespace std::chrono;


#endif