#include <vector>
#include <chrono>
#include <cstdio>
#include <inttypes.h>
#include <memory>
#include <stdio.h>
#include <string>
#include <iostream>
typedef std::string Str_Type;
typedef double Type;

#ifdef _MSC_VER
#define fwrite_unlocked _fwrite_nolock
#define fread_unlocked  _fread_nolock
#endif

#define PRINTF(a) {printf(a); return 1;}

using namespace std;
using namespace std::chrono;
#include <iostream>

int main(int argc, char* argv[])
{	
	std::string name_file_res = "";
	int count_dir = 126;

	if (argc <= 1) //"D:\\Desktop\\FilesCourse\\IllumGrid\\testConfig\\id_defining_faces";	
		name_file_res = "D:\\Desktop\\FilesCourse\\IllumGrid\\testConfig\\ResBound";	
	else
		name_file_res = argv[1];
	if (argc > 2)
		count_dir = std::stoi(argv[2]);

	
	
	std::vector<Type> res;
	std::vector<Type> full_res;


	for (size_t i = 0; i < count_dir; i++)
	{

		int n;
		std::unique_ptr<FILE, int(*)(FILE*)> file_res(fopen((name_file_res + std::to_string(i) + ".bin").c_str(), "rb"), fclose);
		if (!file_res) PRINTF("Error file_res\n")

		fread_unlocked(&n, sizeof(int), 1, file_res.get());
		res.resize(n);
		fread_unlocked(res.data(), sizeof(Type), n, file_res.get());

		full_res.insert(full_res.end(), res.begin(), res.end());
		
		fclose(file_res.get());
	}

	std::unique_ptr<FILE, int(*)(FILE*)> file_res_full(fopen((name_file_res + ".bin").c_str(), "wb"), fclose);
	if (!file_res_full) PRINTF("Error file_res_full\n")

	int n = full_res.size();
	fwrite_unlocked(&n, sizeof(int), 1, file_res_full.get());
	fwrite_unlocked(full_res.data(), sizeof(Type), n, file_res_full.get());
	fclose(file_res_full.get());

	return 0;
}
