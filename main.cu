#include <cuda_runtime.h>
#include "device_launch_parameters.h"

#include <stdio.h>
#include <iostream>
#include <fstream>

#include <string>
#include <vector>

#include <random>
#include <iomanip>

#define BS 32
typedef double Type;
using namespace std;

struct Vector3 {
    Type data[3];
    
    __host__ __device__ Type operator [](const int i)  const {
        if (i >= 3) {
            printf("Error Vector3\n");
            return 0;
        }

        return data[i];
    }
    
    Type& operator [](const int i) {
        if (i >= 3) {
            printf("Error Vector3\n");
            return *data;
        }
        return *(data + i);
    }
};

int CheckError(cudaError_t cudaStatus, const char* my_text = ""){
    if (cudaStatus != cudaSuccess) {
        printf("Error code: %d \n%s \n%s\n", cudaStatus, cudaGetErrorString(cudaStatus), my_text);
        return 1;
    }
    return 0;
}

Type NormError(const std::vector<Type>& A, const std::vector<Type>& B) {

    const int N = A.size();
    Type max = -1;
    Type buf;

    for (size_t i = 0; i < N; i++)
    {
        buf = fabs(A[i] - B[i]);
        if (buf > max)
            max = buf;
    }
    return max;
}

size_t ReadSphereDirectionDecartToSpherical(const std::string name_file_sphere_direction, vector<Vector3>& directions_all, vector<Type>& squares, Type& square_surface) {

    std::ifstream ifile;

    ifile.open(name_file_sphere_direction);
    if (!ifile.is_open()) {
        std::cout << "Error read file sphere direction\n";
        return 1;
    }
    int N = 0;
    ifile >> N;
    directions_all.resize(N);
    squares.resize(N);

    for (int i = 0; i < N; i++) {
        ifile >> squares[i];
        ifile >> directions_all[i][0] >> directions_all[i][1] >> directions_all[i][2];
    }
    ifile >> square_surface;
    ifile.close();

    return 0;
}

void GetS(const int N, const int M, const std::vector<Type>& illum_old, vector<Type>& Integ,
    const vector<Vector3>& directions, const vector<Type>& squares, const Type square_surface) {
    //num_cell equals x

    auto Gamma{ [](const Vector3& direction, const Vector3& direction2) {

        Type sum = 0;
        for (size_t i = 0; i < 3; i++)
        {
            sum += direction[i] * direction2[i];
        }

    return (3. * (1 + sum * sum)) / 4;
    } };


    Integ.resize(N * M, 0);

    for (size_t k = 0; k < M; k++)
        for (size_t i = 0; i < N; i++)
        {

            for (int num_direction = 0; num_direction < M; num_direction++)
            {
                Integ[i * M + k] += Gamma(directions[num_direction], directions[k]) *
                    illum_old[i * M + num_direction] * squares[num_direction];
            }

            Integ[i * M + k] / square_surface;

        }
}

__device__ Type Gamma(const Vector3& direction, const Vector3& direction2) {

    Type sum = 0;
    for (size_t i = 0; i < 3; i++)
    {
        sum += direction[i] * direction2[i];
    }

    return (3. * (1 + sum * sum)) / 4;
}

__global__ void d_GetS(const int N, const int M, Type* illum_old, Type* Integ,
    const Type* d_directions, Type* squares, const Type square_surface) {

    const Vector3* directions = (Vector3*)d_directions;

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= N) return;
    if (k >= M) return;


    for (int num_direction = 0; num_direction < M; num_direction++)
    {
        Integ[i * M + k] += Gamma(directions[num_direction], directions[k]) *
            illum_old[i * M + num_direction] * squares[num_direction];
    }

    Integ[i * M + k] / square_surface;
}

int main()
{       
    vector<Vector3> directions;
    vector<Type> squares;
    Type square_surface;
    
    ReadSphereDirectionDecartToSpherical("D:\\Desktop\\FilesCourse\\Grids\\Precision\\surface_126_dir.txt", directions, squares, square_surface);
    
    int N = 1e5;
    int M = directions.size();

    std::vector<Type> h_Illum(N * M, 0);

    for (size_t i = 0; i < N * M; i++){
        h_Illum[i] = rand();
    }

    std::vector<Type> Integ(N * M, 0);

    GetS(N, M, h_Illum, Integ, directions, squares, square_surface);

    Type* d_Integ;
    Type* d_directions;
    Type* d_squares;
    Type* d_illum;

    cudaEvent_t start, finish;  //для засечки времени

     // выделение памяти и перессылка на device исходного вектора Y
    {
        //проверка наличия карты 
        if (CheckError(cudaSetDevice(0), "cudaSetDevice failed!Do you have a CUDA - capable GPU installed ?")) return 1;


        //======================================buf_device_host==============================================================//

        if (CheckError(cudaMalloc(&d_Integ, N * M * sizeof(Type)), "cudaMalloc failed!")) return 1;
        if (CheckError(cudaMalloc(&d_illum, N * M * sizeof(Type)), "cudaMalloc failed!")) return 1;
        if (CheckError(cudaMalloc(&d_directions, M * 3 * sizeof(Type)), "cudaMalloc failed!")) return 1;
        if (CheckError(cudaMalloc(&d_squares, M * sizeof(Type)), "cudaMalloc failed!")) return 1;

        //======================================CreateEvent==============================================================//

        cudaEventCreate(&start);
        cudaEventCreate(&finish);
    }

    // перессылка данных на device
    {
        if (CheckError(cudaMemcpy(directions.data(), d_directions, M * 3 * sizeof(Type), cudaMemcpyHostToDevice), "cudaMemcpy failed!")) return 1;
        if (CheckError(cudaMemcpy(h_Illum.data(), d_illum, M * N * sizeof(Type), cudaMemcpyHostToDevice), "cudaMemcpy failed!")) return 1;
        if (CheckError(cudaMemcpy(squares.data(), d_squares, M * sizeof(Type), cudaMemcpyDeviceToHost), "cudaMemcpy failed!")) return 1;
    }

    cudaEventRecord(start);
    cudaEventSynchronize(start);

   
    dim3 threads(BS, BS);
    dim3 blocks((N + BS - 1) / BS, (M + BS - 1) / BS);

    d_GetS << <blocks, threads >> > (N, M, d_illum, d_Integ, d_directions, d_squares, square_surface);
     
    
    cudaEventRecord(finish);
    cudaEventSynchronize(finish);

    float time;
    cudaEventElapsedTime(&time, start, finish);


    // Check for any errors launching the kernel
    if (CheckError(cudaGetLastError(), "calculate_n_body failed!")) return 1;

    // cudaDeviceSynchronize waits for the kernel to finish, and returns any errors encountered during the launch.
    if (CheckError(cudaDeviceSynchronize(), "cudaDeviceSynchronize returned error code %d after launching global functions")) return 1;


    std::vector<Type> device_Integ(N * M, 0);
    //копирование на host 
    {
        if (CheckError(cudaMemcpy(device_Integ.data(), d_Integ, N * M * sizeof(Type), cudaMemcpyDeviceToHost), "cudaMemcpy failed!")) return 1;
    }

    // прекращение работы с device
    {
        cudaEventDestroy(start);
        cudaEventDestroy(finish);

        cudaFree(d_Integ);
        cudaFree(d_squares);
        cudaFree(d_directions);
        cudaFree(d_illum);

        // cudaDeviceReset must be called before exiting in order for profiling and tracing tools such as Nsight and Visual Profiler to show complete traces.
        if (CheckError(cudaDeviceReset(), "cudaDeviceReset failed!")) return 1;
    }

    std::cout << "error host and device:= " << NormError(Integ, device_Integ) << '\n';

    printf("GPU time = %f\n", time / 1000.0f);
    return 0;
}