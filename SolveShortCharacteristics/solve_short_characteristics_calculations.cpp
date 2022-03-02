#include "solve_short_characteristics_calculations.h"


size_t ReadDataArray(const size_t class_file_vtk, const std::string name_file_vtk, int& size_grid,
	vtkDataArray*& density, vtkDataArray*& absorp_coef, vtkDataArray*& rad_en_loose_rate, const bool is_print/*=false*/) {

	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	if (ReadFileVtk(class_file_vtk, name_file_vtk, unstructured_grid, density, absorp_coef, rad_en_loose_rate, is_print)) {
		std::cout << "Error reading the file vtk\n";
		return 1;
	}

	size_grid = unstructured_grid->GetNumberOfCells();

	return 0;
}
size_t ReadFileVtk(const size_t class_file_vtk, const std::string name_file_vtk, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid,
	vtkDataArray*& density, vtkDataArray*& absorp_coef, vtkDataArray*& rad_en_loose_rate, const bool is_print/*=false*/) {

	/*������ ��������� ����� � ������ � vtkUnstructuredGrid*/

	vtkSmartPointer<vtkGenericDataObjectReader> reader_vtk =
		vtkSmartPointer<vtkGenericDataObjectReader>::New();
	reader_vtk->ReadAllScalarsOn();
	reader_vtk->SetFileName(name_file_vtk.c_str());
	reader_vtk->Update();

	if (reader_vtk->IsFileUnstructuredGrid()) {
		unstructuredgrid = reader_vtk->GetUnstructuredGridOutput();
		unstructuredgrid->Modified();
	}
	else {
		std::cout << "Error read file\n";
		std::cout << "file_vtk is not UnstructuredGrid\n";
		return 1;
	}

	switch (class_file_vtk) {
	case 0:
		density = NULL;
		absorp_coef = NULL;
		rad_en_loose_rate = NULL;
	case 1:
		density = unstructuredgrid->GetCellData()->GetScalars("alpha");
		absorp_coef = unstructuredgrid->GetCellData()->GetScalars("alpha");
		rad_en_loose_rate = unstructuredgrid->GetCellData()->GetScalars("Q");
		break;
	case 2:
		density = unstructuredgrid->GetCellData()->GetScalars("density");
		absorp_coef = unstructuredgrid->GetCellData()->GetScalars("absorp_coef");
		rad_en_loose_rate = unstructuredgrid->GetCellData()->GetScalars("radEnLooseRate");
		break;
	}

	if (is_print) {
		std::cout << "Grid has " << unstructuredgrid->GetNumberOfPoints() << " points.\n";
		std::cout << "Grid has " << unstructuredgrid->GetNumberOfCells() << " cells.\n";
		if (class_file_vtk) {
			std::cout << "density_Size: " << density->GetSize() << '\n';
			std::cout << "absorp_coef_Size: " << absorp_coef->GetSize() << '\n';
			std::cout << "Q_Size: " << rad_en_loose_rate->GetSize() << '\n';
		}
	}

	reader_vtk->GetOutput()->GlobalReleaseDataFlagOn();
	return 0;
}

size_t ReadSphereDirectionDecartToSpherical(const std::string name_file_sphere_direction, vector<Vector3>& directions_all, vector<Type>& squares, Type& square_surface) {

	ifstream ifile;

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

int InitGlobalValue(Vector3& start_point_plane_coord, Matrix3& transform_matrix, Matrix3& inverse_transform_matrix,
	Matrix3& straight_face, Matrix3& inclined_face) {
	// 3 ���� ������������
		{
			straight_face << 1. / 6, 1. / 6, 1,
				2. / 3, 1. / 6, 1,
				1. / 6, 2. / 3, 1;

			straight_face_inverse = straight_face.inverse(); // � �������
		}

		// 3 ���� ������������ �� ��������� ���������
		{
			inclined_face <<
				0, sqrt(2. / 3), 1,
				sqrt(2) / 4, 1. / (2 * sqrt(6)), 1,
				-sqrt(2) / 4, 1. / (2 * sqrt(6)), 1;

			inclined_face_inverse = inclined_face.inverse();  // � �������
		}


		//������� �������� �� ������������ ��������� � ���������� ��������� ��������� 
		{ transform_matrix <<
			-1. / sqrt(2), 1. / sqrt(2), 0,
			-1. / sqrt(6), -1. / sqrt(6), sqrt(2. / 3),
			1. / sqrt(3), 1. / sqrt(3), 1. / sqrt(3);
		}

		//������� �������� �� ��������� ��������� �  ���������� ������������ ���������
		{
			inverse_transform_matrix <<
				-1. / sqrt(2), -1. / sqrt(6), 1. / sqrt(3),
				1. / sqrt(2), -1. / sqrt(6), 1. / sqrt(3),
				0, sqrt(2. / 3), 1. / sqrt(3);
		}

		// ������ ���������� ���������
		start_point_plane_coord << 0.5, 0.5, 0;
		return 0;
}

int FindNeighborsPairFace(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<int>& all_pairs_face) {

	int count_unique_face = 0;
	const int N = unstructured_grid->GetNumberOfCells();
	all_pairs_face.resize(N * 4);
	for (int i = 0; i < N * 4; i++)
		all_pairs_face[i] = -2;

	vtkSmartPointer<vtkIdList> idp = vtkSmartPointer< vtkIdList>::New();
	vtkSmartPointer<vtkIdList> idc = vtkSmartPointer< vtkIdList>::New();

	int id_a, id_b, id_c;
	for (vtkIdType num_cell = 0; num_cell < N; ++num_cell) {

		for (int num_face = 0; num_face < 4; ++num_face) {
			if (all_pairs_face[num_cell * 4 + num_face] != -2) continue;
			++count_unique_face;

			idp = unstructured_grid->GetCell(num_cell)->GetFace(num_face)->GetPointIds();
			id_a = idp->GetId(0);
			id_b = idp->GetId(1);
			id_c = idp->GetId(2);

			/*����� ���� �������� � ����������� �� ������!*/
			unstructured_grid->GetCellNeighbors(num_cell, idp, idc);

			if (idc->GetNumberOfIds() == 1) {
				int id_neighbor_cell = idc->GetId(0);
				int id_neighbor_face = GetNumberNeighborFace(id_a, id_b, id_c, unstructured_grid->GetCell(id_neighbor_cell));
				all_pairs_face[num_cell * 4 + num_face] = id_neighbor_cell * 4 + id_neighbor_face;
				all_pairs_face[id_neighbor_cell * 4 + id_neighbor_face] = num_cell * 4 + num_face;
			}
			else if (idc->GetNumberOfIds() == 0)
				all_pairs_face[num_cell * 4 + num_face] = -1; // ��������� ������
			else
				std::cout << "More than 1 neighbor????\n";
		}

	}

	return count_unique_face;
}
int GetNumberNeighborFace(const int a, const int b, const int c, vtkCell* neighbor_cell) {

	vtkIdList* idc;

	int x, y, z;
	for (size_t i = 0; i < 4; i++)
	{
		idc = neighbor_cell->GetFace(i)->GetPointIds();
		x = idc->GetId(0);
		y = idc->GetId(1);
		z = idc->GetId(2);

		if (a == x && b == y && c == z) return i;
		else if (a == x && b == z && c == y) return i;
		else if (a == y && b == x && c == z) return i;
		else if (a == y && b == z && c == x) return i;
		else if (a == z && b == x && c == y) return i;
		else if (a == z && b == y && c == x) return i;

	}
	return -2;
}


int InitNodesValue(const std::vector<int>& all_pairs_face, std::vector<cell>& nodes_value) {


	const int n = all_pairs_face.size()/4;
	nodes_value.resize(n);

	for (size_t i = 0; i < n; ++i)
	{
		//nodes_value[i].id = i;	
		for (int j = 0; j < 4; j++) {
			nodes_value[i].neighbours_id_face[j] = all_pairs_face[i * 4 + j];
			nodes_value[i].nodes_value[j] = Vector3(-666, -666, -666);
		}
	}
	return 0;
}


int ResetNodesValue(std::vector<cell>& nodes_value) {


	const int n = nodes_value.size();
	
	for (size_t i = 0; i < n; ++i)
	{
		//nodes_value[i].id = i;	
		for (int j = 0; j < 4; j++) {
			nodes_value[i].nodes_value[j] = Vector3(-666, -666, -666);
		}
	}
	return 0;
}

int ReadNormalFile(std::string& name_file_normals, std::vector<Normals>& normals) {
	std::ifstream ifile;

	ifile.open(name_file_normals);
	if (!ifile.is_open()) {
		std::cout << "Error read file normals\n";
		return 1;
	}

	int N;
	ifile >> N;
	normals.resize(N);

	Normals norm(4);
	for (size_t i = 0; i < N; i++)
	{
		for (int j = 0; j < 4; j++)
			ifile >> norm.n[j][0] >> norm.n[j][1] >> norm.n[j][2];
		normals[i] = norm;
	}

	ifile.close();
	return 0;
}

Type NormIllum(const std::vector<Type>& Illum, const std::vector<Type>& Illum2) {
	Type max = -1;
	Type buf;
	for (size_t i = 0; i < Illum.size(); i++)
	{
		buf = fabs(Illum[i] - Illum2[i]);
		if (buf > max)
			max = buf;
	}
	return max;
}
Type NormIllumOmp(const std::vector<Type>& Illum, const std::vector<Type>& Illum2) {
	Type max = -1;
	
	const int n = Illum.size();
#pragma omp parallel default(none) shared(n, max, Illum, Illum2)
	{
		Type buf;
#pragma omp for
		for (size_t i = 0; i < n; i++)
		{
			buf = fabs(Illum[i] - Illum2[i]);

			if (buf > max)
			{
#pragma omp critical 
				{
					if (buf > max)
						max = buf;
				}
			}
		}
	}
	return max;
}

int ReadGraph(const std::string name_file_graph, std::vector<IntId>& sorted_id_cell) {
	ifstream ifile;
	ifile.open(name_file_graph);
	if (!ifile.is_open()) {
		std::cout << "Error open file\n";
		std::cout << "file_graph is not opened for reading\n";
		return 1;
	}

	IntId i;
	int count = 0;
	while (ifile >> i) {
		sorted_id_cell[count++] = i;
	}
	ifile.close();
	return 0;
}

int ReadGraphBin(const std::string name_file_graph, std::vector<IntId>& sorted_id_cell) {

	std::unique_ptr<FILE, int(*)(FILE*)> file_graph(fopen(name_file_graph.c_str(), "rb"), fclose);
	if (!file_graph) { printf("file_graph is not opened for writing\n"); return 1; }

	const int n = sorted_id_cell.size();
	fread_unlocked(sorted_id_cell.data(), sizeof(IntId), n, file_graph.get());

	fclose(file_graph.get());
	return 0;
}

int SetVertexMatrix(const size_t number_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, Eigen::Matrix4d& vertex_tetra) {

	// 4 ������� ������������(�� �������� � ������� � ������ ������)

	vtkPoints* points = unstructured_grid->GetCell(number_cell)->GetPoints();
	for (size_t j = 0; j < 3; j++) {
		vertex_tetra(j, 2) = points->GetPoint(2)[j];
		vertex_tetra(j, 0) = points->GetPoint(0)[j];
		vertex_tetra(j, 1) = points->GetPoint(1)[j];
		vertex_tetra(j, 3) = points->GetPoint(3)[j];
	}

	for (size_t i = 0; i < 4; i++)
		vertex_tetra(3, i) = 1;

	return 0;
}

int FindInAndOutFaces(const Vector3& direction, const int number_cell, const std::vector<Normals>& normals, int* face_state) {
	//face_state  -0=> ��������� �����,  1=> ��������  face_state.size=4!!!  

	Vector3 normal;

	for (size_t i = 0; i < 4; ++i) {

		normal = normals[number_cell].n[i];

		if (normal.dot(direction) < -eps)
			face_state[i] = 1;
		else
			face_state[i] = 0;  // ���� ����� �����������, �������, ��� ��� �� �������� ������������

	}

	return 0;
}

size_t CenterOfTetra(const int number_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, Vector3& point_in_tetra) {

	auto MakeS{ [](Type* P0,Type* P1,Type* P2) {
		Type Sum = 0;
		Type a[3], b[3];
		for (size_t i = 0; i < 3; i++) {
			a[i] = P1[i] - P0[i];
			b[i] = P2[i] - P0[i];
		}

		Sum = pow(a[1] * b[2] - a[2] * b[1], 2) + pow(a[0] * b[2] - a[2] * b[0], 2) + pow(a[0] * b[1] - a[1] * b[0], 2);
		return 0.5 * sqrt(Sum);
} };

	Type P0[3], P1[3], P2[3], P3[3];
	vtkSmartPointer<vtkIdList> idp = unstructuredgrid->GetCell(number_cell)->GetPointIds();

	unstructuredgrid->GetPoint(idp->GetId(0), P0);
	unstructuredgrid->GetPoint(idp->GetId(1), P1);
	unstructuredgrid->GetPoint(idp->GetId(2), P2);
	unstructuredgrid->GetPoint(idp->GetId(3), P3);

	Type Squr[4] = { MakeS(P1,P2,P3),MakeS(P0,P2,P3), MakeS(P0,P1,P3),MakeS(P0,P1,P2) };


	Type Sum = Squr[0] + Squr[1] + Squr[2] + Squr[3];
	for (size_t i = 0; i < 3; i++) {
		point_in_tetra[i] = (Squr[0] * P0[i] + Squr[1] * P1[i] + Squr[2] * P2[i] + Squr[3] * P3[i]) / Sum;
	}
	return 0;
}
size_t FindAllCenterOfTetra(const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, std::vector<Vector3>& point_in_tetra) {

	const int n = unstructuredgrid->GetNumberOfCells();
	point_in_tetra.resize(n);

	for (size_t i = 0; i < n; ++i)
		CenterOfTetra(i, unstructuredgrid, point_in_tetra[i]);

	return 0;
}

int FromGlobalToLocalTetra(const Eigen::Matrix4d& vertex_tetra, const Eigen::Vector3d& global_coord, Eigen::Vector3d& local_coord)
{
	// vertex_tetra -> [X;Y;Z;1]
	// �������� ���� ����� ������������ transpose ��-�� ������������� ������� ��������

	Eigen::Matrix4d vertex_tetra_inverse = vertex_tetra.inverse();
	// ����� � �����������������
	for (int i = 0; i < 3; ++i)
	{
		local_coord[i] = 0;
		for (int j = 0; j < 3; ++j)
			local_coord[i] += vertex_tetra_inverse(i, j) * global_coord[j];
		local_coord[i] += vertex_tetra_inverse(i, 3);
	}

	//local_coord = vertex_tetra * global_coord;
	return 0;
}
int FromLocalToGlobalTetra(const Eigen::Matrix4d& vertex_tetra, const Eigen::Vector3d& local_coord, Eigen::Vector3d& global_coord) {
	// vertex_tetra -> [X,Y,Z,1]
	// �������� � ������ �.�. �������������� ��������, � 4� ������ ��������� � �������� �� ������

	Type eta4 = 1 - local_coord[0] - local_coord[1] - local_coord[2];
	for (int i = 0; i < 3; ++i)
	{
		global_coord[i] = 0;
		for (int j = 0; j < 3; ++j)
			global_coord[i] += vertex_tetra(i, j) * local_coord[j];
		global_coord[i] += vertex_tetra(i, 3) * eta4;
	}

	//global_coord = vertex_tetra.inverse() * local_coord;

	return 0;
}

int FromTetraToPlane(const Eigen::Matrix3d& transform_matrix, const Eigen::Vector3d& start_point, const Eigen::Vector3d& tetra_coord, Eigen::Vector3d& plane_coord) {
	plane_coord = transform_matrix * (tetra_coord - start_point);
	return 0;
}
int FromPlaneToTetra(const Eigen::Matrix3d& inverse_transform_matrix, const Eigen::Vector3d& start_point, const Eigen::Vector3d& plane_coord,
	Eigen::Vector3d& tetra_coord) {
	tetra_coord = inverse_transform_matrix * plane_coord + start_point;
	return 0;
}

size_t IntersectionWithPlane(vtkCell* face, const Vector3& start_point, const Vector3& direction, Vector3& result) {

	//������� ������������

	Type A[3];
	Type B[3];
	Type C[3];
	face->GetPoints()->GetPoint(0, A);
	face->GetPoints()->GetPoint(1, B);
	face->GetPoints()->GetPoint(2, C);


	Type a, b, c, d;  // ��������� ��������� ���������
	Type t;

	a = A[1] * (B[2] - C[2]) + B[1] * (C[2] - A[2]) + C[1] * (A[2] - B[2]);
	b = A[0] * (C[2] - B[2]) + B[0] * (A[2] - C[2]) + C[0] * (B[2] - A[2]);
	c = A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]);
	d = A[0] * (C[1] * B[2] - B[1] * C[2]) + B[0] * (A[1] * C[2] - C[1] * A[2]) + C[0] * (B[1] * A[2] - A[1] * B[2]);

	t = -(a * start_point[0] + b * start_point[1] + c * start_point[2] + d) / (a * direction[0] + b * direction[1] + c * direction[2]);

	for (size_t i = 0; i < 3; ++i)
		result[i] = (direction[i] * t + start_point[i]);  // ����� ����������� ����  (start->direction) � ����������!!! face

	return 0;
}


size_t SetBasis(const Type* start_point, Vector3& normal, Matrix3& basis) {
	/*�� ��������� ����� � ������� ������ ��������� ����� ��������� ��������� (vec1, vec2).
	  ������� ����. ������ ���� ������ �����������(������������ normal). ������ ������ �� ���������� ������������*/
	Vector3 vec_1;
	Vector3 vec_2;

	if (abs(normal[1]) < 1e-20) {
		vec_1[0] = 0;
		vec_1[1] = 1;
		vec_1[2] = 0;
	}
	else {
		vec_1[0] = 1;
		vec_1[2] = 0;
		vec_1[1] = -(normal[0] * vec_1[0] + normal[2] * vec_1[2]) / normal[1];  //��-�� ���������� ������������ (N, vec1)==0
	}

	// ���������� ���������� ������ ���������
	if (normal[1] < 0)
		for (int i = 0; i < 3; ++i)
			vec_1[i] *= -1;

	// ������� ��������� ���������. Eigen �������� � �� �����!!!
	Eigen::Vector3d c = normal.cross(vec_1);

	for (size_t i = 0; i < 3; ++i)
		vec_2[i] = -c(i);

	vec_1.normalize();
	vec_2.normalize();

	basis.row(0) = vec_1;
	basis.row(1) = vec_2;
	basis.row(2) = normal;

	return 0;
}
size_t Make2dPoint(const Type* start, const Matrix3& local_basis, const Type* point, Vector3& new_point) {



	for (size_t i = 0; i < 3; i++)
		new_point[i] = 0;

	//������� 3d ����� � 2d (� ��������� ������ {start, local_basis}) 
	for (size_t k = 0; k < 3; k++) {
		new_point[0] += (point[k] - start[k]) * local_basis(0, k);
		new_point[1] += (point[k] - start[k]) * local_basis(1, k);
	}
	return 0;
}
size_t NormalToFace(const int num_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, int number_face, Vector3& n) {

	vtkSmartPointer<vtkIdList> idp = unstructured_grid->GetCell(num_cell)->GetFace(number_face)->GetPointIds();

	Type P0[3], P1[3], P2[3];

	unstructured_grid->GetPoints()->GetPoint(idp->GetId(0), P0);
	unstructured_grid->GetPoints()->GetPoint(idp->GetId(1), P1);
	unstructured_grid->GetPoints()->GetPoint(idp->GetId(2), P2);

	Type a[3], b[3];
	for (size_t i = 0; i < 3; i++) {
		a[i] = P1[i] - P0[i];
		b[i] = P2[i] - P0[i];
	}
	n[0] = a[1] * b[2] - a[2] * b[1];
	n[1] = -a[0] * b[2] + a[2] * b[0];
	n[2] = a[0] * b[1] - a[1] * b[0];


	n.normalize();

	vtkSmartPointer<vtkIdList> idp2 = unstructured_grid->GetCell(num_cell)->GetPointIds();

	size_t id;
	for (size_t i = 0; i < 4; i++) {
		int count = 0;
		for (size_t j = 0; j < 3; j++)
			if (idp2->GetId(i) != idp->GetId(j))
				count++;
		if (count == 3) {
			id = i;
			break;
		}

	}

	Type sum = 0;
	Type P3[3];
	unstructured_grid->GetCell(num_cell)->GetPoints()->GetPoint(idp2->GetId(id), P3);


	sum = P1[0] * (P2[1] - P3[1]) * P0[2] + P0[0] * (P3[1] - P2[1]) * P1[2] +
		P0[0] * (P1[1] - P3[1]) * P2[2] + P2[2] * (P1[0] * P3[1] - P1[0] * P0[1]) +
		P3[0] * (P0[2] * (P1[1] - P2[1]) + P1[2] * (P2[1] - P0[1]) + P2[2] * (P0[1] - P1[1]))
		+ P3[2] * (P1[0] * (P0[1] - P2[1]) + P0[0] * (P2[1] - P1[1])) +
		P2[0] * (P0[2] * (P3[1] - P1[1]) + P1[2] * (P0[1] - P3[1]) + P3[2] * (P1[1] - P0[1]));

	if (sum < 0)
		for (size_t i = 0; i < 3; i++)
			n[i] *= -1;
	return 0;
}
bool InTriangle(const int num_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cell_face, int number_face, const Eigen::Vector3d& XX) {
	/*face --- �����������, X --- ����� ��� ��������*/

	// ������� ������������
	Type AA[3];
	Type BB[3];
	Type CC[3];
	cell_face->GetFace(number_face)->GetPoints()->GetPoint(0, AA);
	cell_face->GetFace(number_face)->GetPoints()->GetPoint(1, BB);
	cell_face->GetFace(number_face)->GetPoints()->GetPoint(2, CC);


	Vector3 A, B, C, X;
	{
		Type Xx[3] = { XX[0],XX[1],XX[2] };
		Vector3 n;
		Matrix3 basis;
		NormalToFace(num_cell, unstructuredgrid, number_face, n);
		SetBasis(AA, n, basis);
		Make2dPoint(AA, basis, AA, A);
		Make2dPoint(AA, basis, BB, B);
		Make2dPoint(AA, basis, CC, C);
		Make2dPoint(AA, basis, Xx, X);
	}

	// �������� �������
	Type r1 = (A[0] - X[0]) * (B[1] - A[1]) - (B[0] - A[0]) * (A[1] - X[1]);
	Type r2 = (B[0] - X[0]) * (C[1] - B[1]) - (C[0] - B[0]) * (B[1] - X[1]);
	Type r3 = (C[0] - X[0]) * (A[1] - C[1]) - (A[0] - C[0]) * (C[1] - X[1]);

	if (r1 < 0 && r2 < 0 && r3 < 0)
		return true;
	else if (r1 > 0 && r2 > 0 && r3 > 0)
		return true;
	else return false;
}

Type BoundaryFunction(const int id_cell, const Vector3& x, const Vector3& direction, const std::vector<Type>& illum_old,
	const vector<Type>& directions, const vector<Type>& squares) {
	//Type betta = 0.5;

	//Type I0 = betta * GetSPart(id_cell, direction, illum_old, directions, squares, square_surface);

	return 0;
}
Type GetS(const int num_cell, const Vector3& direction, const std::vector<Type>& illum_old,
	const vector<Vector3>& directions, const vector<Type>& squares) {
	//num_cell equals x
	auto Gamma{ [](const Vector3& direction, const Vector3& direction2) {
	return (3. * (1 + pow(direction.dot(direction2),2))) / 4;
	} };


	Vector3 cur_direction;
	Type S = 0;
	const int N_dir = directions.size();
	const int N_cell = illum_old.size() / N_dir;

	for (int num_direction = 0; num_direction < N_dir; num_direction++)
	{
		S += Gamma(directions[num_direction], direction) * illum_old[num_direction * N_cell + num_cell] * squares[num_direction];
	}
	return S / square_surface;     // ���� *4PI, �� ��-�� ���������� Gamma ��������� �� 4PI
}

int CalculateInt(const int num_cells, const int num_directions, const std::vector<Type>& illum,
	const vector<Vector3>& directions, const vector<Type>& squares, vector<Type>& int_scattering) {

	Vector3 direction;
	for (int num_direction = 0; num_direction < num_directions; ++num_direction) {
		direction = directions[num_direction];
		for (size_t cell = 0; cell < num_cells; cell++)
		{
			int_scattering[num_direction * num_cells + cell] = GetS(cell, direction, illum, directions, squares);
		}
	}

	return 0;
}

int CalculateIntOmp(const int num_cells, const int num_directions, const std::vector<Type>& illum,
	const vector<Vector3>& directions, const vector<Type>& squares, vector<Type>& int_scattering) {

#pragma omp parallel default(none) shared(num_cells, num_directions, illum, directions, squares, int_scattering, square_surface)
		{
			Vector3 direction;

#pragma omp for
			for (int num_direction = 0; num_direction < num_directions; ++num_direction) {
				direction = directions[num_direction];
				for (size_t cell = 0; cell < num_cells; cell++)
				{
					int_scattering[num_direction * num_cells + cell] = GetS(cell, direction, illum, directions, squares);
				}
			}
		}
	return 0;
}


Vector3 GetInterpolationCoef(const Eigen::Matrix3d& interpolation_nodes, const Eigen::Vector3d& function_value) {
	//interpolation_nodes --- ��������� ���� ������������ (������ (� ����������� ������������ ���������) {x_i, y_i, 1})
	
	return interpolation_nodes.partialPivLu().solve(function_value);
	
}
Vector3 GetInterpolationCoefInverse(const Eigen::Matrix3d& interpolation_nodes, const Eigen::Vector3d& function_value) {
	//interpolation_nodes --- ��������� ���� ������������ (������ (� ����������� ������������ ���������) {x_i, y_i, 1})

	return interpolation_nodes * function_value;
}

int Min(const int a, const int b) {
	if (a < b) return a;
	return b;
}

size_t MakeEnergy(const vector<Type>& Illum, const vector<Type>& squares, const Type scuare_surface, vector<Type>& energy) {

	const int n = energy.size();

	for (size_t i = 0; i < n; ++i) {
		energy[i] = IntegarteDirection(i, Illum, squares, scuare_surface);
	}

	return 0;
}
Type IntegarteDirection(const int num_cell, const vector<Type>& Illum, const vector<Type>& squares, const Type scuare_surface) {
	Type res = 0;
	int n = squares.size();  // number of directions
	int m = Illum.size() / n;  // number of cells

	std::vector<pair<Type, Type>> I(n);


	for (size_t i = 0; i < n; i++) {
		I[i] = make_pair(Illum[m * i + num_cell], squares[i]);
	}

	auto cmp{ [](const pair<Type,Type> left, const pair<Type,Type> right) {

		return left.first < right.first;
} };
	std::sort(I.begin(), I.end(), cmp);

	for (size_t i = 0; i < n; i++) {
		res += I[i].first * I[i].second;
	}

	/*for (size_t i = 0; i < n; i++) {
		res += Illum[m * i + num_cell] * squares[i];
		}*/

	return res / scuare_surface;
}


size_t WriteFileSolution(const std::string name_file_out, const std::vector<Type>& vector_illum, const std::vector<Type>& vector_energy,
	 vtkSmartPointer<vtkUnstructuredGrid>& u_grid) {

	int n = u_grid->GetNumberOfCells();

	vtkSmartPointer<vtkDoubleArray> IllumArray =
		vtkSmartPointer<vtkDoubleArray>::New();

	vtkSmartPointer<vtkDoubleArray> EnergyArray =
		vtkSmartPointer<vtkDoubleArray>::New();

	for (size_t i = 0; i < n; i++) {
		EnergyArray->InsertNextTuple1(vector_energy[i]);
		IllumArray->InsertNextTuple1(vector_illum[i]);  // �� ������� �����������
	}

	EnergyArray->SetName("energy");
	u_grid->GetCellData()->AddArray(EnergyArray);
	
	IllumArray->SetName("illum");
	u_grid->GetCellData()->AddArray(IllumArray);


	vtkSmartPointer<vtkGenericDataObjectWriter> writer =
		vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	writer->SetFileTypeToBinary();
	writer->SetFileName(name_file_out.c_str());
	writer->SetInputData(u_grid);
	writer->Write();
	return 0;
}

size_t WriteFileSolution(const std::string name_file_out, const std::vector<Type>& vector_illum, const std::vector<Type>& vector_energy,
	const std::string& file_vtk) {

	vtkSmartPointer<vtkUnstructuredGrid> ugrid = 
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	vtkDataArray* foo;
	ReadFileVtk(0, file_vtk, ugrid, foo, foo, foo, false);

	WriteFileSolution(name_file_out, vector_illum, vector_energy, ugrid);

	return 0;
}



int ReadSizes(const std::string& name_file_size, int& countX, int& countX0, int& countOutC,
				int& countOut, int& countIn, int& countS, int& countRes, int& countTry) {

	std::ifstream ifile(name_file_size+".txt");

	if (!ifile.is_open()) {
		printf("Error read files size\n");
		return 1;
	}

	ifile >> countTry;
	ifile >> countX;
	ifile >> countX0;
	ifile >> countOutC;
	ifile >> countOut;
	ifile >> countIn;
	ifile >> countS;
	ifile >> countRes;

	ifile.close();
	return 0;
}


size_t ReadCompactFastGridData(const int count_dir, const int N,
	Str_Type& name_file_in_faces, Str_Type& name_file_out_faces, Str_Type& name_file_count_out_faces, Str_Type& name_file_local_x0, Str_Type& name_file_x,
	Str_Type& name_file_s, Str_Type& name_file_id_neighbors, Str_Type& name_file_centers,
	Str_Type& name_file_dist_try, Str_Type& name_file_id_try, Str_Type& name_file_res, Str_Type& name_file_sizes, Str_Type& name_file_graph,
	std::vector<cell>& grid, std::vector<int>& OutC,
	std::vector<int>& Out, std::vector<int>& In, std::vector<Type>& S, std::vector<Vector3>& X, std::vector<Vector2>& X0,
	std::vector<Type>& res_inner_bound, std::vector<int>& id_try_surface, vector<IntId>& sorted_id_cell) {


	std::unique_ptr<FILE, int(*)(FILE*)> file_graph(fopen((name_file_graph + ".bin").c_str(), "rb"), fclose);
	if (!file_graph) { printf("file_graph is not opened for writing\n"); return 1; }

	fread_unlocked(sorted_id_cell.data(), sizeof(IntId), sorted_id_cell.size(), file_graph.get());
	fclose(file_graph.get());


	grid.resize(N);
	int countX, countX0, countOutC, countOut, countIn, countS, countRes, countTry;
	ReadSizes(name_file_sizes, countX, countX0, countOutC, countOut, countIn, countS, countRes, countTry);

	//name_file_id_neighbors
	std::unique_ptr<FILE, int(*)(FILE*)> file_id_neighbors(fopen((name_file_id_neighbors + ".bin").c_str(), "rb"), fclose);
	if (!file_id_neighbors) PRINTF("Error file_id_neighbors\n")
		for (int i = 0; i < N; i++) {
			grid[i].neighbours_id_face.resize(4, -5);
			fread_unlocked(grid[i].neighbours_id_face.data(), sizeof(int), 4, file_id_neighbors.get());
		}
	fclose(file_id_neighbors.get());

	res_inner_bound.resize(countRes);
	id_try_surface.resize(countTry);
	OutC.resize(countOutC);
	Out.resize(countOut);
	In.resize(countIn);
	S.resize(countS);
	X.resize(countX);
	X0.resize(countX0);

	std::unique_ptr<FILE, int(*)(FILE*)> file_res(fopen((name_file_res + ".bin").c_str(), "rb"), fclose);
	if (!file_res) PRINTF("Error file_res\n")

		std::unique_ptr<FILE, int(*)(FILE*)> file_id(fopen((name_file_id_try + ".bin").c_str(), "rb"), fclose);
	if (!file_id) PRINTF("Error file_id\n")

		std::unique_ptr<FILE, int(*)(FILE*)> file_count_out(fopen((name_file_count_out_faces + ".bin").c_str(), "rb"), fclose);
	if (!file_count_out) { printf("file_graph is not opened for writing\n"); return 1; }

	std::unique_ptr<FILE, int(*)(FILE*)> file_out(fopen((name_file_out_faces + ".bin").c_str(), "rb"), fclose);
	if (!file_out) { printf("file_graph is not opened for writing\n"); return 1; }

	std::unique_ptr<FILE, int(*)(FILE*)> file_in(fopen((name_file_in_faces + ".bin").c_str(), "rb"), fclose);
	if (!file_in) { printf("file_graph is not opened for writing\n"); return 1; }

	std::unique_ptr<FILE, int(*)(FILE*)> file_x(fopen((name_file_x + ".bin").c_str(), "rb"), fclose);
	if (!file_x) { printf("file_graph is not opened for writing\n"); return 1; }

	std::unique_ptr<FILE, int(*)(FILE*)> file_x0(fopen((name_file_local_x0 + ".bin").c_str(), "rb"), fclose);
	if (!file_x0) { printf("file_graph is not opened for writing\n"); return 1; }

	std::unique_ptr<FILE, int(*)(FILE*)> file_s(fopen((name_file_s + ".bin").c_str(), "rb"), fclose);
	if (!file_s) { printf("file_graph is not opened for writing\n"); return 1; }


	fread_unlocked(res_inner_bound.data(), sizeof(Type), countRes, file_res.get());

	fread_unlocked(id_try_surface.data(), sizeof(int), countTry, file_id.get());

	fread_unlocked(S.data(), sizeof(Type), S.size(), file_s.get());

	fread_unlocked(X0.data(), sizeof(Vector2), X0.size(), file_x0.get());

	fread_unlocked(X.data(), sizeof(Vector3), X.size(), file_x.get());

	fread_unlocked(In.data(), sizeof(int), In.size(), file_in.get());

	fread_unlocked(Out.data(), sizeof(int), Out.size(), file_out.get());
	fread_unlocked(OutC.data(), sizeof(int), OutC.size(), file_count_out.get());


	fclose(file_res.get());
	fclose(file_id.get());
	fclose(file_count_out.get());
	fclose(file_out.get());
	fclose(file_in.get());
	fclose(file_x.get());
	fclose(file_x0.get());
	fclose(file_s.get());


	return 0;
}