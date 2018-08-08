#include "SMatrix.hpp"
#include "Matrix.hpp"
#include "Determinant.hpp"

using namespace std;

int main()
{
	SMatrix<int, 4, 4> mtx = { 	3, -3, -5, 	8,
							   -3,  2,  4, -6,
							    2, -5, -7,  5,
							   -4,  3,  5, -6 };
	Determinant<int> det{mtx};

	cout << det.det() << endl;

	return 0;
}