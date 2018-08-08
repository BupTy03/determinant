#include "SMatrix.hpp"
#include "Matrix.hpp"
#include "Determinant.hpp"

using namespace std;

int main()
{
	SMatrix<int, 5, 5> mtx;

	int count = 0;
	for(auto& it : mtx)
		it = ++count;

	Determinant<int> det{mtx};

	cout << "Matrix: " << det.matrix() << endl;
	cout << "Transposed: " << (transposed(det)).matrix() << endl;
	cout << "Determinant: " << det.det() << endl;
	cout << "Minor(2, 2): " << (det.minor(2, 2)).matrix() << endl;
	cout << "Addition(2, 2): " << det.addition(2, 2) << endl;

	return 0;
}
