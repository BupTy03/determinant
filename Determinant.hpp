#pragma once
#include <iostream>
#include <string>
#include <stdexcept>
#include "Matrix.hpp"
#include "SMatrix.hpp"

struct Determinant_error : std::runtime_error
{
	explicit Determinant_error(const char* q) : std::runtime_error(q){}
	explicit Determinant_error(std::string& s) : std::runtime_error(s){}
};

template<typename T>
class Determinant
{
private:
	Matrix<T> mtx;
	Index rnk{0};
	T result;
	bool solved{false};

	T decomp_of_line_method()
	{
		T sum{};
		for (Index i = 0; i < rnk; ++i)
		{
			sum += mtx[0][i] * addition(1, i + 1);
		}
		return sum;
	}
	inline T star_method() noexcept
	{
		return 	(mtx[0][0] * mtx[1][1] * mtx[2][2] + mtx[0][1] * mtx[1][2] * mtx[2][0] + mtx[0][2] * mtx[1][0] * mtx[2][1]) -
				(mtx[0][2] * mtx[1][1] * mtx[2][0] + mtx[0][1] * mtx[1][0] * mtx[2][2] + mtx[0][0] * mtx[1][2] * mtx[2][1]);
	}

public:
	template<Index R>
	explicit Determinant(const SMatrix<T, R, R>& _smtx) : rnk{R}
	{
		mtx.resize(rnk, rnk);
		for(Index i = 0; i < rnk; ++i)
		{
			for(Index j = 0; j < rnk; ++j)
			{
				this->mtx[i][j] = _smtx[i][j];
			}
		}
	}
	explicit Determinant(const Matrix<T>& _mtx)
	{
		if(_mtx.size_dim1() == 0 || _mtx.size_dim2() == 0)
			throw Determinant_error{"empty matrix"};

		if(_mtx.size_dim1() != _mtx.size_dim2())
			throw Determinant_error{"dimensions should be equal"};

		for(Index i = 0; i < this->rnk; ++i)
		{
			for(Index j = 0; j < this->rnk; ++j)
			{
				this->mtx[i][j] = _mtx[i][j];
			}
		}
	}
	explicit Determinant(Matrix<T>&& _mtx)
	{
		if(_mtx.size_dim1() != _mtx.size_dim2())
			throw Determinant_error{"dimensions should be equal"};

		this->rnk = _mtx.size_dim1();
		this->mtx = std::move(_mtx);
	}
	Determinant minor(Index k, Index m)
	{
		k = k - 1;
		m = m - 1;
		Matrix<T> __mtx{mtx.size_dim1() - 1, mtx.size_dim2() - 1};
		auto dim1 = __mtx.size_dim1();
		auto dim2 = __mtx.size_dim2();
		for(decltype(dim1) i = 0, _i = 0; i <= dim1 && _i < dim1; ++i, ++_i)
		{
			if(i == k)
			{
				--_i;
				continue;
			}
			for(decltype(dim2) j = 0, _j = 0; j <= dim2 && _j < dim2; ++j, ++_j)
			{
				if(j == m)
				{
					--_j;
					continue;
				}

				__mtx[_i][_j] = mtx[i][j];
			}	
		}

		return Determinant{std::move(__mtx)};
	}
	T det()
	{
		if(solved)
		{
			return result;
		}

		switch (this->rnk)
		{
		case 1:
			result = mtx[0][0];
			break;
		case 2:
			result = mtx[0][0] * mtx[1][1] - mtx[0][1] * mtx[1][0];
			break;
		case 3:
			result = star_method();
			break;
		default:
			result = decomp_of_line_method();
			break;
		}
		solved = true;
		return result;
	}

	inline T addition(Index i, Index j)
	{
		T mul(((i + j) % 2) ? -1 : 1);
		return (minor(i, j)).det() * mul;
	}

	inline Matrix<T>& matrix() noexcept { return mtx; }
	inline Matrix<T> const& matrix() const noexcept { return mtx; }
	inline Index rank() const noexcept { return rnk; }
};
