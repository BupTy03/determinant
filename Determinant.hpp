#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include "Matrix.hpp"
#include "SMatrix.hpp"

struct Determinant_error : std::runtime_error
{
	explicit Determinant_error(const char* q) : std::runtime_error(q){}
	explicit Determinant_error(const std::string& s) : std::runtime_error(s){}
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
		for (decltype(rnk) i = 0; i < rnk; ++i)
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
			throw Determinant_error{"wrong matrix"};

		if(_mtx.size_dim1() != _mtx.size_dim2())
			throw Determinant_error{"dimensions should be equal"};

		rnk = _mtx.size_dim1();
		mtx.resize(rnk, rnk);
		for(Index i = 0; i < rnk; ++i)
		{
			for(Index j = 0; j < rnk; ++j)
			{
				mtx[i][j] = _mtx[i][j];
			}
		}
	}
	explicit Determinant(Matrix<T>&& _mtx)
	{
		if(_mtx.size_dim1() != _mtx.size_dim2())
			throw Determinant_error{"dimensions should be equal"};

		rnk = _mtx.size_dim1();
		mtx = std::move(_mtx);
	}
	Determinant minor(Index k, Index m)
	{
		k = k - 1;
		m = m - 1;
		auto minor_rnk = rnk - 1;
		Matrix<T> __mtx(minor_rnk, minor_rnk);
		for(decltype(minor_rnk) i = 0, _i = 0; i <= minor_rnk && _i < minor_rnk; ++i, ++_i)
		{
			if(i == k)
			{
				--_i;
				continue;
			}
			for(decltype(minor_rnk) j = 0, _j = 0; j <= minor_rnk && _j < minor_rnk; ++j, ++_j)
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
	void transpose()
	{
		for(decltype(rnk) i = 1; i < rnk; ++i)
		{
			for(decltype(rnk) j = 1; j <= i; ++j)
			{
				std::swap(mtx[i][i - j], mtx[i - j][i]);
			}
		}
	}
	friend Determinant transposed(const Determinant& det)
	{
		Determinant __det = det;
		__det.transpose();
		return __det;
	}
	inline Matrix<T>& matrix() noexcept { return mtx; }
	inline Matrix<T> const& matrix() const noexcept { return mtx; }
	inline Index rank() const noexcept { return rnk; }
};
