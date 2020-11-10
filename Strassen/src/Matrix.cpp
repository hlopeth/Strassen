#include "Matrix.h"
#include <cassert>
#include <sstream>
#include <iomanip>
#include <random>

Matrix::Matrix():
	_colls(0),
	_rows(0),
	_data(0)
{
}

Matrix::Matrix(uint32_t with, uint32_t height):
	_colls(with),
	_rows(height),
	_data(with * height)
{
}

Matrix::Matrix(const Matrix& m):
	_colls(m.colls()),
	_rows(m.rows()),
	_data(m.colls() * m.rows())
{
	for (uint32_t i = 0; i < _colls * _rows; i++)
	{
		_data[i] = m(i);
	}
}

double& Matrix::operator()(uint32_t row, uint32_t coll)
{
	assert(row >= 0 && row < _rows);
	assert(coll >= 0 && coll < _colls);
	return _data[row * _colls + coll];
}

double& Matrix::operator()(uint32_t i)
{
	assert(i >= 0 && i < _data.size());
	return _data[i];
}

double Matrix::operator()(uint32_t row, uint32_t coll) const
{
	assert(row >= 0 && row < _rows);
	assert(coll >= 0 && coll < _colls);
	return _data[row * _colls + coll];
}

double Matrix::operator()(uint32_t i) const
{
	assert(i >= 0 && i < _data.size());
	return _data[i];
}

Matrix Matrix::operator+(const Matrix& other)
{
	assert(_colls == other.colls() && _rows == other.rows());
	Matrix r(_colls, _rows);
	for (uint32_t i = 0; i < _rows * _colls; i++)
	{
		r(i) = _data[i] + other(i);
	}
	return r;
}

Matrix Matrix::operator-(const Matrix& other)
{
	assert(_colls == other.colls() && _rows == other.rows());
	Matrix r(_colls, _rows);
	for (uint32_t i = 0; i < _rows * _colls; i++)
	{
		r(i) = _data[i] - other(i);
	}
	return r;
}

Matrix Matrix::operator*(const Matrix& b)
{
	assert(_colls == b.rows());
	Matrix res(b.colls(), _rows);
	for (uint32_t r = 0; r < _rows; r++)
		for (uint32_t c = 0; c < b.colls(); c++)
		{
			double sum = 0;
			for (uint32_t i = 0; i < b.rows(); i++)
			{
				double a1 = b(i, c);
				double a2 = this->operator()(r, i);
				sum += this->operator()(r, i) * b(i, c);
			}
			res(r, c) = sum;
		}
	return res;
}

Matrix& Matrix::operator+=(const Matrix& other)
{
	assert(_colls == other.colls() && _rows == other.rows());
	for (uint32_t i = 0; i < _rows * _colls; i++)
	{
		_data[i] += other(i);
	}
	return *this;
}

Matrix& Matrix::operator-=(const Matrix& other)
{
	assert(_colls == other.colls() && _rows == other.rows());
	for (uint32_t i = 0; i < _rows * _colls; i++)
	{
		_data[i] -= other(i);
	}
	return *this;
}

uint32_t Matrix::colls() const
{
	return _colls;
}

uint32_t Matrix::rows() const
{
	return _rows;
}

std::string to_string(const Matrix& val)
{
	std::ostringstream os;
	os << "Matrix " << val.rows() << "x" << val.colls() << std::endl;
	for (uint32_t i = 0; i < val.rows(); i++)
	{
		for (uint32_t j = 0; j < val.colls(); j++ )
		{
			os << std::setw(8) << std::setprecision(2) << std::fixed << val(i, j);
		}
		os << std::endl;
	}
	return os.str();
}

inline double randomDouble() {
	static std::uniform_real_distribution<double> distribution(0.0, 1.0);
	static std::mt19937 generator;
	return distribution(generator);
}

inline double randomDouble(double min, double max) {
	// Returns a random real in [min,max).
	return min + (max - min) * randomDouble();
}

Matrix randomMatrix(uint32_t with, uint32_t height, double min, double max)
{
	Matrix m(with, height);
	for (uint32_t i = 0; i < with * height; i++)
	{
		m(i) = randomDouble(min, max);
	}
	return m;
}
