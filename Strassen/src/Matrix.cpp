#include "Matrix.h"
#include <cassert>
#include <sstream>
#include <iomanip>
#include <random>

Matrix::Matrix():
	_colls(0),
	_rows(0),
	_data(nullptr)
{
}

Matrix::Matrix(uint32_t with, uint32_t height):
	_colls(with),
	_rows(height)
{
	_data = new double[_colls * _rows];
}

Matrix::Matrix(const Matrix& m):
	_colls(m.colls()),
	_rows(m.rows())
{
	_data = new double[m.colls() * m.rows()];
	for (uint32_t i = 0; i < _colls * _rows; i++)
	{
		_data[i] = m(i);
	}
}

Matrix::Matrix(Matrix&& m)
{
	_colls = m._colls;
	_rows = m._rows;
	_data = m._data;

	m._colls = 0;
	m._rows = 0;
	m._data = nullptr;
}

double& Matrix::operator()(uint32_t row, uint32_t coll)
{
	assert(row >= 0 && row < _rows);
	assert(coll >= 0 && coll < _colls);
	return _data[row * _colls + coll];
}

double& Matrix::operator()(uint32_t i)
{
	assert(i >= 0 && i < _colls * _rows);
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
	assert(i >= 0 && i < _colls * _rows);
	return _data[i];
}

Matrix Matrix::operator+(const Matrix& other) const
{
	assert(_colls == other.colls() && _rows == other.rows());
	Matrix r(_colls, _rows);
	for (uint32_t i = 0; i < _rows * _colls; i++)
	{
		r(i) = _data[i] + other(i);
	}
	return r;
}

Matrix Matrix::operator-(const Matrix& other) const
{
	assert(_colls == other.colls() && _rows == other.rows());
	Matrix r(_colls, _rows);
	for (uint32_t i = 0; i < _rows * _colls; i++)
	{
		r(i) = _data[i] - other(i);
	}
	return r;
}

Matrix Matrix::operator*(const Matrix& b) const
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

Matrix& Matrix::operator=(Matrix&& other)
{
	this->_colls = other._colls;
	this->_rows = other._rows;
	this->_data = other._data;

	other._colls = 0;
	other._rows = 0;
	other._data = nullptr;
	return *this;
}

Matrix& Matrix::operator=(const Matrix& other)
{
	_colls = other._colls;
	_rows = other._rows;
	delete[] _data;
	_data = new double[_colls * _rows];
	for (uint32_t i = 0; i < _colls * _rows; i++)
	{
		_data[i] = other(i);
	}
	return *this;
}

bool Matrix::operator==(const Matrix& other)
{
	if (this->colls() != other.colls() || this->rows() != other.rows())
		return false;
	for (uint32_t i = 0; i < colls() * rows(); i++)
	{
		auto a = this->operator()(i);
		auto b = other(i);
		if (abs(this->operator()(i) - other(i)) > 0.0001)
			return false;
	}
	return true;
}

void Matrix::recize(uint32_t newColls, uint32_t newRows)
{
	double* newData = new double[newColls * newRows];
	for(uint32_t i = 0; i < _rows; i++)
		for (uint32_t j = 0; j < _colls; j++)
		{
			newData[i * newColls + j] = this->operator()(i, j);
		}

	_colls = newColls;
	_rows = newRows;
	delete[] _data;
	_data = newData;
	
}

uint32_t Matrix::colls() const
{
	return _colls;
}

uint32_t Matrix::rows() const
{
	return _rows;
}

void Matrix::split(const Matrix& src, Matrix& m11, Matrix& m12, Matrix& m21, Matrix& m22)
{
	uint32_t n = src._colls / 2;
	delete[] m11._data;
	m11._colls = n;
	m11._rows = n;
	m11._data = new double[n * n];

	m12._colls = n;
	m12._rows = n;
	delete[] m12._data;
	m12._data = new double[n * n];

	m21._colls = n;
	m21._rows = n;
	delete[] m21._data;
	m21._data = new double[n * n];

	m22._colls = n;
	m22._rows = n;
	delete[] m22._data;
	m22._data = new double[n * n];

	for(uint32_t i = 0; i < src.rows(); i++)
		for (uint32_t j = 0; j < src.colls(); j++)
		{
			if (i < n && j < n)
				m11(i, j) = src(i, j);
			if (i < n && j >= n)
				m12(i, j - n) = src(i, j);
			if (i >= n && j < n)
				m21(i - n, j) = src(i, j);
			if (i >= n && j >= n)
				m22(i - n, j - n) = src(i, j);
		}
}

Matrix Matrix::merge(const Matrix& m11, const Matrix& m12, const Matrix& m21, const Matrix& m22)
{
	uint32_t n = m11.colls();
	uint32_t m = m11.colls() * 2;
	Matrix dest;
	dest._colls = m;
	dest._rows = m;
	dest._data = new double[m * m];
	for (uint32_t i = 0; i < m; i++)
		for (uint32_t j = 0; j < m; j++)
		{
			if (i < n && j < n)
				dest(i, j) = m11(i, j);
			if (i < n && j >= n)
				dest(i, j) = m12(i, j - n);
			if (i >= n && j < n)
				dest(i, j) = m21(i - n, j);
			if (i >= n && j >= n)
				dest(i, j) = m22(i - n, j - n);
		}
	return dest;
}

Matrix::~Matrix()
{
	_rows = 0;
	_colls = 0;
	delete[] _data;
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

//
//void copy(const Matrix& src, Matrix& dest, uint32_t rowSrcOffs, uint32_t collSrcOffs, uint32_t rows, uint32_t colls)
//{
//	for(uint32_t i = 0; i < rows; i++)
//		for (uint32_t j = 0; j < colls; j++)
//		{
//			dest(i, j) = src(i + rowSrcOffs, j + collSrcOffs);
//		}
//}

