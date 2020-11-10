#pragma once
#include <cstdint>
#include <vector>
#include <string>

class Matrix
{
public:
	Matrix();
	Matrix(uint32_t with, uint32_t height);
	Matrix(const Matrix& m);
	Matrix(Matrix&& m);
	double& operator()(uint32_t row, uint32_t coll);
	double& operator()(uint32_t i);
	double operator()(uint32_t row, uint32_t coll) const;
	double operator()(uint32_t i) const;
	Matrix operator+(const Matrix& other) const;
	Matrix operator-(const Matrix& other) const;
	Matrix operator*(const Matrix& other) const;
	Matrix& operator+=(const Matrix& other);
	Matrix& operator-=(const Matrix& other);
	Matrix& operator=(Matrix&& other);
	Matrix& operator=(const Matrix& other);
	bool operator==(const Matrix& other);
	void recize(uint32_t newColls, uint32_t newRows);
	uint32_t colls() const;
	uint32_t rows() const;
	static void split(const Matrix& src, Matrix& m11, Matrix& m12, Matrix& m21, Matrix& m22);
	static Matrix merge(const Matrix& m11, const Matrix& m12, const Matrix& m21, const Matrix& m22);
	~Matrix();
private:
	uint32_t _colls;
	uint32_t _rows;
	double* _data;
};

std::string to_string(const Matrix& val);
Matrix randomMatrix(uint32_t with, uint32_t height, double min, double max);
//void copy(const Matrix& src, Matrix& dest, uint32_t rowSrcOffs, uint32_t collSrcOffs, uint32_t rows, uint32_t colls);