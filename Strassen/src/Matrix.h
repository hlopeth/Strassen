#pragma once
#include <cstdint>
#include <vector>
#include <string>

class Matrix
{
public:
	Matrix(uint32_t with, uint32_t height);
	double& operator()(uint32_t row, uint32_t coll);
	double& operator()(uint32_t i);
	double operator()(uint32_t row, uint32_t coll) const;
	double operator()(uint32_t i) const;
	Matrix operator+(const Matrix& other);
	Matrix operator-(const Matrix& other);
	Matrix operator*(const Matrix& other);
	Matrix& operator+=(const Matrix& other);
	Matrix& operator-=(const Matrix& other);
	uint32_t colls() const;
	uint32_t rows() const;
private:
	const uint32_t _colls;
	const uint32_t _rows;
	std::vector<double> _data;
};

std::string to_string(const Matrix& val);
Matrix randomMatrix(uint32_t with, uint32_t height, double min, double max);