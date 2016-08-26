#ifndef RANDOM_NUMBER_CLASS_H
#define RANDOM_NUMBER_CLASS_H

#include<random>

template<typename Type>
struct Random
{
	Type operator()(Type left, Type right);
};

// template specialization: int
template<>
struct Random<int>
{
	// randomly return an int in the section [left, right]
	int operator()(int left, int right)
	{
		std::random_device rd;
		std::mt19937 eng(rd());
		std::uniform_int_distribution<int> num(left, right);
		return num(eng);
	}
};

// template specialization: double
template<>
struct Random<double>
{
	// randomly return a double in the section [left, right)
	double operator()(double left, double right)
	{
		std::random_device rd;
		std::mt19937 eng(rd());
		std::uniform_real_distribution<double> num(left, right);
		return num(eng);
	}
};

#endif // !RANDOM_NUMBER_CLASS_H