#include "Array.h"

Array::Array(int size) : Size(size)
{
	Data = new int[Size];
}

Array::Array(const Array & rhs)
{
	Size = rhs.Size;
	Data = new int[Size];
	for(int i = 0; i < Size; ++i)
	{
		Data[i] = rhs.Data[i];
	}
}

Array& Array::operator=(const Array & rhs)
{
	if(this != &rhs)
	{
		Size = rhs.Size;
		Data = new int[Size];
		for(int i = 0; i < Size; ++i)
		{
			Data[i] = rhs.Data[i];
		}
	}
	return *this;
}

Array::~Array()
{
	delete[] Data;
}

int& Array::operator[](int index)
{
	return Data[index];
}

const int& Array::operator[](int index) const
{
	return Data[index];
}

// 将 Data 全部设置为 0
void Array::initialize()
{
	for(int i = 0; i < Size; ++i)
	{
		Data[i] = 0;
	}
}

// 返回数组的大小 Size
int Array::size() const
{
	return Size;
}

// DFS-CAS(6,6)算法
void Array::DFS(std::vector<Array>& seqs, int elec, int orbs, int LastSite)
{
	if(elec == 0)
	{
		seqs.push_back(*this);
		return;
	}
	for(int i = LastSite + 1; i < orbs; ++i)
	{
		Data[i] = 1;
		DFS(seqs, elec - 1, orbs, i);
		Data[i] = 0;
	}
}