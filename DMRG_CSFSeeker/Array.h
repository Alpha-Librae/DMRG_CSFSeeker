#ifndef ARRAY_CLASS_H
#define ARRAY_CLASS_H

#include <vector>

class Array
{
public:
	Array(int size = 0);
	Array(const Array& rhs);
	Array& operator=(const Array& rhs);
	~Array();
	int& operator[](int index);
	const int& operator[](int index) const;
	void initialize(); // 将 Data 全部设置为 0
	int size() const; // 返回数组的大小 Size
	friend class StateFunctions;
private:
	int Size;
	int* Data;
	void DFS(std::vector<Array>& seqs, int elec, int orbs, int LastSite = -1);
};

#endif // !ARRAY_CLASS_H