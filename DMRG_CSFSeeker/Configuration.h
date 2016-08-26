#ifndef CONFIGURATION_CLASS_H
#define CONFIGURATION_CLASS_H

#include <iostream>
#include <string>
#include "RandomNumber.h"
#include "Array.h"

// 策略:
// 组态由 2Norbs 个自旋轨道以 ababab... 这样的序列记录。
// 序列的每一位用 0 表示无电子，用 1 表示有电子
class Configuration
{
public:
	// constructors
	// 构造函数：接受的参数为活性空间的轨道数、电子数以及组态的自旋多重度
	// 自旋多重度 S=2s+1 的默认值为 1，即单重态。
	// 因为增加了根据给定的序列生成组态的功能，故初始设定的电子数也不再重要，从构造函数中可以省去。
	Configuration(int orbs, int elec = 0, int spin = 1);
	Configuration(const Configuration& rhs); // 复制构造函数
	// 用于中心扩展方法的构造函数
	Configuration(const Array& CoreAlpha, const Array& CoreBeta,
		int orbs, int elec, int spin,
		int CASCoreOrbs, int CASCoreElecAlpha, int CASCoreElecBeta);
	// 测试功能
	Configuration(const Array& CoreAlpha, const Array& CoreBeta,
		int orbs, int elec, int spin,
		int CASCoreOrbs, int CASCoreElecAlpha, int CASCoreElecBeta,
		int LargerSpaceOrbs, int NLargerSpaceElecAlpha, int NLargerSpaceElecBeta);
	Configuration(const std::string& rhs); // 由 444111 这种序列构成一个组态
	// the deconstructor
	~Configuration();
	// member functions
	int& alpha(int index); // 返回索引为 index 的 alpha 自旋轨道的状态
	const int& alpha(int index) const; // 右值返回索引为 index 的 alpha 自旋轨道的状态
	int& beta(int index); // 返回索引为 index 的 beta 自旋轨道的状态
	const int& beta(int index) const; // 右值返回索引为 index 的 beta 自旋轨道的状态
	int& operator[](int index); // 返回索引为 index 的自旋轨道的状态
	const int& operator[](int index) const; // 右值返回索引为 index 的自旋轨道的状态
	double coeff() const; // 返回组态系数
	void initialize(); // 将 Sequence 的每一位都置 0
	void SetHartreeFockDet(); // 将组态设定为 Hartree-Fock 行列式表示的最低占据形式
	void SetBySeed(const std::string& rhs); // 由 444111 这种序列构成一个组态
	const Configuration& randomize(); // 随机生成电子排布序列，按照轨道数、电子数和自旋多重度，遵循 Paili 不相容原理
	Configuration& mutate(double rate, double StepSize = 0.1); // 按照给定的概率发生变异
	Configuration& MutateStrategy2(); // 变异策略2
	Configuration& operator=(const Configuration& rhs); // 重载赋值操作符
	// 重载输出操作符
	friend std::ostream& operator<<(std::ostream& os, const Configuration& rhs);
	// 友元类：态函数
	friend class StateFunctions;
	// 用于按照组态系数绝对值从大到小进行排序的函数对象
	friend struct SortByCoefficient;
	// 用于按照序列进行排序的函数对象
	friend struct SortBySequence;
	// 用于执行 unique 操作，定义相等操作符
	friend bool operator==(const Configuration& lhs, const Configuration& rhs);
private:
	int SpinOrbs; // 自旋轨道的个数
	int Nelec; // 总的电子数
	int Spin; // 自旋多重度，S=2s+1
	double Coefficient; // 组态系数，作为其适应度的表征
	// ababab... 这样的序列，每一位用 0 表示无电子，用 1 表示有电子
	int* Sequence;
	// 私有辅助函数
	// 根据已知的序列 Sequence 计算总的电子数 Nelec 和自旋多重度 Spin
	void ComputeNelecSpin();
	// 从指定空间轨道位置 site 开始，按照 status（1,2,3,4）连续填充 Sequence 中的两个位置
	void SetOccupationStatus(int site, int status);
};

// 用于将组态按照系数绝对值从大到小进行排序的函数对象
struct SortByCoefficient
{
	bool operator()(const Configuration& lhs, const Configuration& rhs);
};

// 用于将组态按照序列进行唯一排序的函数对象。顺序并不重要，只需要保证相同序列排序后连续。
struct SortBySequence
{
	bool operator()(const Configuration& lhs, const Configuration& rhs);
};

#endif // !CONFIGURATION_CLASS_H