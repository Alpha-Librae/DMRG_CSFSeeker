#include "Configuration.h"
#include <cmath>
#include <stdexcept>

// constructors
// 构造函数：接受的参数为活性空间的轨道数、电子数以及组态的自旋多重度(默认为 1，即单重态)
Configuration::Configuration(int orbs, int elec, int spin)
{
	SpinOrbs = 2 * orbs; // orbs 为空间轨道的数目
	Nelec = elec;
	Spin = spin;
	Sequence = new int[SpinOrbs];
}

// the copy constructor
Configuration::Configuration(const Configuration& rhs)
{
	SpinOrbs = rhs.SpinOrbs;
	Nelec = rhs.Nelec;
	Spin = rhs.Spin;
	Coefficient = rhs.Coefficient;
	Sequence = new int[SpinOrbs];
	for(int ix = 0; ix != SpinOrbs; ++ix)
	{
		Sequence[ix] = rhs.Sequence[ix];
	}
}

Configuration::Configuration(const Array& CoreAlpha, const Array& CoreBeta,
	int orbs, int elec, int spin, int CASCoreOrbs, int CASCoreElecAlpha, int CASCoreElecBeta)
{
	SpinOrbs = 2 * orbs;
	Sequence = new int[SpinOrbs];
	Nelec = elec;
	Spin = spin;

	int OuterAlpha = (Nelec + Spin - 1) / 2 - CASCoreElecAlpha; // Core 空间外 alpha 电子的数目
	int OuterBeta = (Nelec - Spin + 1) / 2 - CASCoreElecBeta; // Core 空间外 beta 电子的数目

	initialize(); // 先将所有的位置都标记为没有电子
	// 先填充核心 CAS 空间左侧的全占据部分
	int idx = 0;
	while(idx < OuterAlpha)
	{
		alpha(idx) = 1;
		++idx;
	}
	// 然后按照 CoreAlpha 序列填充剩下的 alpha 电子
	for(int i = 0; i != CoreAlpha.size(); ++i)
	{
		alpha(idx) = CoreAlpha[i];
		++idx;
	}
	// 对于 beta 电子，执行类似的操作
	idx = 0;
	while(idx < OuterBeta)
	{
		beta(idx) = 1;
		++idx;
	}
	for(int i = 0; i != CoreBeta.size(); ++i)
	{
		beta(idx) = CoreBeta[i];
		++idx;
	}
}

// 测试功能
Configuration::Configuration(const Array& CoreAlpha, const Array& CoreBeta,
	int orbs, int elec, int spin,
	int CASCoreOrbs, int CASCoreElecAlpha, int CASCoreElecBeta,
	int LargerSpaceOrbs, int NLargerSpaceElecAlpha, int NLargerSpaceElecBeta)
{
#ifdef DEBUG
	if(LargerSpaceOrbs > orbs || NLargerSpaceElecAlpha + NLargerSpaceElecBeta > elec)
	{
		throw std::runtime_error("ERROR: invalid inner space configuration!");
	}
#endif // DEBUG
	SpinOrbs = 2 * orbs;
	Sequence = new int[SpinOrbs];
	Nelec = elec;
	Spin = spin;

	int NElecAlpha = (Nelec + Spin - 1) / 2; // 总的 alpha 电子数
	int NElecBeta = (Nelec - Spin + 1) / 2; // 总的 beta 电子数
	int NREAlpha = NLargerSpaceElecAlpha - CASCoreElecAlpha; // 在次小空间内进行随机填充的 alpha 电子数
	int NREBeta = NLargerSpaceElecBeta - CASCoreElecBeta; // 在次小空间内进行随机填充的 beta 电子数
	int NFEAlpha = NElecAlpha - NLargerSpaceElecAlpha; // 小空间左侧固定填充位置的 alpha 电子数
	int NFEBeta = NElecBeta - NLargerSpaceElecBeta; // 小空间左侧固定填充位置的 beta 电子数
	Random<int> rand; // 确定随机填充位置的随机数

	initialize(); // 先将所有的轨道都填 0

	// alpha 电子
	// 先填充全占据部分
	int idx = 0;
	while(idx < NFEAlpha)
	{
		alpha(idx) = 1;
		++idx;
	}
	// 随机填充的左半部分
	int NRFSites = LargerSpaceOrbs - CASCoreOrbs; // 随机填充的轨道的数目
	int site = 0; // 记录随机到的位置
	Array RFOrbitals(NRFSites);
	RFOrbitals.initialize();
	for(int i = 0; i < NREAlpha; ++i)
	{
		site = rand(0, NRFSites - 1);
		if(RFOrbitals[site] == 1)
		{
			--i;
		}
		else
		{
			RFOrbitals[site] = 1;
		}
	}
	int NRFOrbsLeft = (LargerSpaceOrbs + (LargerSpaceOrbs % 2)) / 2;
	int NCASOrbsLeft = (CASCoreOrbs + (CASCoreOrbs % 2)) / 2;
	int NRFOrbsO = NRFOrbsLeft - NCASOrbsLeft;
	site = 0;
	for(int i = 0; i < NRFOrbsO; ++i)
	{
		alpha(idx) = RFOrbitals[site];
		++site;
		++idx;
	}
	// 小空间 CAS 穷举部分
	for(int i = 0; i < CoreAlpha.size(); ++i)
	{
		alpha(idx) = CoreAlpha[i];
		++idx;
	}
	// 随机填充的右半部分
	int NRFOrbsRight = (LargerSpaceOrbs - (LargerSpaceOrbs % 2)) / 2;
	int NCASOrbsRight = (CASCoreOrbs - (CASCoreOrbs % 2)) / 2;
	int NRFOrbsU = NRFOrbsRight - NCASOrbsRight;
	for(int i = 0; i < NRFOrbsU; ++i)
	{
		alpha(idx) = RFOrbitals[site];
		++site;
		++idx;
	}

	// beta 电子
	// 先填充全占据部分
	idx = 0;
	while(idx < NFEBeta)
	{
		beta(idx) = 1;
		++idx;
	}
	// 随机填充的左半部分
	site = 0; // 记录随机到的位置
	RFOrbitals.initialize();
	for(int i = 0; i < NREBeta; ++i)
	{
		site = rand(0, NRFSites - 1);
		if(RFOrbitals[site] == 1)
		{
			--i;
		}
		else
		{
			RFOrbitals[site] = 1;
		}
	}
	site = 0;
	for(int i = 0; i < NRFOrbsO; ++i)
	{
		beta(idx) = RFOrbitals[site];
		++site;
		++idx;
	}
	// 小空间 CAS 穷举部分
	for(int i = 0; i < CoreBeta.size(); ++i)
	{
		beta(idx) = CoreBeta[i];
		++idx;
	}
	// 随机填充的右半部分
	for(int i = 0; i < NRFOrbsU; ++i)
	{
		beta(idx) = RFOrbitals[site];
		++site;
		++idx;
	}
}

// 由 444111 这种序列构成一个组态
Configuration::Configuration(const std::string & rhs)
{
	SpinOrbs = 2 * rhs.length();
	Sequence = new int[SpinOrbs];
	for(int ix = 0; ix != rhs.length(); ++ix)
	{
		SetOccupationStatus(ix, rhs[ix] - '0');
	}
	// 计算电子个数和自旋多重度
	ComputeNelecSpin();
}

// the deconstructor
Configuration::~Configuration()
{
	delete[] Sequence;
}

// 返回索引为 index 的 alpha 自旋轨道的状态
int& Configuration::alpha(int index)
{
#ifdef _INPUT_CHECK
	if(index < 0 || index > SpinOrbs / 2)
	{
		throw std::out_of_range("ERROR: the requested Alpha spin orbital does not exist!");
	}
#endif // _INPUT_CHECK
	return Sequence[2 * index];
}

// 右值返回索引为 index 的 alpha 自旋轨道的状态
const int& Configuration::alpha(int index) const
{
#ifdef _INPUT_CHECK
	if(index < 0 || index > SpinOrbs / 2)
	{
		throw std::out_of_range("ERROR: the requested Alpha spin orbital does not exist!");
	}
#endif // _INPUT_CHECK
	return Sequence[2 * index];
}

// 返回索引为 index 的 beta 自旋轨道的状态
int& Configuration::beta(int index)
{
#ifdef _INPUT_CHECK
	if(index < 0 || index > SpinOrbs / 2)
	{
		throw std::out_of_range("ERROR: the requested Alpha spin orbital does not exist!");
	}
#endif // _INPUT_CHECK
	return Sequence[2 * index + 1];
}

// 右值返回索引为 index 的 beta 自旋轨道的状态
const int& Configuration::beta(int index) const
{
#ifdef _INPUT_CHECK
	if(index < 0 || index > SpinOrbs / 2)
	{
		throw std::out_of_range("ERROR: the requested Alpha spin orbital does not exist!");
	}
#endif // _INPUT_CHECK
	return Sequence[2 * index + 1];
}

// 返回索引为 index 的自旋轨道的状态
int& Configuration::operator[](int index)
{
#ifdef _INPUT_CHECK
	if(index < 0 || index >= SpinOrbs)
	{
		throw std::out_of_range("ERROR: the requested Alpha spin orbital does not exist!");
	}
#endif // _INPUT_CHECK
	return Sequence[index];
}

// 右值返回索引为 index 的自旋轨道的状态
const int& Configuration::operator[](int index) const
{
#ifdef _INPUT_CHECK
	if(index < 0 || index >= SpinOrbs)
	{
		throw std::out_of_range("ERROR: the requested Alpha spin orbital does not exist!");
	}
#endif // _INPUT_CHECK
	return Sequence[index];
}

// member functions
// 返回组态系数
double Configuration::coeff() const
{
	return Coefficient;
}

// 将 Sequence 的每一位都置 0
void Configuration::initialize()
{
	for(int ix = 0; ix != SpinOrbs; ++ix)
	{
		Sequence[ix] = 0;
	}
}

// 将组态设定为 Hartree-Fock 行列式表示的最低占据形式
void Configuration::SetHartreeFockDet()
{
	int Alpha = (Nelec + Spin - 1) / 2; // alpha 电子的数目
	int Beta = (Nelec - Spin + 1) / 2; // beta 电子的数目
	initialize(); // 先将全部的位置标记为没有电子
	for(int ix = 0; ix != Alpha; ++ix) // 然后填充 alpha 电子
	{
		Sequence[2 * ix] = 1;
	}
	for(int ix = 0; ix != Beta; ++ix) // 最后填充 beta 电子
	{
		Sequence[2 * ix + 1] = 1;
	}
}

// 由 444111 这种字符串序列构成一个组态
void Configuration::SetBySeed(const std::string & rhs)
{
	int Norbs = SpinOrbs / 2;
	int Length;
	if(Norbs <= rhs.length()) // 如果字符序列太长，则忽略长出的部分
	{
		Length = Norbs;
	}
	else // 如果字符序列太短，则先给序列中没有描述的部分补0，表示没有电子
	{
		Length = rhs.length();
		for(int ix = 2 * Length; ix != SpinOrbs; ++ix)
		{
			Sequence[ix] = 0;
		}
	}
	for(int ix = 0; ix != Length; ++ix)
	{
		SetOccupationStatus(ix, rhs[ix] - '0');
	}
	// 计算电子个数和自旋多重度
	ComputeNelecSpin();
}

const Configuration& Configuration::randomize()
{
	Random<int> randSite; // 整数随机数生成器
	int Alpha = (Nelec + Spin - 1) / 2; // alpha 电子的数目
	int Beta = (Nelec - Spin + 1) / 2; // beta 电子的数目
	int Norbs = SpinOrbs / 2; // 空间轨道的数目
	int Site; // 记录随机到的位置
			  // 存储的序列是 ababab... 这种 alpha 与 beta 交替的序列
	initialize(); // 首先将序列初始化为全 0，表示所有的位置都没有电子
				  // 在 Norbs 个空间轨道中随机选择 Alpha 个，填充 alpha 电子
	for(int ix = 0; ix != Alpha; ++ix)
	{
		Site = 2 * randSite(0, Norbs - 1); // alpha 轨道在前，其索引为偶数
		if(Sequence[Site] == 1) // 如果当前位置已经是 1，说明此位置已经填过电子，那么当前循环无效
		{
			--ix; // 使当前循环无效，需要将循环计数器退一步
		}
		else
		{
			Sequence[Site] = 1;
		}
	}
	// 在 Norbs 个空间轨道中随机选择 Beta 个，填充 beta 电子
	for(int ix = 0; ix != Beta; ++ix)
	{
		Site = 2 * randSite(0, Norbs - 1) + 1; // beta 轨道在后，其索引为奇数
		if(Sequence[Site] == 1) // 如果当前位置已经是 1，说明此位置已经填过电子，那么当前循环无效
		{
			--ix; // 使当前循环无效，需要将循环计数器退一步
		}
		else
		{
			Sequence[Site] = 1;
		}
	}
	return *this; // 返回引用是为了方便之后的赋值操作
}

// 按照给定的概率发生变异
// 变异策略：对每个位置进行变异操作，是否发生变异由 StepSize 控制
Configuration& Configuration::mutate(double rate, double StepSize)
{
	Random<double> rand; // 随机浮点数生成器
	// 如果随机数大于或者等于变异几率，则不发生变异
	if(rand(0, 1) < rate) // 如果要发生变异
	{
		// 对每个 alpha 和 beta 轨道进行操作
		for(int i = 0; i < SpinOrbs; ++i)
		{
			for(int j = i + 2; j < SpinOrbs; j += 2)
			{
				if(rand(0, 1) < StepSize) std::swap(Sequence[i], Sequence[j]);
			}
		}
	}
	return *this;
}

// 变异策略2，MonteCarlo 程序文献中的方法
Configuration& Configuration::MutateStrategy2()
{
	Random<int> rand;
	int excitations = rand(0, Nelec); // 确定有多少个激发
	int Site1; // 激发算符的作用位置
	int Site2; // 湮灭算符的作用位置
	// 对于每个激发操作，组成一个激发算符
	for(int idx_excit = 0; idx_excit != excitations; ++idx_excit)
	{
		// 在 Sequence 序列中随机选择电子产生位置与消除位置
		Site1 = rand(0, SpinOrbs - 1);
		do { // 保证电子产生位置与电子消除位置不同，且操作的电子为相同自旋
			Site2 = 2 * rand(0, SpinOrbs / 2 - 1) + (Site1 % 2);
		} while(Site1 == Site2);
		// 激发操作：互换两个位置的状态
		std::swap(Sequence[Site1], Sequence[Site2]);
	}
	return *this;
}

Configuration& Configuration::operator=(const Configuration& rhs)
{
	if(this != &rhs)
	{
		delete[] Sequence;
		SpinOrbs = rhs.SpinOrbs;
		Nelec = rhs.Nelec;
		Spin = rhs.Spin;
		Sequence = new int[SpinOrbs];
		Coefficient = rhs.Coefficient;
		for(int ix = 0; ix != SpinOrbs; ++ix)
		{
			Sequence[ix] = rhs.Sequence[ix];
		}
	}
	return *this;
}

// 根据序列 Sequence 计算电子个数和自旋多重度
void Configuration::ComputeNelecSpin()
{
	Nelec = 0;
	int Nalpha = 0;
	for(int ix = 0; ix != SpinOrbs; ++ix)
	{
		Nelec += Sequence[ix];
		if(ix % 2 == 0 && Sequence[ix] == 1)
		{
			++Nalpha;
		}
	}
	Spin = std::abs(2 * Nalpha - Nelec) + 1;
}

// 从指定空间轨道位置 site 开始，按照 status（1,2,3,4）连续填充 Sequence 中的两个位置
void Configuration::SetOccupationStatus(int site, int status)
{
	int AlphaSite = 2 * site;
	int BetaSite = AlphaSite + 1;
	switch(status)
	{
	case 1: // 无占据
		Sequence[AlphaSite] = 0;
		Sequence[BetaSite] = 0;
		break;
	case 2: // 只有 beta 占据
		Sequence[AlphaSite] = 0;
		Sequence[BetaSite] = 1;
		break;
	case 3: // 只有 alpha 占据
		Sequence[AlphaSite] = 1;
		Sequence[BetaSite] = 0;
		break;
	case 4: // 双占据
		Sequence[AlphaSite] = 1;
		Sequence[BetaSite] = 1;
		break;
	default: // status 出错的情况，抛出异常
		throw std::runtime_error("ERROR: invalid occupation status!");
	}
}

bool operator==(const Configuration& lhs, const Configuration& rhs)
{
#ifdef DEBUG
	if(lhs.SpinOrbs != rhs.SpinOrbs || lhs.Nelec != rhs.Nelec || lhs.Spin != rhs.Spin)
	{
		throw std::runtime_error("Configuration sequence error!");
	}
#endif // DEBUG
	if(lhs.Coefficient != rhs.Coefficient)
	{
		return false;
	}
	else
	{
		for(int ix = 0; ix != lhs.SpinOrbs; ++ix)
		{
			if(*(lhs.Sequence + ix) != *(rhs.Sequence + ix)) return false;
		}
	}
	return true;
}

std::ostream& operator<<(std::ostream& os, const Configuration& rhs)
{
	for(int ix = 0; ix < rhs.SpinOrbs; ix += 2)
	{
		os << 2 * rhs.Sequence[ix] + rhs.Sequence[ix + 1] + 1;
	}
	return os;
}

bool SortByCoefficient::operator()(const Configuration& lhs, const Configuration& rhs)
{
	return (std::abs(lhs.Coefficient) > std::abs(rhs.Coefficient));
}

bool SortBySequence::operator()(const Configuration& lhs, const Configuration& rhs)
{
#ifdef DEBUG
	if(lhs.SpinOrbs != rhs.SpinOrbs)
	{
		throw std::runtime_error("Configuration sequence error!");
	}
#endif // DEBUG
	for(int ix = 0; ix != lhs.SpinOrbs; ++ix)
	{
		if(*(lhs.Sequence + ix) > *(rhs.Sequence + ix)) return true;
	}
	return false;
}