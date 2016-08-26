#include "StateFunctions.h"
#include <stdexcept>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <stdlib.h>

StateFunctions::StateFunctions(int orbs, int elec, int spin)
{
#ifdef _INPUT_CHECK
	if(orbs < 0 || elec < 0 || spin < 0)
	{
		std::cerr << "input error: negative numbers unacceptable!" << std::endl;
		throw std::runtime_error("input error: negative numbers unacceptable!");
	}
#endif // _INPUT_CHECK
	Population = 0;
	Norbs = orbs;
	Nelec = elec;
	Spin = spin;
	Generation = 0; // 初始种群，记为第 0 代
	if(Nelec % 2 == Spin % 2) // 检查自旋多重度是否正确
	{
		// std::cerr 仅用于测试。应该使用异常处理的方式。
		std::cerr << "input error: spin multiplicity invalid!" << std::endl;
		throw std::runtime_error("input error: spin multiplicity invalid!");
	}
}

// 演化 gen 数目的子代，变异概率为 MutateRate，os 为输出过程信息的目标
// 每一次循环中，系数绝对值不小于 CutOff 阈值的组态都会被记录
// 组态记录的文件名为 LargestCoeffFile (默认为 "largecoeffs.txt")
void StateFunctions::evolute(int gen, double MutateRate,
	double CutOff, const std::string& LargestCoeffFile, std::ostream& os)
{
	// 用于记录每次循环过程中组态系数大于给定阈值的组态
	std::ofstream LargeCoeffsRecord(LargestCoeffFile.c_str());
	std::vector<Configuration> Descendants; // 定于用于保存子代的数组
	double CurrentCoeffSquare;
	os << "Generation\tLargestCoeff\tSumCoeffSquare" << std::endl;
	for(std::size_t GenIdx = 0; GenIdx <= gen; ++GenIdx)
	{
		ResetCoeff(); // 重置 LargestCoeff, SumCoeffSquare 和 AverageCoeffSquare 的数值
		ComputeCoefficients(); // 计算所有组态的系数
		for(std::vector<Configuration>::iterator it = CSFs.begin(); it != CSFs.end(); ++it)
		{
			CurrentCoeffSquare = std::pow(it->Coefficient, 2);
			if(CurrentCoeffSquare > std::pow(LargestCoeff, 2))
			{
				LargestCoeff = it->Coefficient;
			}
			if(CurrentCoeffSquare > CutOff)
			{
				LargeCoeffsRecord << *it << "\t" << it->coeff() << std::endl;
			}
			SumCoeffSquare += CurrentCoeffSquare;
		}
		os << Generation << "\t" << LargestCoeff << "\t" << SumCoeffSquare << std::endl;
		std::vector<Configuration>().swap(Descendants); // 清空用于保存子代的数组
		for(int ix = 0; ix != Population; ++ix) // 从本带随机抽取个体进行变异，并存入子代
		{
			Descendants.push_back(Configuration(RandomChoose()).MutateStrategy2());
		}
		CSFs.swap(Descendants); // 将子代保存为本代，进行下一轮循环
		++Generation;
	}
	// 循环结束后，Descendants 保留着最后一代，而 CSFs 则是最后一代的子代（此子代已经不再进行下一次循环）
	// 所以需要将 CSFs 即本代重新保存为最后一代
	CSFs.swap(Descendants);
	// 关闭文件句柄
	LargeCoeffsRecord.close();
}

// 填充 count 个 随机的行列式
void StateFunctions::FillRandomDets(int count)
{
	Configuration csf(Norbs, Nelec, Spin);
	for(int ix = 0; ix < count; ++ix)
	{
		CSFs.push_back(csf.randomize());
	}
	Population = CSFs.size(); // 重新计算当前的个体数目
}

// 填充 count 个 Hartree-Fock 最低占据行列式
void StateFunctions::FillHartreeFockDets(int count)
{
	Configuration csf(Norbs, Nelec, Spin);
	csf.SetHartreeFockDet();
	for(int ix = 0; ix < count; ++ix)
	{
		CSFs.push_back(csf);
	}
	Population = CSFs.size(); // 重新计算当前的个体数目
}

// 测试功能，CAS(6,6) 的扩展
int StateFunctions::SmallCASExpand(int CoreOrbs, int CoreElec, int CoreSpin)
{
	Array SpinOrbitals(CoreOrbs);
	int count = 0;
	int NCoreElecAlpha = (CoreElec + CoreSpin - 1) / 2;
	std::vector<Array> CASCoreAlpha;
	SpinOrbitals.initialize();
	if(CoreSpin == 1) //  CoreSpin == 1 时
	{
		SpinOrbitals.DFS(CASCoreAlpha, CoreElec / 2, CoreOrbs);
		for(std::vector<Array>::iterator i = CASCoreAlpha.begin(); i != CASCoreAlpha.end(); ++i)
		{
			for(std::vector<Array>::iterator j = CASCoreAlpha.begin(); j != CASCoreAlpha.end(); ++j)
			{
				CSFs.push_back(Configuration(*i, *j, Norbs, Nelec, Spin, CoreOrbs, NCoreElecAlpha, NCoreElecAlpha));
				++count;
			}
		}
	}
	else // CoreSpin != 1 时
	{
		SpinOrbitals.DFS(CASCoreAlpha, NCoreElecAlpha, CoreOrbs);
		int NCoreElecBeta = (CoreElec - CoreSpin + 1) / 2;
		std::vector<Array> CASCoreBeta;
		SpinOrbitals.initialize();
		SpinOrbitals.DFS(CASCoreBeta, NCoreElecBeta, CoreOrbs);
		for(std::vector<Array>::iterator i = CASCoreAlpha.begin(); i != CASCoreAlpha.end(); ++i)
		{
			for(std::vector<Array>::iterator j = CASCoreBeta.begin(); j != CASCoreBeta.end(); ++j)
			{
				CSFs.push_back(Configuration(*i, *j, Norbs, Nelec, Spin, CoreOrbs, NCoreElecAlpha, NCoreElecBeta));
				++count;
			}
		}
	}
	Population += count;
	return count;
}

// 测试功能：在小空间CAS穷举的外面进行一些随机填充
int StateFunctions::TEST(int CoreOrbs, int CoreElec, int CoreSpin, int OuterOrbs, int OuterElec, int OuterSpin)
{
	Array SpinOrbitals(CoreOrbs);
	int count = 0;
	int NCoreElecAlpha = (CoreElec + CoreSpin - 1) / 2;
	int NOuterElecAlpha = (OuterElec + OuterSpin - 1) / 2;
	int NOuterElecBeta = (OuterElec - OuterSpin + 1) / 2;
	std::vector<Array> CASCoreAlpha;
	SpinOrbitals.initialize();
	if(CoreSpin == 1) //  CoreSpin == 1 时
	{
		SpinOrbitals.DFS(CASCoreAlpha, CoreElec / 2, CoreOrbs);
		for(std::vector<Array>::iterator i = CASCoreAlpha.begin(); i != CASCoreAlpha.end(); ++i)
		{
			for(std::vector<Array>::iterator j = CASCoreAlpha.begin(); j != CASCoreAlpha.end(); ++j)
			{
				CSFs.push_back(Configuration(*i, *j, Norbs, Nelec, Spin,
					CoreOrbs, NCoreElecAlpha, NCoreElecAlpha,
					OuterOrbs, NOuterElecAlpha, NOuterElecBeta));
				++count;
			}
		}
	}
	else // CoreSpin != 1 时
	{
		SpinOrbitals.DFS(CASCoreAlpha, NCoreElecAlpha, CoreOrbs);
		int NCoreElecBeta = (CoreElec - CoreSpin + 1) / 2;
		std::vector<Array> CASCoreBeta;
		SpinOrbitals.initialize();
		SpinOrbitals.DFS(CASCoreBeta, NCoreElecBeta, CoreOrbs);
		for(std::vector<Array>::iterator i = CASCoreAlpha.begin(); i != CASCoreAlpha.end(); ++i)
		{
			for(std::vector<Array>::iterator j = CASCoreBeta.begin(); j != CASCoreBeta.end(); ++j)
			{
				CSFs.push_back(Configuration(*i, *j, Norbs, Nelec, Spin,
					CoreOrbs, NCoreElecAlpha, NCoreElecBeta,
					OuterOrbs, NOuterElecAlpha, NOuterElecBeta));
				++count;
			}
		}
	}
	Population += count;
	return count;
}

// 填充由 444111 这种序列构成的 count 个组态
void StateFunctions::FillWithSeed(const std::string & seq, int count)
{
	Configuration csf(seq);
#ifdef _INPUT_CHECK
	if(csf.SpinOrbs != 2 * Norbs || csf.Nelec != Nelec || csf.Spin != Spin)
	{
		throw std::runtime_error("ERROR: the given seed invalid!");
	}
#endif // _INPUT_CHECK
	for(int ix = 0; ix < count; ++ix)
	{
		CSFs.push_back(csf);
	}
	Population = CSFs.size(); // 重新计算当前的个体数目
}

// 获取当前种群的个体数
int StateFunctions::population() const
{
	return Population;
}

// 设置 checkpoint 文件夹的名称
void StateFunctions::CheckPointFolder(const std::string & ChkName)
{
	ChkPointFolderName = ChkName;
}

// 将 LargestCoeff, SumCoeffSquare 和 AverageCoeffSquare 重置
void StateFunctions::ResetCoeff()
{
	LargestCoeff = 0;
	SumCoeffSquare = 0;
}

// 计算所有组态的系数
void StateFunctions::ComputeCoefficients()
{
	// 将本组内的所有组态按照 QCmaquis 格式写到一个 dets_list 文件中
	std::ofstream DetsList("dets_list");
	for(std::vector<Configuration>::iterator it = CSFs.begin(); it != CSFs.end(); ++it)
	{
		DetsList << *it << std::endl;
	}
	DetsList.close(); // dets_list 生成完毕，强制刷新缓冲区
					  // 调用 QCmaquis 计算本组所有组态的系数
	system("mps2ci_2u1pg checkpoint_state.0.0.0.h5 dets_list > dets_result");
	std::ifstream DetsCoeffs("dets_result");
	if(!DetsCoeffs) // 检测是否正常打开了文件
	{
		throw std::runtime_error("coefficients computation error: failed to open file dets_result");
	}
	std::string Coeff;
	// skip useless lines
	for(int ix = 0; ix != Population; ++ix) getline(DetsCoeffs, Coeff);
	// read coefficients
	for(int ix = 0; ix != Population; ++ix)
	{
		getline(DetsCoeffs, Coeff);
		CSFs[ix].Coefficient = atof(Coeff.substr(Coeff.rfind(' ') + 1).c_str());
	}
	DetsCoeffs.close(); // coefficients 读取完毕
}

// 按照轮盘赌方法随机选择一个组态
const Configuration& StateFunctions::RandomChoose()
{
	Random<double> rand;
	double Slice = rand(0, 1) * SumCoeffSquare;
	double CoeffSquareSoFar = 0;
	for(std::vector<Configuration>::iterator it = CSFs.begin(); it != CSFs.end(); ++it)
	{
		CoeffSquareSoFar += std::pow(it->Coefficient, 2);
		if(CoeffSquareSoFar >= Slice) return *it;
	}
}

// 对本组的所有组态按照系数绝对值从大到小的顺序进行排序
void StateFunctions::sort()
{
	std::sort(CSFs.begin(), CSFs.end(), SortByCoefficient());
}

// 清除所有组态，并将 Population 和 Generation 置零
void StateFunctions::clear()
{
	std::vector<Configuration>().swap(CSFs);
	Population = 0;
	Generation = 0;
}

// 重载的输出操作符
std::ostream& operator<<(std::ostream& os, const StateFunctions& rhs)
{
	// 逐行输出组态（按照 QCmaquis 的格式）与系数
	for(std::vector<Configuration>::const_iterator it = rhs.CSFs.begin(); it != rhs.CSFs.end(); ++it)
	{
		os << *it << "\t" << it->coeff() << std::endl;
	}
	return os;
}