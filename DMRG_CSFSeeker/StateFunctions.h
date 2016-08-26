#ifndef STATE_FUNCTIONS_CLASS_H
#define STATE_FUNCTIONS_CLASS_H

#include "Array.h"
#include "Configuration.h"
#include <iostream>
#include <vector>
#include <string>

class StateFunctions
{
public:
	StateFunctions(int orbs, int elec, int spin = 1);
	// 演化 gen 数目的子代，变异概率为 MutateRate，os 为输出过程信息的目标
	void evolute(int gen, double MutateRate, double CutOff = 1E-8,
		const std::string& LargestCoeffFile = "largecoeffs.txt", std::ostream& os = std::cout);
	void FillRandomDets(int count);  // 填充 count 个 随机的行列式
	void FillHartreeFockDets(int count); // 填充 count 个 Hartree-Fock 最低占据行列式
	int SmallCASExpand(int CoreOrbs = 6, int CoreElec = 6, int CoreSpin = 1); // 填充 count 个行列式，其中一个是最低占据，另外几个依次升高。
	int TEST(int CoreOrbs, int CoreElec, int CoreSpin, int OuterOrbs, int OuterElec, int OuterSpin); // 测试功能，在CAS核心扩展之外增加一些随机填充
	void FillWithSeed(const std::string& seq, int count = 1); // 填充由 444111 这种序列构成的 count 个组态
	int population() const; // 获取当前种群的个体数
	void CheckPointFolder(const std::string& ChkName); // 设置 checkpoint 文件夹的名称
	void sort(); // 对本组的所有组态按照系数绝对值从大到小的顺序进行排序
	void clear(); // 清除所有组态，并将 Population 和 Generation 置零
	// 重载的输出操作符
	friend std::ostream& operator<<(std::ostream& os, const StateFunctions& rhs);
private:
	int Population; // CSFs 的个数
	int Norbs; // 轨道（空间）个数
	int Nelec; // 电子数
	int Spin; // 自旋多重度 S=2s+1
	int Generation; // 当前的代数(初代为0)
	std::string ChkPointFolderName{ "checkpoint_state.0.0.0.h5" }; // checkpoint 文件夹的名称
	double LargestCoeff; // 本组内的绝对值最大系数
	double SumCoeffSquare; // 本组的系数平方之和
	std::vector<Configuration> CSFs; // 组态的集合
	// 私有辅助函数
	void ResetCoeff(); // 将 LargestCoeff, SumCoeffSquare 和 AverageCoeffSquare 重置
	void ComputeCoefficients(); // 为每个组态计算系数。在此处集中计算而非对每个组态分别计算是为了性能考虑
	const Configuration& RandomChoose(); // 按照轮盘赌方法随机选择一个组态
	//  void DFS_Exhaustion(Configuration& csf, int elec, int LastSite); // 用于 CAS(6,6) 向外扩展取样方法的深度优先遍历算法
};

#endif // STATE_FUNCTIONS_CLASS_H