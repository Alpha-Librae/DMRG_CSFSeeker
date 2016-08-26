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
	// �ݻ� gen ��Ŀ���Ӵ����������Ϊ MutateRate��os Ϊ���������Ϣ��Ŀ��
	void evolute(int gen, double MutateRate, double CutOff = 1E-8,
		const std::string& LargestCoeffFile = "largecoeffs.txt", std::ostream& os = std::cout);
	void FillRandomDets(int count);  // ��� count �� ���������ʽ
	void FillHartreeFockDets(int count); // ��� count �� Hartree-Fock ���ռ������ʽ
	int SmallCASExpand(int CoreOrbs = 6, int CoreElec = 6, int CoreSpin = 1); // ��� count ������ʽ������һ�������ռ�ݣ����⼸���������ߡ�
	int TEST(int CoreOrbs, int CoreElec, int CoreSpin, int OuterOrbs, int OuterElec, int OuterSpin); // ���Թ��ܣ���CAS������չ֮������һЩ������
	void FillWithSeed(const std::string& seq, int count = 1); // ����� 444111 �������й��ɵ� count ����̬
	int population() const; // ��ȡ��ǰ��Ⱥ�ĸ�����
	void CheckPointFolder(const std::string& ChkName); // ���� checkpoint �ļ��е�����
	void sort(); // �Ա����������̬����ϵ������ֵ�Ӵ�С��˳���������
	void clear(); // ���������̬������ Population �� Generation ����
	// ���ص����������
	friend std::ostream& operator<<(std::ostream& os, const StateFunctions& rhs);
private:
	int Population; // CSFs �ĸ���
	int Norbs; // ������ռ䣩����
	int Nelec; // ������
	int Spin; // �������ض� S=2s+1
	int Generation; // ��ǰ�Ĵ���(����Ϊ0)
	std::string ChkPointFolderName{ "checkpoint_state.0.0.0.h5" }; // checkpoint �ļ��е�����
	double LargestCoeff; // �����ڵľ���ֵ���ϵ��
	double SumCoeffSquare; // �����ϵ��ƽ��֮��
	std::vector<Configuration> CSFs; // ��̬�ļ���
	// ˽�и�������
	void ResetCoeff(); // �� LargestCoeff, SumCoeffSquare �� AverageCoeffSquare ����
	void ComputeCoefficients(); // Ϊÿ����̬����ϵ�����ڴ˴����м�����Ƕ�ÿ����̬�ֱ������Ϊ�����ܿ���
	const Configuration& RandomChoose(); // �������̶ķ������ѡ��һ����̬
	//  void DFS_Exhaustion(Configuration& csf, int elec, int LastSite); // ���� CAS(6,6) ������չȡ��������������ȱ����㷨
};

#endif // STATE_FUNCTIONS_CLASS_H