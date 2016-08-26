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
	Generation = 0; // ��ʼ��Ⱥ����Ϊ�� 0 ��
	if(Nelec % 2 == Spin % 2) // ����������ض��Ƿ���ȷ
	{
		// std::cerr �����ڲ��ԡ�Ӧ��ʹ���쳣����ķ�ʽ��
		std::cerr << "input error: spin multiplicity invalid!" << std::endl;
		throw std::runtime_error("input error: spin multiplicity invalid!");
	}
}

// �ݻ� gen ��Ŀ���Ӵ����������Ϊ MutateRate��os Ϊ���������Ϣ��Ŀ��
// ÿһ��ѭ���У�ϵ������ֵ��С�� CutOff ��ֵ����̬���ᱻ��¼
// ��̬��¼���ļ���Ϊ LargestCoeffFile (Ĭ��Ϊ "largecoeffs.txt")
void StateFunctions::evolute(int gen, double MutateRate,
	double CutOff, const std::string& LargestCoeffFile, std::ostream& os)
{
	// ���ڼ�¼ÿ��ѭ����������̬ϵ�����ڸ�����ֵ����̬
	std::ofstream LargeCoeffsRecord(LargestCoeffFile.c_str());
	std::vector<Configuration> Descendants; // �������ڱ����Ӵ�������
	double CurrentCoeffSquare;
	os << "Generation\tLargestCoeff\tSumCoeffSquare" << std::endl;
	for(std::size_t GenIdx = 0; GenIdx <= gen; ++GenIdx)
	{
		ResetCoeff(); // ���� LargestCoeff, SumCoeffSquare �� AverageCoeffSquare ����ֵ
		ComputeCoefficients(); // ����������̬��ϵ��
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
		std::vector<Configuration>().swap(Descendants); // ������ڱ����Ӵ�������
		for(int ix = 0; ix != Population; ++ix) // �ӱ��������ȡ������б��죬�������Ӵ�
		{
			Descendants.push_back(Configuration(RandomChoose()).MutateStrategy2());
		}
		CSFs.swap(Descendants); // ���Ӵ�����Ϊ������������һ��ѭ��
		++Generation;
	}
	// ѭ��������Descendants ���������һ������ CSFs �������һ�����Ӵ������Ӵ��Ѿ����ٽ�����һ��ѭ����
	// ������Ҫ�� CSFs ���������±���Ϊ���һ��
	CSFs.swap(Descendants);
	// �ر��ļ����
	LargeCoeffsRecord.close();
}

// ��� count �� ���������ʽ
void StateFunctions::FillRandomDets(int count)
{
	Configuration csf(Norbs, Nelec, Spin);
	for(int ix = 0; ix < count; ++ix)
	{
		CSFs.push_back(csf.randomize());
	}
	Population = CSFs.size(); // ���¼��㵱ǰ�ĸ�����Ŀ
}

// ��� count �� Hartree-Fock ���ռ������ʽ
void StateFunctions::FillHartreeFockDets(int count)
{
	Configuration csf(Norbs, Nelec, Spin);
	csf.SetHartreeFockDet();
	for(int ix = 0; ix < count; ++ix)
	{
		CSFs.push_back(csf);
	}
	Population = CSFs.size(); // ���¼��㵱ǰ�ĸ�����Ŀ
}

// ���Թ��ܣ�CAS(6,6) ����չ
int StateFunctions::SmallCASExpand(int CoreOrbs, int CoreElec, int CoreSpin)
{
	Array SpinOrbitals(CoreOrbs);
	int count = 0;
	int NCoreElecAlpha = (CoreElec + CoreSpin - 1) / 2;
	std::vector<Array> CASCoreAlpha;
	SpinOrbitals.initialize();
	if(CoreSpin == 1) //  CoreSpin == 1 ʱ
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
	else // CoreSpin != 1 ʱ
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

// ���Թ��ܣ���С�ռ�CAS��ٵ��������һЩ������
int StateFunctions::TEST(int CoreOrbs, int CoreElec, int CoreSpin, int OuterOrbs, int OuterElec, int OuterSpin)
{
	Array SpinOrbitals(CoreOrbs);
	int count = 0;
	int NCoreElecAlpha = (CoreElec + CoreSpin - 1) / 2;
	int NOuterElecAlpha = (OuterElec + OuterSpin - 1) / 2;
	int NOuterElecBeta = (OuterElec - OuterSpin + 1) / 2;
	std::vector<Array> CASCoreAlpha;
	SpinOrbitals.initialize();
	if(CoreSpin == 1) //  CoreSpin == 1 ʱ
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
	else // CoreSpin != 1 ʱ
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

// ����� 444111 �������й��ɵ� count ����̬
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
	Population = CSFs.size(); // ���¼��㵱ǰ�ĸ�����Ŀ
}

// ��ȡ��ǰ��Ⱥ�ĸ�����
int StateFunctions::population() const
{
	return Population;
}

// ���� checkpoint �ļ��е�����
void StateFunctions::CheckPointFolder(const std::string & ChkName)
{
	ChkPointFolderName = ChkName;
}

// �� LargestCoeff, SumCoeffSquare �� AverageCoeffSquare ����
void StateFunctions::ResetCoeff()
{
	LargestCoeff = 0;
	SumCoeffSquare = 0;
}

// ����������̬��ϵ��
void StateFunctions::ComputeCoefficients()
{
	// �������ڵ�������̬���� QCmaquis ��ʽд��һ�� dets_list �ļ���
	std::ofstream DetsList("dets_list");
	for(std::vector<Configuration>::iterator it = CSFs.begin(); it != CSFs.end(); ++it)
	{
		DetsList << *it << std::endl;
	}
	DetsList.close(); // dets_list ������ϣ�ǿ��ˢ�»�����
					  // ���� QCmaquis ���㱾��������̬��ϵ��
	system("mps2ci_2u1pg checkpoint_state.0.0.0.h5 dets_list > dets_result");
	std::ifstream DetsCoeffs("dets_result");
	if(!DetsCoeffs) // ����Ƿ����������ļ�
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
	DetsCoeffs.close(); // coefficients ��ȡ���
}

// �������̶ķ������ѡ��һ����̬
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

// �Ա����������̬����ϵ������ֵ�Ӵ�С��˳���������
void StateFunctions::sort()
{
	std::sort(CSFs.begin(), CSFs.end(), SortByCoefficient());
}

// ���������̬������ Population �� Generation ����
void StateFunctions::clear()
{
	std::vector<Configuration>().swap(CSFs);
	Population = 0;
	Generation = 0;
}

// ���ص����������
std::ostream& operator<<(std::ostream& os, const StateFunctions& rhs)
{
	// ���������̬������ QCmaquis �ĸ�ʽ����ϵ��
	for(std::vector<Configuration>::const_iterator it = rhs.CSFs.begin(); it != rhs.CSFs.end(); ++it)
	{
		os << *it << "\t" << it->coeff() << std::endl;
	}
	return os;
}