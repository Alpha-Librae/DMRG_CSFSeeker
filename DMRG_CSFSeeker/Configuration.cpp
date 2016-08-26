#include "Configuration.h"
#include <cmath>
#include <stdexcept>

// constructors
// ���캯�������ܵĲ���Ϊ���Կռ�Ĺ�������������Լ���̬���������ض�(Ĭ��Ϊ 1��������̬)
Configuration::Configuration(int orbs, int elec, int spin)
{
	SpinOrbs = 2 * orbs; // orbs Ϊ�ռ�������Ŀ
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

	int OuterAlpha = (Nelec + Spin - 1) / 2 - CASCoreElecAlpha; // Core �ռ��� alpha ���ӵ���Ŀ
	int OuterBeta = (Nelec - Spin + 1) / 2 - CASCoreElecBeta; // Core �ռ��� beta ���ӵ���Ŀ

	initialize(); // �Ƚ����е�λ�ö����Ϊû�е���
	// �������� CAS �ռ�����ȫռ�ݲ���
	int idx = 0;
	while(idx < OuterAlpha)
	{
		alpha(idx) = 1;
		++idx;
	}
	// Ȼ���� CoreAlpha �������ʣ�µ� alpha ����
	for(int i = 0; i != CoreAlpha.size(); ++i)
	{
		alpha(idx) = CoreAlpha[i];
		++idx;
	}
	// ���� beta ���ӣ�ִ�����ƵĲ���
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

// ���Թ���
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

	int NElecAlpha = (Nelec + Spin - 1) / 2; // �ܵ� alpha ������
	int NElecBeta = (Nelec - Spin + 1) / 2; // �ܵ� beta ������
	int NREAlpha = NLargerSpaceElecAlpha - CASCoreElecAlpha; // �ڴ�С�ռ��ڽ���������� alpha ������
	int NREBeta = NLargerSpaceElecBeta - CASCoreElecBeta; // �ڴ�С�ռ��ڽ���������� beta ������
	int NFEAlpha = NElecAlpha - NLargerSpaceElecAlpha; // С�ռ����̶����λ�õ� alpha ������
	int NFEBeta = NElecBeta - NLargerSpaceElecBeta; // С�ռ����̶����λ�õ� beta ������
	Random<int> rand; // ȷ��������λ�õ������

	initialize(); // �Ƚ����еĹ������ 0

	// alpha ����
	// �����ȫռ�ݲ���
	int idx = 0;
	while(idx < NFEAlpha)
	{
		alpha(idx) = 1;
		++idx;
	}
	// ���������벿��
	int NRFSites = LargerSpaceOrbs - CASCoreOrbs; // ������Ĺ������Ŀ
	int site = 0; // ��¼�������λ��
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
	// С�ռ� CAS ��ٲ���
	for(int i = 0; i < CoreAlpha.size(); ++i)
	{
		alpha(idx) = CoreAlpha[i];
		++idx;
	}
	// ��������Ұ벿��
	int NRFOrbsRight = (LargerSpaceOrbs - (LargerSpaceOrbs % 2)) / 2;
	int NCASOrbsRight = (CASCoreOrbs - (CASCoreOrbs % 2)) / 2;
	int NRFOrbsU = NRFOrbsRight - NCASOrbsRight;
	for(int i = 0; i < NRFOrbsU; ++i)
	{
		alpha(idx) = RFOrbitals[site];
		++site;
		++idx;
	}

	// beta ����
	// �����ȫռ�ݲ���
	idx = 0;
	while(idx < NFEBeta)
	{
		beta(idx) = 1;
		++idx;
	}
	// ���������벿��
	site = 0; // ��¼�������λ��
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
	// С�ռ� CAS ��ٲ���
	for(int i = 0; i < CoreBeta.size(); ++i)
	{
		beta(idx) = CoreBeta[i];
		++idx;
	}
	// ��������Ұ벿��
	for(int i = 0; i < NRFOrbsU; ++i)
	{
		beta(idx) = RFOrbitals[site];
		++site;
		++idx;
	}
}

// �� 444111 �������й���һ����̬
Configuration::Configuration(const std::string & rhs)
{
	SpinOrbs = 2 * rhs.length();
	Sequence = new int[SpinOrbs];
	for(int ix = 0; ix != rhs.length(); ++ix)
	{
		SetOccupationStatus(ix, rhs[ix] - '0');
	}
	// ������Ӹ������������ض�
	ComputeNelecSpin();
}

// the deconstructor
Configuration::~Configuration()
{
	delete[] Sequence;
}

// ��������Ϊ index �� alpha ���������״̬
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

// ��ֵ��������Ϊ index �� alpha ���������״̬
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

// ��������Ϊ index �� beta ���������״̬
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

// ��ֵ��������Ϊ index �� beta ���������״̬
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

// ��������Ϊ index �����������״̬
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

// ��ֵ��������Ϊ index �����������״̬
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
// ������̬ϵ��
double Configuration::coeff() const
{
	return Coefficient;
}

// �� Sequence ��ÿһλ���� 0
void Configuration::initialize()
{
	for(int ix = 0; ix != SpinOrbs; ++ix)
	{
		Sequence[ix] = 0;
	}
}

// ����̬�趨Ϊ Hartree-Fock ����ʽ��ʾ�����ռ����ʽ
void Configuration::SetHartreeFockDet()
{
	int Alpha = (Nelec + Spin - 1) / 2; // alpha ���ӵ���Ŀ
	int Beta = (Nelec - Spin + 1) / 2; // beta ���ӵ���Ŀ
	initialize(); // �Ƚ�ȫ����λ�ñ��Ϊû�е���
	for(int ix = 0; ix != Alpha; ++ix) // Ȼ����� alpha ����
	{
		Sequence[2 * ix] = 1;
	}
	for(int ix = 0; ix != Beta; ++ix) // ������ beta ����
	{
		Sequence[2 * ix + 1] = 1;
	}
}

// �� 444111 �����ַ������й���һ����̬
void Configuration::SetBySeed(const std::string & rhs)
{
	int Norbs = SpinOrbs / 2;
	int Length;
	if(Norbs <= rhs.length()) // ����ַ�����̫��������Գ����Ĳ���
	{
		Length = Norbs;
	}
	else // ����ַ�����̫�̣����ȸ�������û�������Ĳ��ֲ�0����ʾû�е���
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
	// ������Ӹ������������ض�
	ComputeNelecSpin();
}

const Configuration& Configuration::randomize()
{
	Random<int> randSite; // ���������������
	int Alpha = (Nelec + Spin - 1) / 2; // alpha ���ӵ���Ŀ
	int Beta = (Nelec - Spin + 1) / 2; // beta ���ӵ���Ŀ
	int Norbs = SpinOrbs / 2; // �ռ�������Ŀ
	int Site; // ��¼�������λ��
			  // �洢�������� ababab... ���� alpha �� beta ���������
	initialize(); // ���Ƚ����г�ʼ��Ϊȫ 0����ʾ���е�λ�ö�û�е���
				  // �� Norbs ���ռ��������ѡ�� Alpha ������� alpha ����
	for(int ix = 0; ix != Alpha; ++ix)
	{
		Site = 2 * randSite(0, Norbs - 1); // alpha �����ǰ��������Ϊż��
		if(Sequence[Site] == 1) // �����ǰλ���Ѿ��� 1��˵����λ���Ѿ�������ӣ���ô��ǰѭ����Ч
		{
			--ix; // ʹ��ǰѭ����Ч����Ҫ��ѭ����������һ��
		}
		else
		{
			Sequence[Site] = 1;
		}
	}
	// �� Norbs ���ռ��������ѡ�� Beta ������� beta ����
	for(int ix = 0; ix != Beta; ++ix)
	{
		Site = 2 * randSite(0, Norbs - 1) + 1; // beta ����ں�������Ϊ����
		if(Sequence[Site] == 1) // �����ǰλ���Ѿ��� 1��˵����λ���Ѿ�������ӣ���ô��ǰѭ����Ч
		{
			--ix; // ʹ��ǰѭ����Ч����Ҫ��ѭ����������һ��
		}
		else
		{
			Sequence[Site] = 1;
		}
	}
	return *this; // ����������Ϊ�˷���֮��ĸ�ֵ����
}

// ���ո����ĸ��ʷ�������
// ������ԣ���ÿ��λ�ý��б���������Ƿ��������� StepSize ����
Configuration& Configuration::mutate(double rate, double StepSize)
{
	Random<double> rand; // ���������������
	// �����������ڻ��ߵ��ڱ��켸�ʣ��򲻷�������
	if(rand(0, 1) < rate) // ���Ҫ��������
	{
		// ��ÿ�� alpha �� beta ������в���
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

// �������2��MonteCarlo ���������еķ���
Configuration& Configuration::MutateStrategy2()
{
	Random<int> rand;
	int excitations = rand(0, Nelec); // ȷ���ж��ٸ�����
	int Site1; // �������������λ��
	int Site2; // �������������λ��
	// ����ÿ���������������һ���������
	for(int idx_excit = 0; idx_excit != excitations; ++idx_excit)
	{
		// �� Sequence ���������ѡ����Ӳ���λ��������λ��
		Site1 = rand(0, SpinOrbs - 1);
		do { // ��֤���Ӳ���λ�����������λ�ò�ͬ���Ҳ����ĵ���Ϊ��ͬ����
			Site2 = 2 * rand(0, SpinOrbs / 2 - 1) + (Site1 % 2);
		} while(Site1 == Site2);
		// ������������������λ�õ�״̬
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

// �������� Sequence ������Ӹ������������ض�
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

// ��ָ���ռ���λ�� site ��ʼ������ status��1,2,3,4��������� Sequence �е�����λ��
void Configuration::SetOccupationStatus(int site, int status)
{
	int AlphaSite = 2 * site;
	int BetaSite = AlphaSite + 1;
	switch(status)
	{
	case 1: // ��ռ��
		Sequence[AlphaSite] = 0;
		Sequence[BetaSite] = 0;
		break;
	case 2: // ֻ�� beta ռ��
		Sequence[AlphaSite] = 0;
		Sequence[BetaSite] = 1;
		break;
	case 3: // ֻ�� alpha ռ��
		Sequence[AlphaSite] = 1;
		Sequence[BetaSite] = 0;
		break;
	case 4: // ˫ռ��
		Sequence[AlphaSite] = 1;
		Sequence[BetaSite] = 1;
		break;
	default: // status �����������׳��쳣
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