#ifndef CONFIGURATION_CLASS_H
#define CONFIGURATION_CLASS_H

#include <iostream>
#include <string>
#include "RandomNumber.h"
#include "Array.h"

// ����:
// ��̬�� 2Norbs ����������� ababab... ���������м�¼��
// ���е�ÿһλ�� 0 ��ʾ�޵��ӣ��� 1 ��ʾ�е���
class Configuration
{
public:
	// constructors
	// ���캯�������ܵĲ���Ϊ���Կռ�Ĺ�������������Լ���̬���������ض�
	// �������ض� S=2s+1 ��Ĭ��ֵΪ 1��������̬��
	// ��Ϊ�����˸��ݸ���������������̬�Ĺ��ܣ��ʳ�ʼ�趨�ĵ�����Ҳ������Ҫ���ӹ��캯���п���ʡȥ��
	Configuration(int orbs, int elec = 0, int spin = 1);
	Configuration(const Configuration& rhs); // ���ƹ��캯��
	// ����������չ�����Ĺ��캯��
	Configuration(const Array& CoreAlpha, const Array& CoreBeta,
		int orbs, int elec, int spin,
		int CASCoreOrbs, int CASCoreElecAlpha, int CASCoreElecBeta);
	// ���Թ���
	Configuration(const Array& CoreAlpha, const Array& CoreBeta,
		int orbs, int elec, int spin,
		int CASCoreOrbs, int CASCoreElecAlpha, int CASCoreElecBeta,
		int LargerSpaceOrbs, int NLargerSpaceElecAlpha, int NLargerSpaceElecBeta);
	Configuration(const std::string& rhs); // �� 444111 �������й���һ����̬
	// the deconstructor
	~Configuration();
	// member functions
	int& alpha(int index); // ��������Ϊ index �� alpha ���������״̬
	const int& alpha(int index) const; // ��ֵ��������Ϊ index �� alpha ���������״̬
	int& beta(int index); // ��������Ϊ index �� beta ���������״̬
	const int& beta(int index) const; // ��ֵ��������Ϊ index �� beta ���������״̬
	int& operator[](int index); // ��������Ϊ index �����������״̬
	const int& operator[](int index) const; // ��ֵ��������Ϊ index �����������״̬
	double coeff() const; // ������̬ϵ��
	void initialize(); // �� Sequence ��ÿһλ���� 0
	void SetHartreeFockDet(); // ����̬�趨Ϊ Hartree-Fock ����ʽ��ʾ�����ռ����ʽ
	void SetBySeed(const std::string& rhs); // �� 444111 �������й���һ����̬
	const Configuration& randomize(); // ������ɵ����Ų����У����չ���������������������ضȣ���ѭ Paili ������ԭ��
	Configuration& mutate(double rate, double StepSize = 0.1); // ���ո����ĸ��ʷ�������
	Configuration& MutateStrategy2(); // �������2
	Configuration& operator=(const Configuration& rhs); // ���ظ�ֵ������
	// �������������
	friend std::ostream& operator<<(std::ostream& os, const Configuration& rhs);
	// ��Ԫ�ࣺ̬����
	friend class StateFunctions;
	// ���ڰ�����̬ϵ������ֵ�Ӵ�С��������ĺ�������
	friend struct SortByCoefficient;
	// ���ڰ������н�������ĺ�������
	friend struct SortBySequence;
	// ����ִ�� unique ������������Ȳ�����
	friend bool operator==(const Configuration& lhs, const Configuration& rhs);
private:
	int SpinOrbs; // ��������ĸ���
	int Nelec; // �ܵĵ�����
	int Spin; // �������ضȣ�S=2s+1
	double Coefficient; // ��̬ϵ������Ϊ����Ӧ�ȵı���
	// ababab... ���������У�ÿһλ�� 0 ��ʾ�޵��ӣ��� 1 ��ʾ�е���
	int* Sequence;
	// ˽�и�������
	// ������֪������ Sequence �����ܵĵ����� Nelec ���������ض� Spin
	void ComputeNelecSpin();
	// ��ָ���ռ���λ�� site ��ʼ������ status��1,2,3,4��������� Sequence �е�����λ��
	void SetOccupationStatus(int site, int status);
};

// ���ڽ���̬����ϵ������ֵ�Ӵ�С��������ĺ�������
struct SortByCoefficient
{
	bool operator()(const Configuration& lhs, const Configuration& rhs);
};

// ���ڽ���̬�������н���Ψһ����ĺ�������˳�򲢲���Ҫ��ֻ��Ҫ��֤��ͬ���������������
struct SortBySequence
{
	bool operator()(const Configuration& lhs, const Configuration& rhs);
};

#endif // !CONFIGURATION_CLASS_H