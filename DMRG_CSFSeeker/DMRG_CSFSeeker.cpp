// DMRG_CSFSeeker.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include "StateFunctions.h"

using namespace std;

int main(int argc, char** argv)
{
	time_t start, finish;

	time(&start);

	StateFunctions testobj(22, 22);
	testobj.TEST(6, 6, 1, 22, 22, 1);
	testobj.evolute(0, 1.0);

	cout << "--------------------------------------------" << endl;
	time(&finish);
	double duration = difftime(finish, start);
	cout << "--> time: " << duration << " s" << endl;
	cout << "--------------------------------------------" << endl;

	return 0;
}
