#include "bracket.h"
#include "FiniteElementMatrix.h"
#include "VectFunctions.h"
#include "BracketFunctions.h"
#include "InterplNode.h"


#ifdef _QTGUI_
#include <QtCore/QCoreApplication>
#endif
#include <iostream>
#include <algorithm>
#include <iterator>

#include <tchar.h>
#include <windows.h>


using namespace std;

void test_Matrixs();

void display(unsigned rows, unsigned columns, double** arr);

void StartTimeMeasure(LARGE_INTEGER& StartPerformCount);
double StopTimeMeasure(LARGE_INTEGER& StartPerformCount);

void test_DefFaceInd();
void test_FormFaceFunc();

//------------------------------------------------------------------

int main(int argc, char *argv[])
{
	test_Matrixs();
	//test_DefFaceInd();
	//test_FormFaceFunc();
	system("pause");
    return 0;
}

//------------------------------------------------------------------


void test_Matrixs()
{
	//Создаю матрицу, содержащую вершины элемента
	double** simplex_peaks;
	simplex_peaks = new double*[4];

	for (unsigned i=0;i<4;i++)
	{
		simplex_peaks[i] = new double[3];
	}

    simplex_peaks[0][0]=2;
    simplex_peaks[0][1]=0;
    simplex_peaks[0][2]=0;

    simplex_peaks[1][0]=0;
    simplex_peaks[1][1]=1;
    simplex_peaks[1][2]=0;

    simplex_peaks[2][0]=0;
    simplex_peaks[2][1]=0;
    simplex_peaks[2][2]=1;

    simplex_peaks[3][0]=0;
    simplex_peaks[3][1]=0;
    simplex_peaks[3][2]=0;
    //создаю матрицы диэлектрической и магнитной проницаемости
    double** Mu;
	double** Eps;
	Mu = new double* [3];
	Eps = new double* [3];

    for (unsigned i=0;i<3;i++)
    {
		Mu[i] = new double [3];
		Eps[i] = new double [3];
        for (unsigned j=0;j<3;j++)
        {
            if (i==j)
            {
                Mu[i][j]=1;
                Eps[i][j]=1;
            }
            else
            {
                Eps[i][j]=0;
                Mu[i][j]=0;
            }
        }
    }
	Eps[0][1]=1;
	Eps[0][2]=2;

	cout << "----------Peaks-----------" << endl;
	display(4, 3, simplex_peaks);
	cout << "----------Mu-----------" << endl;
	display(3, 3, Mu);
	cout << "----------Eps-----------" << endl;
	display(3, 3, Eps);

	vector<double> a,b;
	a.push_back(1.0);
	a.push_back(2.0);
	a.push_back(3.0);

	b.push_back(-1.0);
	b.push_back(-2.0);
	b.push_back(-3.0);

	//double pr=NumericalVectorTensorVectorProduct(a,b,Eps);
    FiniteElementMatrix f(1,simplex_peaks,Eps,Mu);
}



void display(unsigned rows, unsigned columns, double** arr)
{
	for (unsigned i=0;i<rows;i++)
	{
		for (unsigned j=0;j<columns; j++)
			cout << arr[i][j] << "  ";
		cout << endl;
	}
}

void StartTimeMeasure(LARGE_INTEGER& StartPerformCount)
{
	QueryPerformanceCounter (&StartPerformCount);
}

double StopTimeMeasure(LARGE_INTEGER& StartPerformCount)
{
	LARGE_INTEGER Frequency, StopPerformCount;
	QueryPerformanceCounter (&StopPerformCount);
	QueryPerformanceFrequency(&Frequency);
	return ((double)(StopPerformCount.QuadPart - StartPerformCount.QuadPart) / (double)Frequency.QuadPart * 1.E3);
}

void test_DefFaceInd()
{
	unsigned* arr;
	arr = new unsigned[3];

	unsigned* non_zero_arr;
	non_zero_arr = new unsigned[3];

	DefFaceInd(0,12,11,15,arr,non_zero_arr);

	for (unsigned i=0;i<3;i++)
	cout << arr[i] << " ";

	for (unsigned i=0;i<3;i++)
	cout << non_zero_arr[i] << " ";
}

void test_FormFaceFunc()
{
	//Создаю матрицу, содержащую вершины элемента
	double** simplex_peaks;
	simplex_peaks = new double*[4];

	for (unsigned i=0;i<4;i++)
	{
		simplex_peaks[i] = new double[3];
	}

    simplex_peaks[0][0]=2;
    simplex_peaks[0][1]=0;
    simplex_peaks[0][2]=0;

    simplex_peaks[1][0]=0;
    simplex_peaks[1][1]=1;
    simplex_peaks[1][2]=0;

    simplex_peaks[2][0]=0;
    simplex_peaks[2][1]=0;
    simplex_peaks[2][2]=1;

    simplex_peaks[3][0]=0;
    simplex_peaks[3][1]=0;
    simplex_peaks[3][2]=0;
    //создаю матрицы диэлектрической и магнитной проницаемости
    double** Mu;
	double** Eps;
	Mu = new double* [3];
	Eps = new double* [3];

    for (unsigned i=0;i<3;i++)
    {
		Mu[i] = new double [3];
		Eps[i] = new double [3];
        for (unsigned j=0;j<3;j++)
        {
            if (i==j)
            {
                Mu[i][j]=1;
                Eps[i][j]=1;
            }
            else
            {
                Eps[i][j]=0;
                Mu[i][j]=0;
            }
        }
    }
	Eps[0][1]=1;
	Eps[0][2]=2;

	cout << "----------Peaks-----------" << endl;
	display(4, 3, simplex_peaks);
	cout << "----------Mu-----------" << endl;
	display(3, 3, Mu);
	cout << "----------Eps-----------" << endl;
	display(3, 3, Eps);

    FiniteElementMatrix f(0,simplex_peaks,Eps,Mu);
	
	GainPower_t trm1={1.0, 1, 1, 0, 0};
	GainPower_t trm2={3.0, 0, 0, 0, 2};

	vector<GainPower_t> terms;
	terms.push_back(trm1);
	terms.push_back(trm2);

	Bracket br(terms);

	unsigned* index_array;
	unsigned* non_zero_array;
	index_array = new unsigned[3];
	non_zero_array = new unsigned[3];

	index_array[0] = 1;
	index_array[1] = 2;
	index_array[2] = 3;

	non_zero_array[0] = 1;
	non_zero_array[1] = 2;
	non_zero_array[2] = 4;

	GainPower_t trm3={1.0, 0, 0, 0, 0};
	vector<Bracket> vect_br;
	vector<GainPower_t> vect_t;
	vect_t.push_back(trm3);
	Bracket unit_br(vect_t);
	vect_br.push_back(unit_br);
	vect_br.push_back(unit_br);
	vect_br.push_back(unit_br);

	f.FormFaceFunc(&br, index_array, non_zero_array, &vect_br, 2);
}



