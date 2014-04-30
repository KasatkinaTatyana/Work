#include "bracket.h"
#include "FiniteElementMatrix.h"
#include "VectFunctions.h"
#include "BracketFunctions.h"


#ifdef _QTGUI_
#include <QtCore/QCoreApplication>
#endif
#include <iostream>
#include <algorithm>
#include <iterator>


using namespace std;

void test_bracket(unsigned N, unsigned p_i, double g_i);
void test_simplify(unsigned N, unsigned p_i, double g_i);
void test_prod();
//void test_vectBracket(unsigned N, unsigned p_i, double g_i);
void test_Matrixes();
//void test_integrate();
void display(unsigned rows, unsigned columns, double** arr);

//------------------------------------------------------------------

int main(int argc, char *argv[])
{
    //test_bracket(2,1,1);
	//test_simplify(2,1,1);
    //test_vectBracket(2,1,1);
	//test_prod();
	test_Matrixes();
	//test_integrate();
	system("pause");
    return 0;
}

//------------------------------------------------------------------

/*void test_bracket(unsigned N, unsigned p_i, double g_i)
{
	Bracket br(3);

	std::vector<double> gains;
	std::vector<Power_t> powers;

	gains.push_back(1.0);
	gains.push_back(2.0);
	gains.push_back(3.0);
	gains.push_back(4.0);

	Power_t pw={1, 0, 0, 0};
	powers.push_back(pw);
	
	pw.p1=0;
	pw.p2=1;
	powers.push_back(pw);

	pw.p2=0;
	pw.p3=1;
	powers.push_back(pw);

	pw.p3=0;
	pw.p4=1;
	powers.push_back(pw);

	br.SetGains(gains);
	br.SetPowers(powers);

	std::vector<double> gains1;
	std::vector<Power_t> powers1;

	gains1=br.GetGains();
	powers1=br.GetPowers();

	Bracket b1(2);

	b1 = br;

	Bracket b2(1);
	pw.p1=1;
	pw.p2=1;
	pw.p3=0;
	pw.p4=0;

	gains.clear();
	powers.clear();

	gains.push_back(1.0);
	powers.push_back(pw);
	b2.SetGains(gains);
	b2.SetPowers(powers);

	Bracket b3(2);
	gains.clear();
	powers.clear();

	gains.push_back(1.0);
	gains.push_back(1.0);

	pw.p2=0;
	powers.push_back(pw);


	pw.p3=1;
	powers.push_back(pw);

	b3.SetGains(gains);
	b3.SetPowers(powers);

	//Bracket b4=b2*b3;
	Bracket b4=b2+b3;

	double aa=1;
}

void test_simplify(unsigned N, unsigned p_i, double g_i)
{
	Bracket br(4);
	std::vector<double> gains;
	std::vector<Power_t> powers;

	Power_t pw={1, 1, 0, 0};
	gains.push_back(1.0);
	gains.push_back(2.0);

	gains.push_back(5.0);
	gains.push_back(0.0);
	gains.push_back(3.0);

	powers.push_back(pw);
	powers.push_back(pw);

	Power_t pw1={1, 1, 1, 1};
	powers.push_back(pw1);

	Power_t pw2={0, 1, 0, 0};
	powers.push_back(pw2);

	powers.push_back(pw);

	br.SetGains(gains);
	br.SetPowers(powers);

	br.ShowElements();

	Bracket br1(2);

	gains.clear();
	powers.clear();

	gains.push_back(2.0);
	gains.push_back(3.0);

	Power_t pw3={0, 2, 0, 0};

	powers.push_back(pw3);
	powers.push_back(pw1);

	br1.SetGains(gains);
	br1.SetPowers(powers);

	Bracket br2=br - br1;
	double a=0;
}

void test_prod()
{
	Bracket br1(2);

	std::vector<double> gains;
	std::vector<Power_t> powers;

	gains.clear();
	powers.clear();

	gains.push_back(2.0);
	gains.push_back(3.0);

	Power_t pw3={0, 2, 0, 0};
	Power_t pw1={1, 1, 1, 1};
	Power_t pw2={1, 0, 0, 0};

	powers.push_back(pw3);
	powers.push_back(pw1);

	br1.SetGains(gains);
	br1.SetPowers(powers);

	Bracket br2(1);

	gains.clear();
	powers.clear();

	gains.push_back(1.0);
	powers.push_back(pw2);

	br2.SetGains(gains);
	br2.SetPowers(powers);

	br1*=br2;

	br1.ShowElements();
}*/

/*void Test_AddVTVProduct()
{
    std::vector<double> a, b, c, d;

    a.push_back(1);
    a.push_back(0);
    a.push_back(0);

    b.push_back(0);
    b.push_back(1);
    b.push_back(0);

    c.push_back(0);
    c.push_back(0);
    c.push_back(1);

    d.push_back(1);
    d.push_back(1);
    d.push_back(1);

    //Создаю матрицу, содержащую вершины элемента
    double simplex_peaks[4][3];

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
    double Mu[3][3], Eps[3][3], M[3][3];

    for (unsigned i=0;i<3;i++)
    {
        for (unsigned j=0;j<3;j++)
        {
            if (i==j)
            {
                Mu[i][j]=1;
                Eps[i][j]=1;
                M[i][j]=1;
            }
            else
            {
                Eps[i][j]=0;
                Mu[i][j]=0;
                M[i][j]=0;
            }
        }
    }
    FiniteElementMatrix f(1,simplex_peaks,Eps,Mu);
}*/

void test_Matrixes()
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

/*void test_integrate()
{
	GainPower_t trm={1.0, 1, 0, 1, 0};
	std::vector<GainPower_t> terms(1,trm);
	Bracket br(terms);
	//Создаю матрицу, содержащую вершины элемента
    double simplex_peaks[4][3];

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
    double Mu[3][3], Eps[3][3];

    for (unsigned i=0;i<3;i++)
    {
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

	FiniteElementMatrix f(0,simplex_peaks,Eps,Mu);

	double Int = f.Integrate(br);
	std::vector<Bracket> vect_br(1,br);
	std::vector<double> vect_numb;
	f.NumIntegration(vect_br,vect_numb);
	std::cout << Int << std::endl;
	std::cout << vect_numb[0] << std::endl;
	system("pause");
}*/

void display(unsigned rows, unsigned columns, double** arr)
{
	for (unsigned i=0;i<rows;i++)
	{
		for (unsigned j=0;j<columns; j++)
			cout << arr[i][j] << "  ";
		cout << endl;
	}
}