#include "bracket.h"

#include <vector>
#include <cstdlib>
#include <ctime>
#include <iostream>

using namespace std;

std::vector<double> VectProduct(std::vector<double> a, std::vector<double> b)
{
    std::vector<double> result;
    result.push_back(a[1]*b[2]-a[2]*b[1]);
    result.push_back(a[2]*b[0]-a[0]*b[2]);
    result.push_back(a[0]*b[1]-a[1]*b[0]);
    return result;
}
//--------------------------------------------------------------------------------------------------------
void VectProduct(std::vector<double>* a, std::vector<double>* b, std::vector<double>* result)
{
	result->assign(3, 0.0);
	result->at(0) = a->at(1)*b->at(2)-a->at(2)*b->at(1);
	result->at(1) = a->at(2)*b->at(0)-a->at(0)*b->at(2);
    result->at(2) = a->at(0)*b->at(1)-a->at(1)*b->at(0);
}
//---------------¬екторное произведение двух векторов, элементами кот. €вл€ютс€ скобки--------------------
void VectBracketProduct(vector<Bracket>& a, vector<Bracket>& b, vector<Bracket>& result)
{

    Bracket b1=a[1]*b[2];
    Bracket b2=a[2]*b[1];

    result.push_back(b1-b2);

    b1.BracketCleanUp();
    b2.BracketCleanUp();

    b1=a[2]*b[0];
    b2=a[0]*b[2];

    result.push_back(b1-b2);

    b1.BracketCleanUp();
    b2.BracketCleanUp();

    b1=a[0]*b[1];
    b2=a[1]*b[0];

    result.push_back(b1-b2);

    //b1.BracketCleanUp();
    //b2.BracketCleanUp();
}
//----------------------------------------------------------------------------------------------------
//-------------—кал€рное произведение двух векторов---------------------------------------------------
double ScalarProduct(std::vector<double> a, std::vector<double> b)
{
    double result;
    result=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    return result;
}

//--------------ѕроизведение вектора, тензора и вектора в общем виде----------------------------------
void GeneralVectorTensorVectorProduct(std::vector<Bracket>* vect1,
									  //std::vector<Bracket>& vect2, double** M, Bracket& result)
									  std::vector<Bracket>* vect2, double** M, Bracket* result,
									  Bracket* br_sum, Bracket* br_prod)
{
    double m_Dim=3;
	unsigned unit = 1;

	result->SetSizePtr(&unit);

    for (unsigned i=0;i<m_Dim;i++)
    {
		br_sum->SetSizePtr(&unit);

        for (unsigned k=0;k<m_Dim;k++)
        {
			Mult(&(vect1->at(k)), &M[k][i], br_prod);
			br_sum->Plus(br_prod);
        }
        Mult(br_sum,&(vect2->at(i)),br_prod);

		result->Plus(br_prod);
    }
}

//„исловой вектор vect умножаетс€ на число number
void MultNumber(std::vector<double>& vect, double number)
{
	for (unsigned i=0; i<vect.size();i++)
	{
		vect[i]=vect[i]*number;
	}
}

//---------------------------------  числовому вектору sum прибавл€етс€ вектор added----------------------------------------
void SumVector(std::vector<double>& sum,const std::vector<double> added)
{
	if ((sum.size())==(added.size()))
		for (unsigned i=0; i<sum.size();i++)
		{
			sum[i]=sum[i]+added[i];
		}
	else
		std::cout << "–азмерности векторов не совпадают!" << std::endl;
}

//-----------------------------ѕроизведение числового вектора, матрицы и числового вектора------------------------------------
double NumericalVectorTensorVectorProduct(std::vector<double>& vect1,
                                          std::vector<double>& vect2, double** M)
{
	unsigned m_Dim=3;
	double s;
	double result=0.0;

	for (unsigned k=0;k<m_Dim;k++)
	{
		s=0.0;
		for (unsigned i=0;i<m_Dim;i++)
		{
			s+=M[i][k]*vect1[i];
		}
		result+=s*vect2[k];
	}
	return result;
}

//-----------------------------ѕроизведение числового вектора, матрицы и числового вектора------------------------------------
double NumericalVectorTensorVectorProduct(double* arr1,
                                          double* arr2, double** M)
{
	unsigned m_Dim=3;
	double s;
	double result=0.0;

	for (unsigned k=0;k<m_Dim;k++)
	{
		s=0.0;
		for (unsigned i=0;i<m_Dim;i++)
		{
			s+=M[i][k]*arr1[i];
		}
		result+=s*arr2[k];
	}
	return result;
}



