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
//---------------��������� ������������ ���� ��������, ���������� ���. �������� ������--------------------
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
//-------------��������� ������������ ���� ��������---------------------------------------------------
double ScalarProduct(std::vector<double> a, std::vector<double> b)
{
    double result;
    result=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    return result;
}

//--------------������������ �������, ������� � ������� � ����� ����----------------------------------
Bracket GeneralVectorTensorVectorProduct(std::vector<Bracket>& vect1,
                                                              std::vector<Bracket>& vect2, double** M)
{
    double m_Dim=3;

    //������ ������, ������� ������ �� ��������
	std::vector<GainPower_t> zero_terms;
	GainPower_t trm={0.0, 0, 0, 0, 0};
	zero_terms.push_back(trm);

    Bracket result(1);

    Bracket br_product(1);

    Bracket br_sum(1);

    for (unsigned i=0;i<m_Dim;i++)
    {
        br_sum.SetTerms(zero_terms);

        for (unsigned k=0;k<m_Dim;k++)
        {
            br_product=vect1[k]*(M[k][i]);
            br_sum=br_sum+br_product;
        }
        br_sum=br_sum*vect2[i];

        result=result+br_sum;
    }
    return result;
}

//�������� ������ vect ���������� �� ����� number
void MultNumber(std::vector<double>& vect, double number)
{
	for (unsigned i=0; i<vect.size();i++)
	{
		vect[i]=vect[i]*number;
	}
}

//---------------------------------� ��������� ������� sum ������������ ������ added----------------------------------------
void SumVector(std::vector<double>& sum,const std::vector<double> added)
{
	if ((sum.size())==(added.size()))
		for (unsigned i=0; i<sum.size();i++)
		{
			sum[i]=sum[i]+added[i];
		}
	else
		std::cout << "����������� �������� �� ���������!" << std::endl;
}

//-----------------------------������������ ��������� �������, ������� � ��������� �������------------------------------------
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
