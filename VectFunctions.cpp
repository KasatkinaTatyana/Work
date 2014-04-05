#include "bracket.h"

#include <vector>
#include <cstdlib>
#include <ctime>

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
//---------------Векторное произведение двух векторов, элементами кот. являются скобки--------------------
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
//-------------Скалярное произведение двух векторов---------------------------------------------------
double ScalarProduct(std::vector<double> a, std::vector<double> b)
{
    double result;
    result=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    return result;
}

//--------------Произведение вектора, тензора и вектора в общем виде----------------------------------
Bracket GeneralVectorTensorVectorProduct(std::vector<Bracket>& vect1,
                                                              std::vector<Bracket>& vect2, double M[][3])
{
    double m_Dim=3;

    //создаю скобку, которая ничего не содержит
    std::vector<double> zero_gains;
    zero_gains.push_back(0.0);
    std::vector<Power_t> zero_powers;
    Power_t p={0, 0, 0, 0};
    zero_powers.push_back(p);

    Bracket result(1);
    result.SetGains(zero_gains);
    result.SetPowers(zero_powers);

    Bracket br_product=vect1[0]*(M[0][0]);

    Bracket br_sum(1);

    for (unsigned i=0;i<m_Dim;i++)
    {
        br_sum.SetGains(zero_gains);
        br_sum.SetPowers(zero_powers);

        for (unsigned k=0;k<m_Dim;k++)
        {
            if ((i!=0)||(k!=0))
            {
                br_product=vect1[k]*(M[k][i]);
            }
            br_sum=br_sum+br_product;
        }
        br_sum=br_sum*vect2[i];

        result=result+br_sum;

        br_sum.BracketCleanUp();
    }
    return result;
}
