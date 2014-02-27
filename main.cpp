#include <QtCore/QCoreApplication>
#include <iostream>

#include "bracket.h"
#include "FiniteElementMatrix.h"

void test_bracket(unsigned N, unsigned p_i, double g_i);
void test_vectBracket(unsigned N, unsigned p_i, double g_i);
void Test_nm();
void Test_DefVector();
void Test_AddVTVProduct();

//------------------------------------------------------------------

int main(int argc, char *argv[])
{
    //test_bracket(2,1,1);
    //test_vectBracket(2,1,1);
    //Test_nm();
    //Test_DefVector();
    Test_AddVTVProduct();
    return 0;
}

//------------------------------------------------------------------

void test_bracket(unsigned N, unsigned p_i, double g_i)
{
    Bracket b1(N);
    Bracket b2(N);

    std::vector<Power_t> powers;
    Power_t p = {p_i,p_i,p_i,p_i};

    for (unsigned i=0;i<N;i++)
    {
        powers.push_back(p);
    }

    b1.SetPowers(powers);
    b2.SetPowers(powers);

    std::vector<double> gains;
    for (unsigned i=0;i<N;i++)
    {
        gains.push_back(g_i);
    }

    b1.SetGains(gains);
    b2.SetGains(gains);

    Bracket b3 = b1 * b2;



    b1.ShowElements();

    b2.ShowElements();

    b3.ShowElements();

    Bracket b4 = b3 * b2*b1;

    b4.ShowElements();
}

/*void test_vectBracket(unsigned N, unsigned p_i, double g_i)
{
    double simplex_peaks[4][3];

    for (unsigned i=0;i<4;i++)
    {
        for (unsigned j=0;j<3;j++)
        {
            simplex_peaks[i][j]=i;
        }
    }

    FiniteElementMatrix f(1,simplex_peaks);

    Bracket b1(N);
    Bracket b2(N);

    std::vector<Power_t> powers;
    Power_t p = {p_i,p_i,p_i,p_i};

    for (unsigned i=0;i<N;i++)
    {
        powers.push_back(p);
    }

    b1.SetPowers(powers);
    powers.clear();

    Power_t p2 = {2*p_i,2*p_i,2*p_i,2*p_i};

    for (unsigned i=0;i<N;i++)
    {
        powers.push_back(p2);
    }
    b2.SetPowers(powers);

    std::vector<double> gains;
    for (unsigned i=0;i<N;i++)
    {
        gains.push_back(g_i);
    }

    b1.SetGains(gains);
    b2.SetGains(gains);

    f.AddToVectBracket(b1.GetGains(),b1.GetPowers());
    f.AddToVectBracket(b2.GetGains(),b2.GetPowers());

    f.ShowVectBracket();
}*/

/*void Test_nm()
{
    double simplex_peaks[4][3];

    for (unsigned i=0;i<4;i++)
    {
        for (unsigned j=0;j<3;j++)
        {
            simplex_peaks[i][j]=i;
        }
    }

    FiniteElementMatrix f(1,simplex_peaks);
    std::vector<unsigned> vect;
    vect=f.Def_nm(2,3);
    unsigned a=1;
}*/

/*void Test_DefVector()
{
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

    FiniteElementMatrix f(1,simplex_peaks);
    std::vector<double> vect;
    vect=f.DefVector(1);
    unsigned a=1;
}*/

void Test_AddVTVProduct()
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

    f.AddVTVProduct(a,b,c,d,M,3,4,3,4);
    unsigned var=1;
}
