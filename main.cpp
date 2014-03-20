#include <QtCore/QCoreApplication>
#include <iostream>

#include "bracket.h"
#include "FiniteElementMatrix.h"

using namespace std;

void test_bracket(unsigned N, unsigned p_i, double g_i);
void test_vectBracket(unsigned N, unsigned p_i, double g_i);
void Test_nm();
void Test_DefVector();
void Test_AddVTVProduct();
void Test_AddSilvester();
void Test_MetrMatrix();
void test_GeneralVectorTensorVectorProduct(unsigned N, unsigned p_i, double g_i);
void test_ProductGradKsi();
void test_RotorCalc(unsigned N, unsigned p_i, double g_i);

/*double Fact(unsigned N);
double m_ArrayFact[26];

inline double CalcFact(unsigned N)
{
    return m_ArrayFact[N];
}

double Fact(unsigned N)
{
    double f = 1.0;

    for (unsigned i=0;i<N;i++)
    {
        f *= (i+1);
    }

    return f;
}

void PrintFact(unsigned TO)
{
    for (unsigned i=0;i<TO+1;i++)
    {
        cout<<i<<": "<<Fact(i)<<endl;
    }
}*/


//------------------------------------------------------------------

int main(int argc, char *argv[])
{
    //test_bracket(2,1,1);
    //test_vectBracket(2,1,1);
    //Test_nm();
    //Test_DefVector();
    //Test_AddVTVProduct();
    //Test_AddSilvester();
    //PrintFact(20);
    //Test_MetrMatrix();
    //test_GeneralVectorTensorVectorProduct(2,1,1);
    //test_ProductGradKsi();
    test_RotorCalc(2,1,1);
    return 0;
}

//------------------------------------------------------------------

void test_bracket(unsigned N, unsigned p_i, double g_i)
{
    Bracket b1(N);
    Bracket b2(N);

    std::vector<Power_t> powers, powers1;
    Power_t p = {p_i, p_i, p_i, p_i};
    Power_t p1 = {2*p_i, 2*p_i, 2*p_i, 2*p_i};

    for (unsigned i=0;i<N;i++)
    {
        powers.push_back(p);
        powers1.push_back(p1);
    }

    b1.SetPowers(powers);
    b2.SetPowers(powers1);

    std::vector<double> gains, gains1;
    for (unsigned i=0;i<N;i++)
    {
        gains.push_back(g_i);
        gains1.push_back(g_i*2.0);
    }

    b1.SetGains(gains);
    b2.SetGains(gains1);


    Bracket zero_br(1);
    std::vector<double> zero_gains;
    zero_gains.push_back(0.0);
    std::vector<Power_t> zero_powers;
    Power_t zero_p={0, 0, 0, 0};
    zero_powers.push_back(zero_p);

    zero_br.SetGains(zero_gains);
    zero_br.SetPowers(zero_powers);

    Bracket bbb=b1+zero_br;

    Bracket br_minus=b2-zero_br;

    Bracket b3 = b2*(10.0*2.0);
    b3=b3+b1;




    b1.ShowElements();

    b2.ShowElements();

    b3.ShowElements();

    Bracket b4 = b3 * b2*b1;

    b4.ShowElements();

    Bracket b5=b1*10.0;

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

    //������ �������, ���������� ������� ��������
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
    //������ ������� ��������������� � ��������� �������������
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

    //f.AddVTVProduct(a,b,c,d,M,3,4,3,4);
}

void Test_AddSilvester()
{
    //������ �������, ���������� ������� ��������
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
    //������ ������� ��������������� � ��������� �������������
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
    //f.AddSilvester(1,2,2,1,f.);
}

void Test_MetrMatrix()
{
    //������ �������, ���������� ������� ��������
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
    //������ ������� ��������������� � ��������� �������������
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

    f.MatrixInit();
}

void test_GeneralVectorTensorVectorProduct(unsigned N, unsigned p_i, double g_i)
{
    /*Bracket b1(N);
    Bracket b2(N);

    std::vector<Power_t> powers, powers1;
    Power_t p = {p_i, p_i, p_i, p_i};
    Power_t p1 = {2*p_i, 2*p_i, 2*p_i, 2*p_i};

    for (unsigned i=0;i<N;i++)
    {
        powers.push_back(p);
        powers1.push_back(p1);
    }

    b1.SetPowers(powers);
    b2.SetPowers(powers1);

    std::vector<double> gains, gains1;
    for (unsigned i=0;i<N;i++)
    {
        gains.push_back(g_i);
        gains1.push_back(g_i*2.0);
    }

    b1.SetGains(gains);
    b2.SetGains(gains1);

    Bracket b3 = b2*10.0;

    Bracket vect1[3]={b1, b2, b3};
    Bracket vect2[3]={b3, b1, b2};



    //������ �������, ���������� ������� ��������
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
    //������ ������� ��������������� � ��������� �������������
    double Mu[3][3], Eps[3][3], M[3][3];

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

    M[0][0]=0;
    M[0][1]=1;
    M[0][2]=2;

    M[1][0]=3;
    M[1][1]=4;
    M[1][2]=5;

    M[2][0]=6;
    M[2][1]=7;
    M[2][2]=8;

    FiniteElementMatrix f(1,simplex_peaks,Eps,Mu);
    Bracket b4=f.GeneralVectorTensorVectorProduct(vect1,vect2,M);
    b4.ShowElements();

    //b4.ShowElements();*/
}

void test_ProductGradKsi()
{
    //������ �������, ���������� ������� ��������
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
    //������ ������� ��������������� � ��������� �������������
    double Mu[3][3], Eps[3][3], M[3][3];

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

    M[0][0]=0;
    M[0][1]=1;
    M[0][2]=2;

    M[1][0]=3;
    M[1][1]=4;
    M[1][2]=5;

    M[2][0]=6;
    M[2][1]=7;
    M[2][2]=8;

    FiniteElementMatrix f(1,simplex_peaks,Eps,Mu);
    std::vector<double> vect=f.ProductGradKsi(1,3);

    for (unsigned i=0;i<vect.size();i++)
    {
        std::cout<< " " << vect[i];
    }

}

void test_RotorCalc(unsigned N, unsigned p_i, double g_i)
{
    //������ �������, ���������� ������� ��������
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
    //������ ������� ��������������� � ��������� �������������
    double Mu[3][3], Eps[3][3], M[3][3];

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

    M[0][0]=0;
    M[0][1]=1;
    M[0][2]=2;

    M[1][0]=3;
    M[1][1]=4;
    M[1][2]=5;

    M[2][0]=6;
    M[2][1]=7;
    M[2][2]=8;

    FiniteElementMatrix f(1,simplex_peaks,Eps,Mu);

    Bracket b1(N);
    Bracket b2(N);

    std::vector<Power_t> powers, powers1;
    Power_t p = {p_i, p_i, p_i, p_i};
    Power_t p1 = {2*p_i, 2*p_i, 2*p_i, 2*p_i};

    for (unsigned i=0;i<N;i++)
    {
        powers.push_back(p);
        powers1.push_back(p1);
    }

    b1.SetPowers(powers);
    b2.SetPowers(powers1);

    std::vector<double> gains, gains1;
    for (unsigned i=0;i<N;i++)
    {
        gains.push_back(g_i);
        gains1.push_back(g_i*2.0);
    }

    b1.SetGains(gains);
    b2.SetGains(gains1);

    //--------------------------------------------test_VectBracketProduct--------------------------------
    powers.clear();
    gains.clear();

    std::vector<Bracket> a_left;
    std::vector<Bracket> b_right;

    Bracket brack(1);

    Power_t pow1={1,0,0,0};
    Power_t pow0={0,0,0,0};
    Power_t pow2={0,1,0,0};

    powers.push_back(pow1);
    gains.push_back(1);

    brack.SetGains(gains);
    brack.SetPowers(powers);

    a_left.push_back(brack);

    gains.clear();
    powers.clear();
    brack.BracketCleanUp();

    powers.push_back(pow0);
    gains.push_back(0);

    brack.SetGains(gains);
    brack.SetPowers(powers);

    a_left.push_back(brack);

    gains.clear();
    powers.clear();
    brack.BracketCleanUp();

    powers.push_back(pow0);
    gains.push_back(0);

    brack.SetGains(gains);
    brack.SetPowers(powers);

    a_left.push_back(brack);
    //---------------------------------------------------------

    powers.clear();
    gains.clear();
    brack.BracketCleanUp();

    powers.push_back(pow0);
    gains.push_back(0);

    brack.SetGains(gains);
    brack.SetPowers(powers);

    b_right.push_back(brack);

    powers.clear();
    gains.clear();
    brack.BracketCleanUp();

    powers.push_back(pow2);
    gains.push_back(1);

    brack.SetGains(gains);
    brack.SetPowers(powers);

    b_right.push_back(brack);

    powers.clear();
    gains.clear();
    brack.BracketCleanUp();

    powers.push_back(pow0);
    gains.push_back(0);

    brack.SetGains(gains);
    brack.SetPowers(powers);

    b_right.push_back(brack);
    //-------------------------------------------------------------------------
    std::vector<Bracket> Brack_vect=f.VectBracketProduct(a_left,b_right);

    f.RotorCalc(b1,1,3);

    /*for (unsigned i=0;i<Brack_vect.size();i++)
    {
        ((Bracket)(Brack_vect[i])).ShowElements();
    }*/
}
