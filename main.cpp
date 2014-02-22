#include <QtCore/QCoreApplication>
#include <iostream>

#include "bracket.h"
#include "FiniteElementMatrix.h"

void test_bracket(unsigned N, unsigned p_i, double g_i);
void test_vectBracket(unsigned N, unsigned p_i, double g_i);

//------------------------------------------------------------------

int main(int argc, char *argv[])
{
    //test_bracket(2,1,1);
    test_vectBracket(2,1,1);
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

void test_vectBracket(unsigned N, unsigned p_i, double g_i)
{
    FiniteElementMatrix f(1);

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
}
