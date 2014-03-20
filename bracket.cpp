#include "bracket.h"

#include <iostream>

Bracket::Bracket(unsigned N) : m_N(N)
{
    BracketInit(N);
}

Bracket::Bracket(std::vector<double>& gains, std::vector<Power_t>& powers)
{
    BracketInit(gains, powers);
}

Bracket::~Bracket()
{
    BracketCleanUp();
}

void Bracket::BracketInit(unsigned N)
{
    if (N!=0)
    {
        m_N=N;
    }
    for (unsigned i=0;i<m_N;i++)
    {
        m_Gains.push_back(0);

        Power_t powers = {0,0,0,0};

        m_Powers.push_back(powers);
    }
}

void Bracket::BracketInit(std::vector<double>& gains, std::vector<Power_t>& powers)
{
    unsigned N=gains.size();
    m_N=N;
    for (unsigned i=0;i<m_N;i++)
    {
        m_Gains.push_back(gains[i]);

        m_Powers.push_back(powers[i]);
    }

}

void Bracket::BracketCleanUp()
{
    m_Gains.clear();
    m_Powers.clear();
}

void Bracket::SetGains(std::vector<double>& gains)
{
    m_Gains.clear();
    for (unsigned i=0;i<gains.size();i++)
    {
        m_Gains.push_back(gains[i]);
    }
    m_N=gains.size();
}

void Bracket::SetPowers(std::vector<Power_t>& powers)
{
    m_Powers.clear();
    for (unsigned i=0;i<powers.size();i++)
    {
        m_Powers.push_back(powers[i]);
    }
    m_N=powers.size();
}

void Bracket::ShowElements()
{
    std::cout << "======================" << std::endl;
    std::cout << "m_N = " << m_N << std::endl;

    std::cout << "m_Gains: "<<std::endl;
    for(unsigned i=0;i<m_N;i++)
    {
        if((i%8) == 0) {
            std::cout<< std::endl;
        }
        std::cout<< " " << m_Gains[i];
    }
    std::cout<< std::endl;

    std::cout << "m_Powers: "<<std::endl;
    for(unsigned i=0;i<m_N;i++)
    {
        if((i%2) == 0) {
            std::cout<< std::endl;
        }
        std::cout<<" {"<<m_Powers[i].p1<<", ";
        std::cout<<m_Powers[i].p2<<", ";
        std::cout<<m_Powers[i].p3<<", ";
        std::cout<<m_Powers[i].p4<<"} ";
    }
    std::cout<< std::endl;
    std::cout << "======================" << std::endl;
}

Bracket Bracket::operator*(Bracket& b)
{
    unsigned M = b.BracketSize()*BracketSize();
    unsigned L = b.BracketSize();
    unsigned N = BracketSize();

    Bracket result(M);

    for(unsigned i = 0; i < N; i++)
    {

        unsigned offset = i*L;

        float Ai = m_Gains.at(i);
        Power_t q1 = m_Powers.at(i);

        for(unsigned j = 0; j < L; j++)
         {

            float Bj = b.m_Gains.at(j);

            float Cij = Ai * Bj;

            result.m_Gains[offset + j] = Cij;

            Power_t q2 = b.m_Powers.at(j);

            Power_t q3 = {0,0,0,0};

            q3.p1 = q1.p1+q2.p1;
            q3.p2 = q1.p2+q2.p2;
            q3.p3 = q1.p3+q2.p3;
            q3.p4 = q1.p4+q2.p4;

            result.m_Powers[offset+j]=q3;

        }
    }
    return result;
}

Bracket Bracket::operator*(double number)
{
    Bracket result(m_N);

    std::vector<double> gains;
    for (unsigned i=0;i<m_N;i++)
    {
        gains.push_back(m_Gains[i]*number);
    }
    result.SetGains(gains);
    result.SetPowers(m_Powers);
    return result;
}

Bracket Bracket::operator*=(Bracket& b)
{
    unsigned M = b.BracketSize()*BracketSize();

    Bracket result(M);

    result=(*this)*b;
    return result;
}

Bracket Bracket::operator+(Bracket& b)
{
    //если конструкция вида 0 + скобка, то возвращается значение скобка
    unsigned flag=0;
    for (unsigned i=0;i<m_Gains.size();i++)
    {
        if (m_Gains[i]!=0.0)
            flag=1;
    }

    if (flag==0)
        return b;
    else
    {
        unsigned M = BracketSize();

        std::vector<double> gains=m_Gains;
        std::vector<Power_t> powers=m_Powers;

        flag=0;
        for (unsigned i=0;i<b.BracketSize();i++)
        {
            if ((double)((b.GetGains())[i])!=(0.0))
            {
                M++;
                flag=1;
                gains.push_back((b.GetGains())[i]);
                powers.push_back((b.GetPowers())[i]);
            }
        }

        if (flag==0)
        {
            return (*this);
        }
        else
        {
            Bracket result(M);
            result.SetGains(gains);
            result.SetPowers(powers);
            return result;
        }
    }
}

Bracket Bracket::operator-(Bracket& b)
{
    Bracket b1=b*(-1.0);
    return (*this+b1);
}

/*Bracket Bracket::operator+=(Bracket& b)
{
    unsigned flag=0;
    for (unsigned i=0;i<m_N;i++)
    {
        if ((double)(m_Gains[i])!=(0.0))
        {
            flag=1;
        }
    }

    if (flag==0)
    {
        return b;
    }
    else
    {
        //unsigned M = b.BracketSize()+BracketSize();

        //Bracket result(M);

        //result=(*this)+b;
        //return result;
        return (*this)+b;
    }
}*/
