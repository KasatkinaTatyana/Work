#include "bracket.h"

Bracket::Bracket(unsigned N) : m_N(N)
{
    BracketInit(N);
}

Bracket::Bracket(std::vector<double> gains, std::vector<Power_t> powers)
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

void Bracket::BracketInit(std::vector<double> gains, std::vector<Power_t> powers)
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


