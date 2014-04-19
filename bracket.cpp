#include "bracket.h"
#include "BracketFunctions.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cstdlib>
#include <ctime>

#include <math.h>


using namespace std;

Bracket::Bracket(unsigned N) : m_N(N)
{
    BracketInit(N);
}

/*Bracket::Bracket(std::vector<double>& gains, std::vector<Power_t>& powers)
{
    BracketInit(gains, powers);
}*/

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
        GainPower_t GP = {0.0, 0, 0, 0, 0};

        m_Terms.push_back(GP);
    }
}

/*void Bracket::BracketInit(std::vector<double>& gains, std::vector<Power_t>& powers)
{
    unsigned N=gains.size();
    m_N=N;
    for (unsigned i=0;i<m_N;i++)
    {
        m_Gains.push_back(gains[i]);

        m_Powers.push_back(powers[i]);
    }
}*/

void Bracket::BracketCleanUp()
{
    m_Terms.clear();
}

std::vector<double> Bracket::GetGains()
{
	std::vector<double> gains(m_N);
	for (unsigned i=0;i<m_N;i++)
		gains[i]=m_Terms[i].g;
	return gains;
}

std::vector<Power_t> Bracket::GetPowers()
{
	vector<Power_t> powers(m_N);
	Power_t pw={0, 0, 0, 0};
	for (unsigned i=0;i<m_N;i++)
	{
		powers[i].p1=m_Terms[i].p1;
		powers[i].p2=m_Terms[i].p2;
		powers[i].p3=m_Terms[i].p3;
		powers[i].p4=m_Terms[i].p4;
	}
	return powers;
}

void Bracket::SetGains(std::vector<double>& gains)
{
	unsigned N=gains.size();
	if (N == m_N)
		for (unsigned i=0;i<N;i++)
		{
			m_Terms[i].g=gains[i];
		}
	else
	{
		m_Terms.clear();
		GainPower_t term={0.0, 0, 0, 0, 0};
		for (unsigned i=0;i<N;i++)
		{
			term.g=gains[i];
			m_Terms.push_back(term);
		}
		m_N=N;
	}
}

std::vector<GainPower_t> Bracket::GetTerms()
{
	std::vector<GainPower_t> result;
	//copy(m_Terms.begin(),m_Terms.end(),result.begin());
	for (unsigned i=0;i<m_N;i++)
		result.push_back(m_Terms.at(i));
	return result;
}

void Bracket::SetTerms(std::vector<GainPower_t>& terms)
{
	m_N=terms.size();
	m_Terms.clear();
	for (unsigned i=0;i<m_N;i++)
		m_Terms.push_back(terms.at(i));
}

void Bracket::SetPowers(std::vector<Power_t>& powers)
{
    unsigned N=powers.size();
	if (N == m_N)
		for (unsigned i=0;i<N;i++)
		{
			m_Terms[i].p1=powers[i].p1;
			m_Terms[i].p2=powers[i].p2;
			m_Terms[i].p3=powers[i].p3;
			m_Terms[i].p4=powers[i].p4;
		}
	else
	{
		m_Terms.clear();
		GainPower_t term={0.0, 0, 0, 0, 0};
		for (unsigned i=0;i<N;i++)
		{
			term.p1=powers[i].p1;
			term.p2=powers[i].p2;
			term.p3=powers[i].p3;
			term.p4=powers[i].p4;
			m_Terms.push_back(term);
		}
		m_N=N;
	}
}

void Bracket::ShowElements()
{
    std::cout << "======================" << std::endl;
    std::cout << "m_N = " << m_N << std::endl;

	std::cout << "m_Terms = " << std::endl;

	for (unsigned i=0;i<m_N;i++)
	{
		std::cout << "  gain = "<< m_Terms.at(i).g  << "   ";
		std::cout << "  powers = "<< m_Terms.at(i).p1 << "   ";
		std::cout << "  "<< m_Terms.at(i).p2 << "   ";
		std::cout << "  "<< m_Terms.at(i).p3 << "   ";
		std::cout << "  "<< m_Terms.at(i).p4 << "   " << std::endl;
	}
    std::cout << "======================" << std::endl;
}

Bracket Bracket::operator*(const Bracket& b)
{
    unsigned M = b.m_N*m_N;
    unsigned L = b.m_N;
    unsigned N = m_N;

    Bracket result(M);

	unsigned offset;
	double Ai, Bj, Cij;

    for(unsigned i = 0; i < N; i++)
    {
        offset = i*L;

        Ai = m_Terms.at(i).g;

        for(unsigned j = 0; j < L; j++)
         {

            Bj = b.m_Terms.at(j).g;

            Cij = Ai * Bj;

			result.m_Terms.at(offset+j).g=Cij;
			result.m_Terms.at(offset+j).p1=m_Terms.at(i).p1+b.m_Terms.at(j).p1;
			result.m_Terms.at(offset+j).p2=m_Terms.at(i).p2+b.m_Terms.at(j).p2;
			result.m_Terms.at(offset+j).p3=m_Terms.at(i).p3+b.m_Terms.at(j).p3;
			result.m_Terms.at(offset+j).p4=m_Terms.at(i).p4+b.m_Terms.at(j).p4;
        }
    }

	SimplifyBracket(result);
    return result;
}


Bracket Bracket::operator*(double number)
{
    if (number == 0.0)
    {
       Bracket result(1);
       return result;
    }
    else
    {
        Bracket result = *this;
        for (unsigned i=0;i<m_N;i++)
			result.m_Terms.at(i).g *= number;

		SimplifyBracket(result);
        return result;
    }
}

Bracket& Bracket::operator*=(const Bracket& b)
{
    *this=(*this)*b;
    return *this;
}

Bracket Bracket::operator+(const Bracket& b)
{
	Bracket result = *this;
	for (unsigned i=0;i<b.m_N;i++)
		result.m_Terms.push_back(b.m_Terms.at(i));

	result.m_N=m_N+b.m_N;

	SimplifyBracket(result);
	return result;
}

Bracket Bracket::operator-(const Bracket& b)
{
    Bracket result = b;
	result = result*(-1.0);

	for (unsigned i=0;i<m_N;i++)
		result.m_Terms.push_back(m_Terms.at(i));

	result.m_N=m_N+b.m_N;

	SimplifyBracket(result);
	return result;
}

Bracket& Bracket::operator+=(const Bracket& b)
{
    *this=*this+b;
    return *this;
}

Bracket& Bracket::operator=(const Bracket& right)
{
    if (this == &right)
    {
        return *this;
    }
	else
	{
		m_N=right.m_N;

		m_Terms.clear();
		for (unsigned i=0;i<right.m_N;i++)
			m_Terms.push_back(right.m_Terms.at(i));

		return *this;
	}
}

/*Bracket::Bracket(Bracket& obj)
{
    for (unsigned i=0;i<(obj.m_Gains).size();i++)
        this->m_Gains.push_back((obj.m_Gains)[i]);

    for (unsigned i=0;i<(obj.m_Powers).size();i++)
        this->m_Powers.push_back((obj.m_Powers)[i]);
    this->m_N=obj.m_N;
}*/
