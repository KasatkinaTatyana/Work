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

Bracket::Bracket(std::vector<GainPower_t>& terms)
{
	BracketInit(terms);
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
		GainPower_t GP = {0.0, 0, 0, 0, 0};

		m_Terms.push_back(GP);
	}
}

void Bracket::BracketInit(std::vector<GainPower_t>& terms)
{
	unsigned N=terms.size();
	m_N=N;
	m_Terms.clear();
	for (unsigned i=0;i<m_N;i++)
	{
		m_Terms.push_back(terms.at(i));
	}
}

void Bracket::BracketCleanUp()
{
	m_Terms.clear();
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

void Bracket::SetTermsPtr(std::vector<GainPower_t>* terms)
{
	m_N=terms->size();
	m_Terms.clear();

	//copy(terms->begin(), terms->begin() + m_N, inserter(m_Terms, m_Terms.begin()));
	for (unsigned i=0;i<m_N;i++)
		m_Terms.push_back(terms->at(i));
}


void Bracket::ShowElements()
{
	std::cout << "==================================================================" << std::endl;
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
	std::cout << "==================================================================" << std::endl;
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

Bracket::Bracket(const Bracket& obj)
{
	if (this != &obj)
	{
		m_Terms.clear();
		for (unsigned i=0;i<obj.m_N;i++)
			m_Terms.push_back((obj.m_Terms)[i]);

		m_N=obj.m_N;
	}
	//cout << "Copy!" << endl;
}

void Mult(Bracket* a, Bracket* b, Bracket* result)
{
	unsigned M = a->BracketSize()*b->BracketSize();
	unsigned L = a->BracketSize();
	unsigned N = b->BracketSize();
	
	//Bracket *result = new Bracket(M);

	unsigned offset;
	double Ai, Bj, Cij;

	std::vector<GainPower_t>* aTerms = a->Terms();
	std::vector<GainPower_t>* bTerms = b->Terms();
	std::vector<GainPower_t>* rTerms = result->Terms();

	for(unsigned i = 0; i < L; i++)
	{
		offset = i*N;

		Ai = aTerms->at(i).g;

		for(unsigned j = 0; j < N; j++)
		{

			Bj = bTerms->at(j).g;

			Cij = Ai * Bj;

			rTerms->at(offset+j).g=Cij;
			rTerms->at(offset+j).p1=aTerms->at(i).p1+bTerms->at(j).p1;
			rTerms->at(offset+j).p2=aTerms->at(i).p2+bTerms->at(j).p2;
			rTerms->at(offset+j).p3=aTerms->at(i).p3+bTerms->at(j).p3;
			rTerms->at(offset+j).p4=aTerms->at(i).p4+bTerms->at(j).p4;
		}
	}

	SimplifyBracketPtr(result);
}
