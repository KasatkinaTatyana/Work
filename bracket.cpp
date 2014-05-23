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

const GainPower_t Bracket::zero_term = {0.0, 0, 0, 0, 0};
const unsigned Bracket::unity = 1;

Bracket::Bracket(unsigned N) : m_N(N)
{
	BracketInit(N);
}

Bracket::Bracket(unsigned N, unsigned M)
{
	m_Terms.reserve(M);

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
		m_Terms.assign(N, zero_term);
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
	m_N=0;
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

//Текущая скобка превращается в скобку с нулевыми элементами размера s
void Bracket::SetSizePtr(unsigned* s)
{
	m_Terms.clear();
	m_Terms.assign(*s,zero_term);
	m_N = *s;
}
//К текущей скобка прибавляется скобка added
void Bracket::Plus(Bracket* added)
{
	unsigned L=added->BracketSize();
	m_N = m_N + L;
	std::vector<GainPower_t>* aTerms = added->GetTermsPtr();
	copy(aTerms->begin(), aTerms->end(), inserter(m_Terms,m_Terms.end()));
	SimplifyBracketPtr(this);
}

// В переменную result записывается значение a*b. Выражения вида (x, y, x) в таком варианте вычисляются некорректно
// все переменные a, b, result должны быть разными объектами
void Mult(Bracket* a, Bracket* b, Bracket* result)
{
	unsigned L = a->BracketSize();
	unsigned N = b->BracketSize();
	unsigned M = L*N;
	
	//Bracket *result = new Bracket(M);
	result->SetSizePtr(&M);

	unsigned offset;
	double Ai, Bj, Cij;

	std::vector<GainPower_t>* aTerms = a->GetTermsPtr();
	std::vector<GainPower_t>* bTerms = b->GetTermsPtr();
	std::vector<GainPower_t>* rTerms = result->GetTermsPtr();

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

void Plus(Bracket* a, Bracket* b, Bracket* result)
{
	unsigned M = a->BracketSize() + b->BracketSize();
	unsigned L = a->BracketSize();
	unsigned N = b->BracketSize();
	
	//Bracket *result = new Bracket(M);
	result->SetSizePtr(&M);

	std::vector<GainPower_t>* aTerms = a->GetTermsPtr();
	std::vector<GainPower_t>* bTerms = b->GetTermsPtr();
	std::vector<GainPower_t>* rTerms = result->GetTermsPtr();

	for (unsigned i=0;i<L;i++)
	{
		rTerms->at(i).g=aTerms->at(i).g;
		rTerms->at(i).p1=aTerms->at(i).p1;
		rTerms->at(i).p2=aTerms->at(i).p2;
		rTerms->at(i).p3=aTerms->at(i).p3;
		rTerms->at(i).p4=aTerms->at(i).p4;
	}

	for (unsigned i=0;i<N;i++)
	{
		rTerms->at(i+L).g=bTerms->at(i).g;
		rTerms->at(i+L).p1=bTerms->at(i).p1;
		rTerms->at(i+L).p2=bTerms->at(i).p2;
		rTerms->at(i+L).p3=bTerms->at(i).p3;
		rTerms->at(i+L).p4=bTerms->at(i).p4;
	}
	SimplifyBracketPtr(result);
}

void Mult(Bracket* br, double* numb, Bracket* result)
{
	unsigned N=br->BracketSize();

	result->SetSizePtr(&N);

	std::vector<GainPower_t>* brTerms = br->GetTermsPtr();
	std::vector<GainPower_t>* rTerms = result->GetTermsPtr();


	for (unsigned i=0;i<N;i++)
	{
		rTerms->at(i).g=brTerms->at(i).g*(*numb);
		rTerms->at(i).p1=brTerms->at(i).p1;
		rTerms->at(i).p2=brTerms->at(i).p2;
		rTerms->at(i).p3=brTerms->at(i).p3;
		rTerms->at(i).p4=brTerms->at(i).p4;
	}
}

