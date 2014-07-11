#ifndef BRACKET_H
#define BRACKET_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cstdlib>
#include <ctime>

#include <math.h>

struct _GainPower 
{
	double g;             //коэффициент одночлена
    int p1;               //степень переменной ksi_1  
    int p2;               //степень переменной ksi_2  
    int p3;               //степень переменной ksi_3  
    int p4;               //степень переменной ksi_4  
};

/*struct _Power {
    int p1;               //степень переменной ksi_1  
    int p2;               //степень переменной ksi_2  
    int p3;               //степень переменной ksi_3  
    int p4;               //степень переменной ksi_4  
};*/

typedef struct _GainPower GainPower_t;
//typedef struct _Power Power_t;


class Bracket
{
public:
	static const GainPower_t zero_term;
	static const unsigned unity;

    Bracket(unsigned N = 0);
    Bracket(std::vector<GainPower_t>& terms);
	Bracket::Bracket(unsigned N, unsigned M);   //конструктор, который резервирует M позиций под вектор, 
	//а инициализирует вектор из N компонентов

    virtual ~Bracket();

    void BracketInit(unsigned N);
	void BracketInit(std::vector<GainPower_t>& terms);
    void BracketCleanUp();
    unsigned BracketSize() {return m_N;}
	void SetBracketSize(unsigned N) {m_N = N;}

	std::vector<GainPower_t> GetTerms();

	GainPower_t& TermsElement(unsigned pos)
	{
		return m_Terms.at(pos);
	}

	std::vector<GainPower_t>* GetTermsPtr()
	{
		return &m_Terms;
	}

    void SetTerms(std::vector<GainPower_t>& terms);
	void SetTermsPtr(std::vector<GainPower_t>* terms);

    void ShowElements();

    Bracket operator*(const Bracket& b);

    Bracket operator*(double number);

    Bracket& operator*=(const Bracket& b);

    Bracket operator+(const Bracket& b);

    Bracket& operator+=(const Bracket& b);

    Bracket operator-(const Bracket& b);

    Bracket& operator=(const Bracket& right);

	Bracket(const Bracket& obj);

	void SetSizePtr(unsigned* s);

	void Plus(Bracket* added);

	void test_ref(std::vector<Bracket>& vect_br);


private:
	std::vector<GainPower_t> m_Terms;  //скобка представляет собой вектор одночленов (m_Terms)

    unsigned m_N;
};

void Mult(Bracket* a, Bracket* b, Bracket* result); 
void Plus(Bracket* a, Bracket* b, Bracket* result);
void Mult(Bracket* br, double* numb, Bracket* result);


typedef std::vector<Bracket> bp_t;
typedef std::vector<bp_t> bpp_t;

/*void f(bp_t *p)
{
	p->at(0);
}

void f1(bpp_t *p)
{
	bp_t pv = p->at(0);
	pv.at(0);
}*/

#endif // BRACKET_H
