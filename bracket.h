#ifndef BRACKET_H
#define BRACKET_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cstdlib>
#include <ctime>

#include <math.h>

struct _GainPower {
	double g;             //коэффициент одночлена
    int p1;               //степень переменной ksi_1  
    int p2;               //степень переменной ksi_2  
    int p3;               //степень переменной ksi_3  
    int p4;               //степень переменной ksi_4  
};

struct _Power {
    int p1;               //степень переменной ksi_1  
    int p2;               //степень переменной ksi_2  
    int p3;               //степень переменной ksi_3  
    int p4;               //степень переменной ksi_4  
};

typedef struct _GainPower GainPower_t;
typedef struct _Power Power_t;


class Bracket
{
public:
    Bracket(unsigned N = 0);
    Bracket(std::vector<GainPower_t>& terms);

    virtual ~Bracket();
    void BracketInit(unsigned N);
	void BracketInit(std::vector<GainPower_t>& terms);
    void BracketCleanUp();
    unsigned BracketSize() {return m_N;}

    std::vector<double> GetGains();
    void SetGains(std::vector<double>& gains);

    std::vector<Power_t> GetPowers(); 
    void SetPowers(std::vector<Power_t>& powers);

	std::vector<GainPower_t> GetTerms();
    void SetTerms(std::vector<GainPower_t>& terms);

    void ShowElements();

    Bracket operator*(const Bracket& b);

    Bracket operator*(double number);

    Bracket& operator*=(const Bracket& b);

    Bracket operator+(const Bracket& b);

    Bracket& operator+=(const Bracket& b);

    Bracket operator-(const Bracket& b);

    Bracket& operator=(const Bracket& right);

private:
	std::vector<GainPower_t> m_Terms;  //скобка представляет собой вектор одночленов (m_Terms)

    unsigned m_N;
};

#endif // BRACKET_H
