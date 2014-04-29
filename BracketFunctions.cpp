
#include "bracket.h"
#include "BracketFunctions.h"

#include <vector>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <algorithm>
#include <iterator>

using namespace std;

double ProdVectBracketValue (vector<Bracket>& br, double ksi1, double ksi2, double ksi3)
{
	std::vector<GainPower_t> terms;

	double result=1.0;
	double S;

	for (unsigned i=0;i<br.size();i++)
	{
		S=0;
		terms=br[i].GetTerms();

		for (unsigned j=0;j<br[i].BracketSize();j++)
		{
			S=S+terms[j].g*pow(ksi1,terms[j].p1)*
				pow(ksi2,terms[j].p2)*
				pow(ksi3,terms[j].p3)*
				pow((1.0-ksi1-ksi2-ksi3),terms[j].p4);
		}
		result*=S;

		terms.clear();
	}
	return result;
}
//----------В переменной local_terms присвоить полю p с индексом ind значение value---------------
void LocalTermsChange(GainPower_t& local_terms, unsigned ind,unsigned value)
{
	if (ind==1)
	{
		local_terms.p1=value;
	}
	if (ind==2)
	{
		local_terms.p2=value;
	}
	if (ind==3)
	{
		local_terms.p3=value;
	}
	if (ind==4)
	{
		local_terms.p4=value;
	}
}
//--------------------------------------------------------------------------------------------
//------Исходная скобка упрощается: сортируется по возрастанию степеней, удаляются слагаемые
//------с нулевыми коэффициентами, коэффициенты при одинаковых степенях суммируются-----------
bool comparefun(GainPower_t x, GainPower_t y)
{
	return ( (x.p1*1000+x.p2*100+x.p3*10+x.p4) > (y.p1*1000+y.p2*100+y.p3*10+y.p4) );
}

bool equality(GainPower_t x, GainPower_t y)
{
	return ( (x.p1*1000+x.p2*100+x.p3*10+x.p4) == (y.p1*1000+y.p2*100+y.p3*10+y.p4) );
}
bool cond(GainPower_t x)
{
	return (x.g==0.0);
}
void SimplifyBracket(Bracket& br)
{
	vector<GainPower_t> GP = br.GetTerms();

	sort(GP.begin(),GP.end(),comparefun); //сортировка массива
	//удаление слагаемых с нулевыми коэффициентами

	//Слагаемые с одинаковыми степенями суммируются

	unsigned flag=0, position=0; 
	double s=0;

	for (unsigned i=0;i<(GP.size()-1);i++)
	{
		if (equality(GP[i],GP[i+1]))
		{
			s=s+GP[i].g;
			if (flag==1)
				GP[i].g=0;
			if (flag==0)
				position=i;
			flag=1;
		}
		else
		{
			if (flag==1)
			{
				s=s+GP[i].g;
				GP[i].g=0;
				GP[position].g=s;
				flag=0;
				s=0;
			}
		}
	}
	if (flag==1)
	{
		s=s+GP[GP.size()-1].g;
		GP[GP.size()-1].g=0;
		GP[position].g=s;
	}
	//удаление слагаемых с нулевыми коэффициентами
	vector<GainPower_t>::iterator new_end;
	new_end=remove_if(GP.begin(),GP.end(),cond);
	GP.erase(new_end, GP.end());

	//Проверка на то, чтобы в итоге скобка содержала хотя бы один элемент
	if (GP.size() > 0)
		br.SetTerms(GP);
	else
	{
		GainPower_t t = {0.0, 0, 0, 0, 0};
		std::vector<GainPower_t> term;
		term.push_back(t);

		br.SetTerms(term);
	}
}

//----------------------Возвращает ksi c номером n--------------------------------------
double DefKsi(unsigned n, double ksi1, double ksi2, double ksi3)
{
	double result;
	if (n==1)
		result=ksi1;
	if (n==2)
		result=ksi2;
	if (n==3)
		result=ksi3;
	if (n==4)
		result=1-ksi1-ksi2-ksi3;
	return result;
}