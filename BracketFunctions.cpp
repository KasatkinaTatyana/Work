
#include "bracket.h"
#include "BracketFunctions.h"

#include <vector>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <algorithm>
#include <iterator>

using namespace std;
//-----Возвращает значение произведения скобок при соответсвующих значениях ksi1, ksi2, ksi3------
double ProdVectBracketValue (vector<Bracket>& br, double ksi1, double ksi2, double ksi3)
{
	double result=1.0;
	double S;

	for (unsigned i=0;i<br.size();i++)
	{
		S = BracketValue(br[i], ksi1, ksi2, ksi3);
		result*=S;
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

void SimplifyBracketPtr(Bracket* br)
{
	vector<GainPower_t>* GP = br->GetTermsPtr();

	sort(GP->begin(),GP->end(),comparefun); //сортировка массива
	//удаление слагаемых с нулевыми коэффициентами

	//Слагаемые с одинаковыми степенями суммируются

	unsigned flag=0, position=0; 
	double s=0;

	for (unsigned i=0;i<(GP->size()-1);i++)
	{
		if (equality(GP->at(i),GP->at(i+1)))
		{
			s=s+GP->at(i).g;
			if (flag==1)
				GP->at(i).g=0;
			if (flag==0)
				position=i;
			flag=1;
		}
		else
		{
			if (flag==1)
			{
				s=s+GP->at(i).g;
				GP->at(i).g=0;
				GP->at(position).g=s;
				flag=0;
				s=0;
			}
		}
	}
	if (flag==1)
	{
		s=s+GP->at(GP->size()-1).g;
		GP->at(GP->size()-1).g=0;
		GP->at(position).g=s;
	}
	//удаление слагаемых с нулевыми коэффициентами
	vector<GainPower_t>::iterator new_end;
	new_end=remove_if(GP->begin(),GP->end(),cond);
	GP->erase(new_end, GP->end());

	//Проверка на то, чтобы в итоге скобка содержала хотя бы один элемент
/*
	if (GP->size() > 0)
		br->SetTerms(GP);
	else
	{
		GainPower_t t = {0.0, 0, 0, 0, 0};
		std::vector<GainPower_t> term;
		term.push_back(t);

		br->SetTerms(term);
	}
*/
	if (GP->empty())
	{
		GainPower_t t = {0.0, 0, 0, 0, 0};
		GP->push_back(t);
	}
	br->SetBracketSize(GP->size());
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
//------------Возвращает значение числовое скобки при ksi1, ksi2, ksi3 ------------------
double BracketValue(Bracket& br, double ksi1, double ksi2, double ksi3)
{
	vector<GainPower_t> terms;
	double S=0.0;
	terms=br.GetTerms();

	for (unsigned j=0;j<br.BracketSize();j++)
	{
		S=S+terms[j].g*pow(ksi1,terms[j].p1)*
			pow(ksi2,terms[j].p2)*
			pow(ksi3,terms[j].p3)*
			pow((1.0-ksi1-ksi2-ksi3),terms[j].p4);
	}
	return S;
}
//----В векторе result возвращается числовое значения вектора, компонентами кот. являются 
// скобки при соответствующих значениях ksi1, ksi2, ksi3
void VectBracketValue(std::vector<Bracket>& br, double ksi1, double ksi2, double ksi3, std::vector<double>& result)
{
	result.clear();
	for (unsigned i=0;i<br.size();i++)
		result.push_back(BracketValue(br[i], ksi1, ksi2, ksi3));
}

//-----------------Определение индексов n, m если заданы gamma и beta---------------------------------
// эта процедура используется для формирования базисных функций
void Def_nm(unsigned gamma, unsigned beta, vector<unsigned>& vect)
{
	unsigned n, m;
	if ((gamma==1)&&(beta==2))
	{
		n=3;
		m=4;
	}
	if ((gamma==2)&&(beta==1))
	{
		n=4;
		m=3;
	}
	if ((gamma==1)&&(beta==3))
	{
		n=2;
		m=4;
	}
	if ((gamma==3)&&(beta==1))
	{
		n=4;
		m=2;
	}
	if ((gamma==1)&&(beta==4))
	{
		n=2;
		m=3;
	}
	if ((gamma==4)&&(beta==1))
	{
		n=3;
		m=2;
	}
	if ((gamma==2)&&(beta==3))
	{
		n=1;
		m=4;
	}
	if ((gamma==2)&&(beta==4))
	{
		n=1;
		m=3;
	}
	if ((gamma==4)&&(beta==2))
	{
		n=3;
		m=1;
	}
	if ((gamma==3)&&(beta==4))
	{
		n=1;
		m=2;
	}
	if ((gamma==4)&&(beta==3))
	{
		n=2;
		m=1;
	}
	vect.push_back(n);
	vect.push_back(m);
}
//возвращает вектор, vect[0]=n, vect[1]=m; n<m;
//----------------------------------------------------------------------------------------------------

