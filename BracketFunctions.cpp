#include "bracket.h"
#include "BracketFunctions.h"

#include <vector>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <algorithm>
#include <iterator>

using namespace std;

void VectBracketValue (vector<Bracket>& br, vector<double>& vect, double ksi1, double ksi2, double ksi3)
{
	vect.clear();
	std::vector<double> gains;
	std::vector<Power_t> powers;
	double S;
	for (unsigned i=0;i<br.size();i++)
	{
		S=0;
		gains=br[i].GetGains();
		powers=br[i].GetPowers();
		for (unsigned j=0;j<br[i].BracketSize();j++)
		{
			S=S+gains[j]*pow(ksi1,powers[j].p1)*
				pow(ksi2,powers[j].p2)*
				pow(ksi3,powers[j].p3)*
				pow((1-ksi1-ksi2-ksi3),powers[j].p4);
		}
		vect.push_back(S);

		gains.clear();
		powers.clear();
	}
}

//----------В переменной local_powers присвоить полю с номером ind значение value--------------
void LocalPowersChange(Power_t& local_powers, unsigned ind, unsigned value)
{
	if (ind==1)
	{
		local_powers.p1=value;
	}
	if (ind==2)
	{
		local_powers.p2=value;
	}
	if (ind==3)
	{
		local_powers.p3=value;
	}
	if (ind==4)
	{
		local_powers.p4=value;
	}
}
//--------------------------------------------------------------------------------------------
//------Исходная скобка упрощается: сортируется по возрастанию степеней, удаляются слагаемые
//------с нулевыми коэффициентами, коэффициенты при одинаковых степенях суммируются-----------
bool comparefun(GainsPowers_t x, GainsPowers_t y)
{
	return ( (x.p1*1000+x.p2*100+x.p3*10+x.p4) > (y.p1*1000+y.p2*100+y.p3*10+y.p4) );
}

bool equality(GainsPowers_t x, GainsPowers_t y)
{
	return ( (x.p1*1000+x.p2*100+x.p3*10+x.p4) == (y.p1*1000+y.p2*100+y.p3*10+y.p4) );
}
bool cond(GainsPowers_t x)
{
	return (x.g==0.0);
}
void SimplifyBracket(Bracket& br)
{
	vector<double> gains=br.GetGains();
	vector<Power_t> powers=br.GetPowers();
	vector<GainsPowers_t> GP;
	GainsPowers_t unit_gp={0, 0, 0, 0, 0};
	for (unsigned i=0;i<br.BracketSize();i++)
	{
		unit_gp.g=gains[i];
		unit_gp.p1=powers[i].p1;
		unit_gp.p2=powers[i].p2;
		unit_gp.p3=powers[i].p3;
		unit_gp.p4=powers[i].p4;
		GP.push_back(unit_gp);
	}
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
	vector<GainsPowers_t>::iterator new_end;
	new_end=remove_if(GP.begin(),GP.end(),cond);
	GP.erase(new_end, GP.end());
	
//	copy(GP.begin(),GP.end(),ostream_iterator<int> (cout," "));
	cout << "!";
}