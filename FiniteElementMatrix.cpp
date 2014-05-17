#include "FiniteElementMatrix.h"
#include "bracket.h"
#include "VectFunctions.h"
#include "BracketFunctions.h"

#include <math.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>

using namespace std;

FiniteElementMatrix::FiniteElementMatrix(unsigned p, double** simplex_peaks, double** Eps, double** Mu) : m_P(p)
{
	m_Peaks = new double*[m_CountPeaks];
	m_MatrixMu = new double*[m_Dim];
	m_MatrixEps = new double*[m_Dim];

	for (unsigned i=0;i<m_CountPeaks;i++)
	{
		m_Peaks[i] = new double[m_Dim];
		if (i<m_Dim)
		{
			m_MatrixMu[i] = new double[m_Dim];
			m_MatrixEps[i] = new double[m_Dim];
		}
		for (unsigned j=0;j<m_Dim;j++)
		{
			m_Peaks[i][j]=simplex_peaks[i][j];
			if (i<m_Dim)
			{
				m_MatrixMu[i][j]=Mu[i][j];
				m_MatrixEps[i][j]=Eps[i][j];
			}
		}
	}
	//Инициализация массива факториалами
	//максимальный порядок 20; максимальный необходимый факториал 25

	for (unsigned i=0;i<m_MaxOrderFact+1;i++)
	{
		m_ArrayFact[i]=Fact(i);
	}
	//--------------------------
	//m_MatrixSize=(unsigned)((p+1)*(p+3)*(p+4)/2);

	m_MatrixSize=6*(p+1)+6*p*(p+1)+p*(p*p-1);

	//Инициализация унитарных базисных векторов
	std::vector<double> r1, r2, r3, r4;

	for (unsigned i=0;i<m_Dim;i++)
	{
		r1.push_back(m_Peaks[0][i]);
		r2.push_back(m_Peaks[1][i]);
		r3.push_back(m_Peaks[2][i]);
		r4.push_back(m_Peaks[3][i]);
	}
	for (unsigned i=0;i<m_Dim;i++)
	{
		m_UnVect_1.push_back(r1[i]-r4[i]);
		m_UnVect_2.push_back(r2[i]-r4[i]);
		m_UnVect_3.push_back(r3[i]-r4[i]);
	}

	//Вычисление якобиана
	m_Icob=ScalarProduct(VectProduct(m_UnVect_1,m_UnVect_2),m_UnVect_3);

	//Инициализация градиентов ksi
	m_GradKsi_1=VectProduct(m_UnVect_2,m_UnVect_3);
	for (unsigned i=0;i<m_Dim;i++)
		m_GradKsi_1[i]=m_GradKsi_1[i]/m_Icob;

	m_GradKsi_2=VectProduct(m_UnVect_3,m_UnVect_1);
	for (unsigned i=0;i<m_Dim;i++)
		m_GradKsi_2[i]=m_GradKsi_2[i]/m_Icob;

	m_GradKsi_3=VectProduct(m_UnVect_1,m_UnVect_2);
	for (unsigned i=0;i<m_Dim;i++)
		m_GradKsi_3[i]=m_GradKsi_3[i]/m_Icob;

	for (unsigned i=0;i<m_Dim;i++)
	{
		m_GradKsi_4.push_back(-m_GradKsi_1[i]-m_GradKsi_2[i]-m_GradKsi_3[i]);
	}

	//Чтение из файла корней полиномов Лежандра и весов

	double value;

	m_QuadOrder=3;
	string s_Roots="Roots3Nodes.txt";
	string s_Weights="Weights3Nodes.txt";
	Q3 = (int)pow(m_QuadOrder,3);
	Q2 = (int)pow(m_QuadOrder,2);

	// m_QuadOrder=15;
	// s_Roots="Roots.txt";
	// s_Weights="Weights.txt";

	ifstream tfile(s_Roots);
	while (!tfile.eof())
	{
		tfile>>value;
		m_Roots.push_back(value);
	}
	tfile.close();

	tfile.open(s_Weights);
	while (!tfile.eof())
	{
		tfile>>value;
		m_Weights.push_back(value);
	}
	tfile.close();

	FormArrayNum();
	FormArrayAnalyt();

	MatrixInit();
	ShowMatrixes();
	NumMatrixInit();

	CompareMatrixs();
}

FiniteElementMatrix::~FiniteElementMatrix()
{
	delete[] m_MetrMatrix;
	delete[] m_EulerMatrix;
	delete[] m_NumEulerMatrix;
	delete[] m_NumMetrMatrix;

	for (unsigned i=0; i < m_CountPeaks; i++)
		delete[] m_Peaks[i];

	for (unsigned i=0; i < m_Dim; i++)
	{
		delete[] m_MatrixEps[i];
		delete[] m_MatrixMu[i];
	}

	for (unsigned i=0; i < m_MatrixSize*Q3; i++)
	{
		delete[] m_Arr_AllNodes[i];
		delete[] m_Arr_RotAllNodes[i];
	}
}


void FiniteElementMatrix::findIndex(unsigned gamma, unsigned beta, unsigned order, unsigned& idxLow, unsigned& idxHigh)
{
	if ((gamma==order)||(beta==order))
	{
		idxLow=0;
		idxHigh=m_P;
	}
	else
	{
		idxLow=1;
		idxHigh=m_P+1;
	}
}

void FiniteElementMatrix::MatrixInit()
{
	//Выделяю память под матрицу (метрическая и Эйлера)

	m_MetrMatrix = new double[m_MatrixSize*m_MatrixSize];
	m_EulerMatrix = new double[m_MatrixSize*m_MatrixSize];

	//--------------------------------------------------
	Bracket euler_br, metr_br;

	for (unsigned i=0;i<m_MatrixSize;i++)
		for (unsigned j=0;j<m_MatrixSize;j++)
		{
			euler_br=GeneralVectorTensorVectorProduct(m_ArrAnalyt_EigFunc[i],m_ArrAnalyt_EigFunc[j],m_MatrixEps);

			*(m_EulerMatrix+m_MatrixSize*i+j)=Integrate(euler_br);

			metr_br=GeneralVectorTensorVectorProduct(m_ArrAnalyt_RotEigFunc[i],m_ArrAnalyt_RotEigFunc[j],m_MatrixMu);

			*(m_MetrMatrix+m_MatrixSize*i+j)=Integrate(metr_br);
		}
}//MatrixInit

void FiniteElementMatrix::AddToVectBracket(std::vector<GainPower_t>& terms, std::vector<Bracket>& vect_bracket)
{
	vect_bracket.push_back(Bracket(terms));
}


//-----------------------Определение вектора----------------------------------------------------------
std::vector<double> FiniteElementMatrix::DefVector(unsigned ind)
{
	switch(ind) {
	case 1: return m_GradKsi_1;
	case 2: return m_GradKsi_2;
	case 3: return m_GradKsi_3;
	case 4: return m_GradKsi_4;
	default: return m_GradKsi_1;
	}

}
//----------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------
//-------------Добавление многочлена Сильвестра-------------------------------------------------------
//--ind - индекс (i, j, k или l);
//--numb - число 1, 2, 3, 4
//если gamma или beta равны числу numb, то добавляемый многочлен Сильвестра будет несмещенным,
//иначе - смещенным
//если берется i-й многочлен, то numb=1. В этом случае будет многочлен по степени ksi_1 (p1).
//Если j-й многочлен, то numb=2, многочлен относительно степени ksi_2 (p2);
//numb=3 - p3; numb=4 - p4.
void FiniteElementMatrix::AddSilvester(unsigned gamma, unsigned beta, unsigned numb, unsigned ind,
									   std::vector<Bracket>& vect_bracket)
{
	GainPower_t term={1.0, 0, 0, 0, 0};
	std::vector<GainPower_t> terms;
	/*Power_t local_powers={0,0,0,0};
	std::vector<Power_t> powers;
	std::vector<double> gains;*/

	if ((gamma == numb)||(beta==numb)) //добавление несмещенного многочлена Сильвестра
	{
		if (ind==0)
		{
			if (vect_bracket.empty())
			{
				terms.push_back(term);
				AddToVectBracket(terms,vect_bracket);
			}
		}
		else
		{
			double g;
			for (unsigned s=0;s<=ind-1;s++)
			{
				g=(m_P+0.0)/(s+1.0);
				if (g!=0)
				{
					term.g = g;
					LocalTermsChange(term,numb,1);
					terms.push_back(term);
				}

				g=(-1.*s)/(s+1.0);
				if (g!=0)
				{
					term.g = g;
					LocalTermsChange(term,numb,0);
					terms.push_back(term);
				}

				AddToVectBracket(terms,vect_bracket);
				terms.clear();
			}
		}
	}
	else                    //добавление смещенного многочлена Сильвестра
	{
		if (ind==1)
		{
			if (vect_bracket.empty())
			{
				terms.push_back(term);
				AddToVectBracket(terms,vect_bracket);
			}
		}
		else
		{
			double g;
			for (unsigned s=0;s<=ind-2;s++)
			{
				g=(m_P+0.0)/(s+1.0)/(s+1.0);
				if (g!=0)
				{
					term.g = g;
					LocalTermsChange(term,numb,1);
					terms.push_back(term);
				}

				g=(-1.*s-1.0)/(s+1.0)/(s+1.0);
				if (g!=0)
				{
					term.g = g;
					LocalTermsChange(term,numb,0);
					terms.push_back(term);
				}

				AddToVectBracket(terms,vect_bracket);
				terms.clear();
			}
		}
	}
}
//----------------------------------------------------------------------------------------------
//-----------------------Вычисление факториала-------------------------------------------------
double FiniteElementMatrix::Fact(unsigned N)
{
	double f = 1.0;

	for (unsigned i=0;i<N;i++)
	{
		f *= (i+1);
	}

	return f;
}
//----------------------------------------------------------------------------------------------
inline double FiniteElementMatrix::CalcFact(unsigned N)
{
	return m_ArrayFact[N];
}
//--------------------------Интегрирование скобки-----------------------------------------------
double FiniteElementMatrix::Integrate (Bracket& br)
{
	double I=0;
	unsigned pow1, pow2, pow3, pow4;
	std::vector<GainPower_t> terms;
	terms=br.GetTerms();

	double i1, i2, i3, i4, i_s;
	for (unsigned i=0;i<terms.size();i++)
	{
		pow1=terms.at(i).p1;
		pow2=terms.at(i).p2;
		pow3=terms.at(i).p3;
		pow4=terms.at(i).p4;
		i1=CalcFact(pow1);
		i2=CalcFact(pow2);
		i3=CalcFact(pow3);
		i4=CalcFact(pow4);
		i_s=CalcFact(pow1+pow2+pow3+pow4+3);
		//I=I+terms.at(i).g*i1*i2*i3*i4*6.0/i_s;
		I=I+terms.at(i).g*i1*i2*i3*i4/i_s;
	}
	return I;
}

//----------Вычисление ротора от вектора, который задан в виде: скобка * (ksi_n*nabla(ksi_m) - ksi_m*nabla(ksi_n))
// n < m
std::vector<Bracket> FiniteElementMatrix::RotorCalc(Bracket& br, unsigned n, unsigned m)
{
	std::vector<Bracket> result;
	//Скалярная функция * градиент вектора

	std::vector<double> v_1=DefVector(m);
	std::vector<double> v_2=DefVector(n);
	std::vector<double> vect=VectProduct(v_2,v_1);
	for (unsigned i=0;i<m_Dim;i++)
	{
		result.push_back(br*(vect[i]*2.0));
	}

	//Градиент скалярной функции
	Bracket local_br(1);

	//Трехмерный вектор, состоящий из нулевых скобок
	std::vector<Bracket> vect_nabla_phi;
	for (unsigned i=0;i<m_Dim;i++)
		vect_nabla_phi.push_back(local_br);

	/*std::vector<Power_t> powers, current_powers;
	std::vector<double> gains, current_gains;
	Power_t pw;
	double gn;*/
	std::vector<GainPower_t> terms, current_terms;
	GainPower_t trm;

	/*current_powers=br.GetPowers();
	current_gains=br.GetGains();*/
	current_terms = br.GetTerms();

	for (unsigned i=0;i<br.BracketSize();i++)
	{
		if ((current_terms[i]).p1>0)
		{
			//Добавление градиента ksi_1 с соответствующим коэффициентом,
			//если текущий элемент скобки зависит от переменной ksi_1
			for (unsigned j=0;j<m_Dim;j++)
			{
				trm.p1=(unsigned)(current_terms[i].p1-1);
				trm.p2=current_terms[i].p2;
				trm.p3=current_terms[i].p3;
				trm.p4=current_terms[i].p4;

				trm.g=m_GradKsi_1[j]*current_terms[i].g*((double)(current_terms[i]).p1);

				terms.push_back(trm);

				local_br.SetTerms(terms);

				terms.clear();

				vect_nabla_phi[j]=vect_nabla_phi[j]+local_br;
			}
		}//ksi_1

		//Аналогично градиенты ksi_2, ksi_3, ksi_4
		if ((current_terms[i]).p2>0)
		{
			//Добавление градиента ksi_1 с соответствующим коэффициентом,
			//если текущий элемент скобки зависит от переменной ksi_1
			for (unsigned j=0;j<m_Dim;j++)
			{
				trm.p1=current_terms[i].p1;
				trm.p2=(unsigned)(current_terms[i].p2-1);;
				trm.p3=current_terms[i].p3;
				trm.p4=current_terms[i].p4;

				trm.g=m_GradKsi_2[j]*current_terms[i].g*((double)(current_terms[i]).p2);

				terms.push_back(trm);

				local_br.SetTerms(terms);

				terms.clear();

				vect_nabla_phi[j]=vect_nabla_phi[j]+local_br;
			}
		}//ksi_2

		if ((current_terms[i]).p3>0)
		{
			//Добавление градиента ksi_1 с соответствующим коэффициентом,
			//если текущий элемент скобки зависит от переменной ksi_1
			for (unsigned j=0;j<m_Dim;j++)
			{
				trm.p1=current_terms[i].p1;
				trm.p2=current_terms[i].p2;
				trm.p3=(unsigned)(current_terms[i].p3-1);;
				trm.p4=current_terms[i].p4;

				trm.g=m_GradKsi_3[j]*current_terms[i].g*((double)(current_terms[i]).p3);

				terms.push_back(trm);

				local_br.SetTerms(terms);

				terms.clear();

				vect_nabla_phi[j]=vect_nabla_phi[j]+local_br;
			}
		}//ksi_3

		if ((current_terms[i]).p4>0)
		{
			//Добавление градиента ksi_1 с соответствующим коэффициентом,
			//если текущий элемент скобки зависит от переменной ksi_1
			for (unsigned j=0;j<m_Dim;j++)
			{
				trm.p1=current_terms[i].p1;
				trm.p2=current_terms[i].p2;
				trm.p3=current_terms[i].p3;
				trm.p4=(unsigned)(current_terms[i].p4-1);

				trm.g=m_GradKsi_4[j]*current_terms[i].g*((double)(current_terms[i]).p4);

				terms.push_back(trm);

				local_br.SetTerms(terms);

				terms.clear();

				vect_nabla_phi[j]=vect_nabla_phi[j]+local_br;
			}
		}//ksi_4
	}//vect_nabla_phi готов

	trm.p1=0;
	trm.p2=0;
	trm.p3=0;
	trm.p4=0;
	trm.g=1.0;

	terms.push_back(trm);

	local_br.SetTerms(terms);

	//(ksi_n*nabla(ksi_m) - ksi_m*nabla(ksi_n)) = vect_ksi_nm
	std::vector<Bracket> vect_ksi_nm;

	LocalTermsChange(trm,n,1);
	terms.push_back(trm);
	local_br.SetTerms(terms);

	for (unsigned i=0;i<m_Dim;i++)
		vect_ksi_nm.push_back(local_br*v_1[i]);

	LocalTermsChange(trm,n,0);
	LocalTermsChange(trm,m,1);

	terms.clear();
	terms.push_back(trm);
	local_br.SetTerms(terms);

	Bracket br_1=local_br*v_2[0];
	vect_ksi_nm[0]=vect_ksi_nm[0]-br_1;
	br_1=local_br*v_2[1];
	vect_ksi_nm[1]=vect_ksi_nm[1]-br_1;
	br_1=local_br*v_2[2];
	vect_ksi_nm[2]=vect_ksi_nm[2]-br_1;//vect_ksi_nm

	std::vector<Bracket> res_1;
	VectBracketProduct(vect_nabla_phi,vect_ksi_nm,res_1);

	for (unsigned i=0;i<m_Dim;i++)
		result[i]=result[i]+res_1[i];

	/*for (unsigned i=0;i<m_Dim;i++)
	{
	((Bracket)(result[i])).ShowElements();
	}*/

	return result;
}

void FiniteElementMatrix::ShowMatrixes()
{
	std::cout << "============Analytical EulerMatrix==============" << std::endl;
	for (unsigned i=0;i<m_MatrixSize;i++)
	{
		for (unsigned j=0;j<m_MatrixSize;j++)
		{
			std::cout<< " " << *(m_EulerMatrix+i*m_MatrixSize+j);

		}
		std::cout << std::endl;
	}

	std::cout << "============Analytical MetrMatrix==============" << std::endl;
	for (unsigned i=0;i<m_MatrixSize;i++)
	{
		for (unsigned j=0;j<m_MatrixSize;j++)
		{
			std::cout<< " " << *(m_MetrMatrix+i*m_MatrixSize+j);

		}
		std::cout << std::endl;
	}
	system("pause");

}

std::vector<Bracket> FiniteElementMatrix::FormVectEigFunc(Bracket& br, unsigned n, unsigned m)
{
	std::vector<Bracket> result;
	std::vector<double> a=DefVector(m);
	std::vector<double> b=DefVector(n);

	Bracket local_br(1);

	GainPower_t trm={1.0, 0, 0, 0, 0};
	std::vector<GainPower_t> terms;
	LocalTermsChange(trm,n,1);
	terms.push_back(trm);
	local_br.SetTerms(terms);

	Bracket res_br=local_br*(a[0]);
	result.push_back(res_br);

	res_br=local_br*(a[1]);
	result.push_back(res_br);

	res_br=local_br*(a[2]);
	result.push_back(res_br);

	terms.clear();

	LocalTermsChange(trm,n,0);
	LocalTermsChange(trm,m,1);

	terms.push_back(trm);

	local_br.SetTerms(terms);
	for (unsigned i=0;i<m_Dim;i++)
	{
		res_br=local_br*b[i];
		result[i]=result[i]-res_br;
	}

	for (unsigned i=0;i<m_Dim;i++)
	{
		result[i]=result[i]*br;
	}

	return result;
}

//-------------------------Формирование матрицы Эйлера и метрической матрицы с использованием численного интегрирования-------
void FiniteElementMatrix::NumMatrixInit()
{
	//Выделяю память под матрицу Эйлера и метрическую матрицу
	m_NumEulerMatrix = new double[m_MatrixSize*m_MatrixSize];
	m_NumMetrMatrix = new double[m_MatrixSize*m_MatrixSize];
	//--------------------------------------------------

	double I, R, elem; 
	double u, v;

	for (unsigned i=0;i<m_MatrixSize;i++)
	{
		for (unsigned j=0;j<m_MatrixSize;j++)
		{
			I=0;
			R=0;
			for (unsigned j_u=0; j_u<m_QuadOrder;j_u++)
			{
				for (unsigned j_v=0; j_v<m_QuadOrder;j_v++)
				{
					for (unsigned j_w=0; j_w<m_QuadOrder;j_w++)
					{
						u=m_Roots[j_u];
						v=m_Roots[j_v];

						elem=NumericalVectorTensorVectorProduct(m_Arr_AllNodes[i*Q3+j_u*Q2+j_v*m_QuadOrder+j_w],
																m_Arr_AllNodes[j*Q3+j_u*Q2+j_v*m_QuadOrder+j_w],
																m_MatrixEps);

						I+=elem*m_Weights[j_u]*m_Weights[j_v]*m_Weights[j_w]*pow(u,2.0)*v;

						elem=NumericalVectorTensorVectorProduct(m_Arr_RotAllNodes[i*Q3+j_u*Q2+j_v*m_QuadOrder+j_w],
																m_Arr_RotAllNodes[j*Q3+j_u*Q2+j_v*m_QuadOrder+j_w],
																m_MatrixMu);

						R+=elem*m_Weights[j_u]*m_Weights[j_v]*m_Weights[j_w]*pow(u,2.0)*v;
					}//j_w
				}//j_v
			}//j_u
		I*=pow(0.5,3);
		R*=pow(0.5,3);
		*(m_NumEulerMatrix+m_MatrixSize*i+j)=I;
		*(m_NumMetrMatrix+m_MatrixSize*i+j)=R;
		}
	}
}//NumMatrixInit

void FiniteElementMatrix::NumFormVectEigFunc(std::vector<Bracket>& vect_br, unsigned n, unsigned m, 
											 double ksi1, double ksi2, double ksi3, 
											 std::vector<double>& result)
{
	result.clear();
	double numb=DefKsi(n,ksi1,ksi2,ksi3);
	result=DefVector(m);
	MultNumber(result,numb);

	vector<double> vect=DefVector(n);
	numb=(-1.0)*DefKsi(m,ksi1,ksi2,ksi3);
	MultNumber(vect,numb);

	SumVector(result,vect);

	numb=ProdVectBracketValue(vect_br, ksi1, ksi2, ksi3);
	MultNumber(result,numb);
}

void FiniteElementMatrix::CompareMatrixs()
{
	std::cout << "============EulerMatrix(i,j) - NumEulerMatrix(i,j)==============" << std::endl;
	double max=0;
	for (unsigned i=0;i<m_MatrixSize;i++)
	{
		for (unsigned j=0;j<m_MatrixSize;j++)
		{
			std::cout<< "  " << *(m_EulerMatrix+i*m_MatrixSize+j) - *(m_NumEulerMatrix+i*m_MatrixSize+j);
			if (max < (*(m_EulerMatrix+i*m_MatrixSize+j) - *(m_NumEulerMatrix+i*m_MatrixSize+j)))
				max = *(m_EulerMatrix+i*m_MatrixSize+j) - *(m_NumEulerMatrix+i*m_MatrixSize+j);
		}
		std::cout << std::endl;
	}
	std::cout << "============Max value of difference==============" << std::endl;
	std::cout << max << endl;

	std::cout << "============MetrMatrix(i,j) - NumMetrMatrix(i,j)==============" << std::endl;
	max=0;
	for (unsigned i=0;i<m_MatrixSize;i++)
	{
		for (unsigned j=0;j<m_MatrixSize;j++)
		{
			std::cout<< "  " << *(m_MetrMatrix+i*m_MatrixSize+j) - *(m_NumMetrMatrix+i*m_MatrixSize+j);
			if (max < (*(m_MetrMatrix+i*m_MatrixSize+j) - *(m_NumMetrMatrix+i*m_MatrixSize+j)))
				max = *(m_MetrMatrix+i*m_MatrixSize+j) - *(m_NumMetrMatrix+i*m_MatrixSize+j);
		}
		std::cout << std::endl;
	}
	std::cout << "============Max value of difference==============" << std::endl;
	std::cout << max << endl;
}
//------------------------------Массив m_ArrAnalytEigFunc заполняется элементами------------------------
//------------элементами массива являются векторные собственные функции---------------------------------
void FiniteElementMatrix::FormArrayAnalyt()
{
	unsigned ilow, ihigh,
		jlow, jhigh,
		klow, khigh,
		llow, lhigh;

	std::vector<Bracket> vect_bracket;

	unsigned n1, m1;
	vector<unsigned> nm1;

	Bracket cur_bracket;

	std::vector<Bracket> eig_func, rot_eig_func;

	for (unsigned gamma=1;gamma<=4;gamma++)
	{
		for (unsigned beta = gamma+1;beta<=4;beta++)
		{
			findIndex(gamma, beta, 1, ilow, ihigh);
			findIndex(gamma, beta, 2, jlow, jhigh);
			findIndex(gamma, beta, 3, klow, khigh);
			findIndex(gamma, beta, 4, llow, lhigh);

			for (unsigned i=ilow;i<=ihigh;i++)
			{
				for (unsigned j=jlow;j<=jhigh;j++)
				{
					for (unsigned k=klow;k<=khigh;k++)
					{
						for (unsigned l=llow;l<=lhigh;l++)
						{
							if ((i+j+k+l)==(m_P+2))
							{
								//В список добавляются 4 многочлена Сильвестра
								//В список будут добавляться 4 многочлена Сильвестра
								AddSilvester(gamma, beta, 1, i, vect_bracket);
								AddSilvester(gamma, beta, 2, j, vect_bracket);
								AddSilvester(gamma, beta, 3, k, vect_bracket);
								AddSilvester(gamma, beta, 4, l, vect_bracket);

								Def_nm(gamma,beta,nm1);
								n1=nm1[0];
								m1=nm1[1];

								//Формирую скобку cur_bracket, а именно перемножаю
								//скобки массива vect_bracket

								cur_bracket=vect_bracket[0];

								for (unsigned counter=1; counter<vect_bracket.size();counter++)
								{
									cur_bracket=cur_bracket*vect_bracket[counter];
								}

								eig_func=FormVectEigFunc(cur_bracket,n1,m1);

								rot_eig_func=RotorCalc(cur_bracket,n1,m1);

								m_ArrAnalyt_EigFunc.push_back(eig_func);
								m_ArrAnalyt_RotEigFunc.push_back(rot_eig_func);
							}//условие
							nm1.clear();
							vect_bracket.clear();
						}//l
					}//k
				}//j
			}//i
		}//beta
	}//gamma

	//---------------------Вывод на консоль---------------------------------
	/*cout << "Array of eigfunctions:" << endl;
	cout << endl;

	for (unsigned i=0;i<m_MatrixSize;i++)
	{
	cout << "Eigfunction " << "[" << i << "] = " << endl;
	for (unsigned j=0;j<m_Dim;j++)
	m_ArrAnalyt_EigFunc[i][j].ShowElements();
	}

	cout << "Array of eigfunctions rotors" << endl;
	cout << endl;
	for (unsigned i=0;i<m_MatrixSize;i++)
	{
	cout << "Rotors of eigfunction " << "[" << i << "] = " << endl;
	for (unsigned j=0;j<m_Dim;j++)
	m_ArrAnalyt_RotEigFunc[i][j].ShowElements();
	}*/
}

//----------------Массив m_Arr_AllNodes заполняется элементами---------------------------
//----------------Элементами являются значения собственных функций, вычисленные во------
// всех узлах кубатурных формул. Массив 3 (=m_Dim) строки и m_MatrixSize*m_QuadOrder. Сначала
//идет первая собственная функция, вычисленная во всех узлах, затем вторая и т. д.
// Код метода сильно повторяет код метода FormArrAnalytical.
void FiniteElementMatrix::FormArrayNum()
{
	//Выделение памяти
	m_Arr_AllNodes = new double*[m_MatrixSize*Q3];
	m_Arr_RotAllNodes = new double*[m_MatrixSize*Q3];
	for (unsigned i=0;i<m_MatrixSize*Q3;i++)
	{
		m_Arr_AllNodes[i] = new double[m_Dim];
		m_Arr_RotAllNodes[i] = new double[m_Dim];
	}
	//----------------------------------------------------
	double u, v, w;
	double x, y, z;

	unsigned ilow, ihigh,
		jlow, jhigh,
		klow, khigh,
		llow, lhigh;

	std::vector<Bracket> vect_bracket;    

	vector<Bracket> vect_RotEigFunc;       // ротор векторной функции
	//--------------------Числовые векторы, соответсвующие собственным функциям и их роторам--------------------
	vector<double> num_EigFunc;

	vector<double> num_RotEigFunc;
	//---------------------------------------------------------------------------------------------------------
	Bracket cur_bracket;

	unsigned n1, m1;
	vector<unsigned> nm1;

	int elems_count=0;

	for (unsigned gamma=1;gamma<=4;gamma++)
	{
		for (unsigned beta = gamma+1;beta<=4;beta++)
		{
			findIndex(gamma, beta, 1, ilow, ihigh);
			findIndex(gamma, beta, 2, jlow, jhigh);
			findIndex(gamma, beta, 3, klow, khigh);
			findIndex(gamma, beta, 4, llow, lhigh);

			for (unsigned i=ilow;i<=ihigh;i++)
			{
				for (unsigned j=jlow;j<=jhigh;j++)
				{
					for (unsigned k=klow;k<=khigh;k++)
					{
						for (unsigned l=llow;l<=lhigh;l++)
						{
							if ((i+j+k+l)==(m_P+2))
							{
								//В список добавляются 4 многочлена Сильвестра
								//В список будут добавляться 4 многочлена Сильвестра
								AddSilvester(gamma, beta, 1, i, vect_bracket);
								AddSilvester(gamma, beta, 2, j, vect_bracket);
								AddSilvester(gamma, beta, 3, k, vect_bracket);
								AddSilvester(gamma, beta, 4, l, vect_bracket);

								Def_nm(gamma,beta,nm1);
								n1=nm1[0];
								m1=nm1[1];

								//Формирую скобку cur_bracket для текущего элемента, а именно перемножаю
								//скобки массива vect_bracket

								cur_bracket=vect_bracket[0];

								for (unsigned counter=1; counter<vect_bracket.size();counter++)
								{
									cur_bracket=cur_bracket*vect_bracket[counter];
								}

								vect_RotEigFunc=RotorCalc(cur_bracket,n1,m1);

								for (unsigned j_u=0; j_u<m_QuadOrder;j_u++)
								{
									for (unsigned j_v=0; j_v<m_QuadOrder;j_v++)
									{
										for (unsigned j_w=0; j_w<m_QuadOrder;j_w++)
										{
											u=m_Roots[j_u];
											v=m_Roots[j_v];
											w=m_Roots[j_w];
											x = u*v*w;
											y = u*v*(1.0 - w);
											z = u*(1.0 - v);

											NumFormVectEigFunc(vect_bracket, n1, m1, x, y, z, num_EigFunc);
											VectBracketValue(vect_RotEigFunc, x, y, z, num_RotEigFunc);

											for (unsigned d=0;d<m_Dim;d++)
											{
												m_Arr_AllNodes[elems_count*Q3+j_u*Q2+j_v*m_QuadOrder+j_w][d]=num_EigFunc[d];
												m_Arr_RotAllNodes[elems_count*Q3+j_u*Q2+j_v*m_QuadOrder+j_w][d]=num_RotEigFunc[d];
											}//d
										}//j_w
									}//j_v
								}//j_u
							}//проверка условия (i+j+k+l==m_P+2)
							elems_count++;
							nm1.clear();
							vect_bracket.clear();
						}//l
					}//k
				}//j
			}//i
		}//beta
	}//gamma
	display(m_MatrixSize*Q3,m_Dim,m_Arr_AllNodes);
}//FormArrayNum

void FiniteElementMatrix::display(unsigned rows, unsigned columns, double** arr)
{
	for (unsigned i=0;i<rows;i++)
	{
		for (unsigned j=0;j<columns; j++)
		{
			cout << arr[i][j] << "   ";
		}
		cout << endl;
	}
	system("pause");
}
