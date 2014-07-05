#include "FiniteElementMatrix.h"
#include "bracket.h"
#include "VectFunctions.h"
#include "BracketFunctions.h"
#include "InterplNode.h"

#include <math.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>

#include <tchar.h>
#include <windows.h>

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
	//--------------------------------------------------------------------------------

	//FormArrayAnalyt();
	//FormArrayNum();

	FormArrayAnalyt_LinComb();
	FormArrayNum_LinComb();

	MatrixInit();
	NumMatrixInit();

	ShowMatrixes();	
	CompareMatrixs();
	//---------------------------------------------------------------------------------
	//Визуализация получившихся векторных полей
	string s1 = "F:\\TestBracket\\Array of tetrahedron nodes\\Ksi1Data.txt";
	string s2 = "F:\\TestBracket\\Array of tetrahedron nodes\\Ksi2Data.txt";
	string s3 = "F:\\TestBracket\\Array of tetrahedron nodes\\Ksi3Data.txt";
	ExportFuncValues(s1, s2, s3);
}

FiniteElementMatrix::~FiniteElementMatrix()
{
	//delete[] m_MetrMatrix;
	delete[] m_EulerMatrix;
	delete[] m_NumEulerMatrix;
	//delete[] m_NumMetrMatrix;

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
		//delete[] m_Arr_RotAllNodes[i];
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
	LARGE_INTEGER Frequency, StartPerformCount, StopPerformCount;
	int bHighRes = QueryPerformanceFrequency (&Frequency);
	QueryPerformanceCounter (&StartPerformCount);

	//-----------------------------------------------------------------------------------------------
	//Выделяю память под матрицу (метрическая и Эйлера)

	m_MetrMatrix = new double[m_MatrixSize*m_MatrixSize];
	m_EulerMatrix = new double[m_MatrixSize*m_MatrixSize];

	//--------------------------------------------------
	int N_reserve = 36*pow(m_P,2)*1000;
	Bracket euler_br(1, N_reserve), metr_br(1, N_reserve);   //в эти скобки будет записываться результат произведения собственная функция 
	// на тензор на собсвенную функцию 

	Bracket br_sum(1, N_reserve), br_prod(1, N_reserve); //вспомогательные скобки

	for (unsigned i=0;i<m_MatrixSize;i++)
		for (unsigned j=0;j<m_MatrixSize;j++)
		{
			GeneralVectorTensorVectorProduct(&m_ArrAnalyt_EigFunc[i],&m_ArrAnalyt_EigFunc[j],m_MatrixEps,
				&euler_br, &br_sum, &br_prod);

			*(m_EulerMatrix+m_MatrixSize*i+j)=Integrate(&euler_br);


			// метрическая матрица
			/*
			GeneralVectorTensorVectorProduct(&m_ArrAnalyt_RotEigFunc[i],&m_ArrAnalyt_RotEigFunc[j],m_MatrixMu,
			&metr_br, &br_sum, &br_prod);

			*(m_MetrMatrix+m_MatrixSize*i+j)=Integrate(&metr_br);
			*/
		}
		//-------------------------------------------------------------------------------------------

		QueryPerformanceCounter (&StopPerformCount);
		double msTime = (double)(StopPerformCount.QuadPart - StartPerformCount.QuadPart) / (double)Frequency.QuadPart * 1.E3;
		cout << "MatrixInit: ellapsed time = " << msTime << endl;
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
				g=(m_P+2+0.0)/(s+1.0);     
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
				g=(m_P+2+0.0)/(s+1.0);     // change!       было деление на 2 скобки
				if (g!=0)
				{
					term.g = g;
					LocalTermsChange(term,numb,1);
					terms.push_back(term);
				}

				g=(-1.*s-1.0)/(s+1.0);     // change!     аналогично
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
double FiniteElementMatrix::Integrate (Bracket* br)
{
	double I=0;
	unsigned pow1, pow2, pow3, pow4;
	std::vector<GainPower_t>* terms;
	terms=br->GetTermsPtr();

	unsigned N = terms->size();

	double i1, i2, i3, i4, i_s;
	for (unsigned i=0;i<N;i++)
	{
		pow1=terms->at(i).p1;
		pow2=terms->at(i).p2;
		pow3=terms->at(i).p3;
		pow4=terms->at(i).p4;
		i1=CalcFact(pow1);
		i2=CalcFact(pow2);
		i3=CalcFact(pow3);
		i4=CalcFact(pow4);
		i_s=CalcFact(pow1+pow2+pow3+pow4+3);
		//I=I+terms.at(i).g*i1*i2*i3*i4*6.0/i_s;
		I=I+terms->at(i).g*i1*i2*i3*i4/i_s;
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

	/*
	std::cout << "============Analytical MetrMatrix==============" << std::endl;
	for (unsigned i=0;i<m_MatrixSize;i++)
	{
	for (unsigned j=0;j<m_MatrixSize;j++)
	{
	std::cout<< " " << *(m_MetrMatrix+i*m_MatrixSize+j);

	}
	std::cout << std::endl;
	} */
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
	LARGE_INTEGER Frequency, StartPerformCount, StopPerformCount;
	int bHighRes = QueryPerformanceFrequency (&Frequency);
	QueryPerformanceCounter (&StartPerformCount);
	//-------------------------------------------------------------------------------------------------
	//Выделяю память под матрицу Эйлера и метрическую матрицу
	m_NumEulerMatrix = new double[m_MatrixSize*m_MatrixSize];
	// m_NumMetrMatrix = new double[m_MatrixSize*m_MatrixSize];
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

						/*elem=NumericalVectorTensorVectorProduct(m_Arr_RotAllNodes[i*Q3+j_u*Q2+j_v*m_QuadOrder+j_w],
							m_Arr_RotAllNodes[j*Q3+j_u*Q2+j_v*m_QuadOrder+j_w],
							m_MatrixMu);

						R+=elem*m_Weights[j_u]*m_Weights[j_v]*m_Weights[j_w]*pow(u,2.0)*v;*/
					}//j_w
				}//j_v
			}//j_u
			I*=pow(0.5,3);
			R*=pow(0.5,3);
			*(m_NumEulerMatrix+m_MatrixSize*i+j)=I;
			//*(m_NumMetrMatrix+m_MatrixSize*i+j)=R;
		}
	}
	//-------------------------------------------------------------------------------------------
	QueryPerformanceCounter (&StopPerformCount);
	double msTime = (double)(StopPerformCount.QuadPart - StartPerformCount.QuadPart) / (double)Frequency.QuadPart * 1.E3;

	cout << "NumMatrixInit: ellapsed time = " << msTime << endl;
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
	double max=0;
	std::cout << "============EulerMatrix(i,j) - NumEulerMatrix(i,j)==============" << std::endl;
	for (unsigned i=0;i<m_MatrixSize;i++)
	{
		for (unsigned j=0;j<m_MatrixSize;j++)
		{
			//std::cout<< "  " << *(m_EulerMatrix+i*m_MatrixSize+j) - *(m_NumEulerMatrix+i*m_MatrixSize+j);
			if (max < abs( (*(m_EulerMatrix+i*m_MatrixSize+j) - *(m_NumEulerMatrix+i*m_MatrixSize+j)) ) )
				max = abs(*(m_EulerMatrix+i*m_MatrixSize+j) - *(m_NumEulerMatrix+i*m_MatrixSize+j));
		}
		//std::cout << std::endl;
	}
	std::cout << "============Max value of difference==============" << std::endl;
	std::cout << max << endl;

	/*max=0;
	std::cout << "============MetrMatrix(i,j) - NumMetrMatrix(i,j)==============" << std::endl;
	for (unsigned i=0;i<m_MatrixSize;i++)
	{
		for (unsigned j=0;j<m_MatrixSize;j++)
		{
			//std::cout<< "  " << *(m_MetrMatrix+i*m_MatrixSize+j) - *(m_NumMetrMatrix+i*m_MatrixSize+j);
			if (max < abs( (*(m_MetrMatrix+i*m_MatrixSize+j) - *(m_NumMetrMatrix+i*m_MatrixSize+j)) ) )
				max = abs(*(m_MetrMatrix+i*m_MatrixSize+j) - *(m_NumMetrMatrix+i*m_MatrixSize+j));
		}
		//std::cout << std::endl;
	}
	std::cout << "============Max value of difference==============" << std::endl;
	std::cout << max << endl;*/
}

//------------------------------Массив m_ArrAnalytEigFunc заполняется элементами------------------------
//------------элементами массива являются векторные собственные функции---------------------------------
void FiniteElementMatrix::FormArrayAnalyt()
{
	LARGE_INTEGER Frequency, StartPerformCount, StopPerformCount;
	int bHighRes = QueryPerformanceFrequency (&Frequency);
	QueryPerformanceCounter (&StartPerformCount);
	//-------------------------------------------------------------------------------------------------
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
	//-------------------------------------------------------------------------------------------
	QueryPerformanceCounter (&StopPerformCount);
	double msTime = (double)(StopPerformCount.QuadPart - StartPerformCount.QuadPart) / (double)Frequency.QuadPart * 1.E3;

	cout << "FormArrayAnalyt: ellapsed time = " << msTime << endl;
}

//Другой способ формирования массива всех базисных функций при помощи линейных комбинаций.
//линейно-зависимые базисные функции исключаются
void FiniteElementMatrix::FormArrayAnalyt_LinComb()
{
	unsigned low = 0;
	unsigned high = m_P + 1;

	unsigned flag=0;
	unsigned gamma=0, beta=0;
	unsigned* index_array;
	unsigned* non_zero_arr;

	index_array = new unsigned[3];
	non_zero_arr = new unsigned[3];

	int N_reserve = 36*pow(m_P,2)*1000;

	Bracket cur_bracket;
	Bracket sum1_bracket(1,N_reserve), sum2_bracket(1,N_reserve); 

	std::vector<Bracket> vect_bracket;
	std::vector<Bracket> eig_func_1, eig_func_2, eig_func_3;
	vector<Bracket> eig_func;

	unsigned n1, m1;
	vector<unsigned> nm1;

	for (unsigned i = low;i <= high;i++)
	{
		for (unsigned j = low;j <= high;j++)
		{
			for (unsigned k = low;k <= high;k++)
			{
				for (unsigned l = low;l <= high;l++)
				{
					if ((i+j+k+l)==(m_P+2))
					{
						//определяю, где расположен элемент: на ребре, грани или 
						//во внутренней области
						WhereIsNode(i, j, k, l, flag);
						if (flag==1)
						{
							EdgeIndDefine(i, j, k, l, gamma, beta);
							//стандартная процедура формирования собственной функции
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
							m_ArrAnalyt_EigFunc.push_back(eig_func);
							nm1.clear();
							vect_bracket.clear();
						}//ребро
						//////////////////////////////////////////////////////////////////////////////////
						if (flag==2)
						{
							DefFaceInd(i,j,k,l,index_array,non_zero_arr);
							//добавляю многочлены Сильвестра. Чтобы не думать о gamma и beta
							//с учетом того
							//все многочлены будут несмещенными, один из
							//индексов i,j,k,l = 0, он увеличивается на 1. 
							if (i==0)
								AddSilvester(0, 0, 1, i+1, vect_bracket);
							else
								AddSilvester(0, 0, 1, i, vect_bracket);
							if (j==0)
								AddSilvester(0, 0, 2, j+1, vect_bracket);
							else
								AddSilvester(0, 0, 2, j, vect_bracket);
							if (k==0)
								AddSilvester(0, 0, 3, k+1, vect_bracket);
							else
								AddSilvester(0, 0, 3, k, vect_bracket);
							if (l==0)
								AddSilvester(0, 0, 4, l+1, vect_bracket);
							else
								AddSilvester(0, 0, 4, l, vect_bracket);
							//Формирую скобку cur_bracket, а именно перемножаю
							//скобки массива vect_bracket

							cur_bracket=vect_bracket[0]*(m_P+2.0);

							for (unsigned counter=1; counter<vect_bracket.size();counter++)
							{
								cur_bracket=cur_bracket*vect_bracket[counter];
							}

							FormFaceFunc(&cur_bracket, index_array,
								non_zero_arr, &eig_func, 1);
							m_ArrAnalyt_EigFunc.push_back(eig_func);
							FormFaceFunc(&cur_bracket, index_array,
								non_zero_arr, &eig_func, 2);
							m_ArrAnalyt_EigFunc.push_back(eig_func);

							vect_bracket.clear();
						}//грань
						//////////////////////////////////////////////////////////////////////////////
						if (flag==3)
						{
							//В список будут добавляться 4 многочлена Сильвестра
							AddSilvester(0, 0, 1, i, vect_bracket);
							AddSilvester(0, 0, 2, j, vect_bracket);
							AddSilvester(0, 0, 3, k, vect_bracket);
							AddSilvester(0, 0, 4, l, vect_bracket);

							cur_bracket=vect_bracket[0]*pow(m_P+2.0,2);

							for (unsigned counter=1; counter<vect_bracket.size();counter++)
							{
								cur_bracket=cur_bracket*vect_bracket[counter];
							}

							index_array[0] = i;
							index_array[1] = j;
							index_array[2] = k;
							index_array[3] = l;

							FormInsideFunc(&cur_bracket, index_array, &eig_func,1);
							m_ArrAnalyt_EigFunc.push_back(eig_func);
							FormInsideFunc(&cur_bracket, index_array, &eig_func,2);
							m_ArrAnalyt_EigFunc.push_back(eig_func);
							FormInsideFunc(&cur_bracket, index_array, &eig_func,3);
							m_ArrAnalyt_EigFunc.push_back(eig_func);
						}//внутренняя область
					}//условие
				}//l
			}//k
		}//j
	}//i
	m_MatrixSize = (unsigned) (6*(m_P+1) + m_P*(m_P+1) + m_P*(m_P*m_P-1)/2);   //размерности метрической матрицы
	// и матрицы Эйлера, если базинсые функции формируются как линейные комбинации
}

// Формирование числового массива базисных функций через линейные комбинации
// сначала базисная функция получается аналитически, затем в нее подставляеются
// нужные значения переменных ksi_1, ksi_2, ksi_3, в результате чего получается
// соответсвующая числовая базисная функция
// !!! Эту функцию можно запускать только после того, как отработал метод FormArrayAnalyt_LinComb()
void FiniteElementMatrix::FormArrayNum_LinComb()
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
	//--------------------Числовые векторы, соответсвующие собственным функциям и их роторам--------------------
	vector<double> num_EigFunc;

	vector<double> num_RotEigFunc;
	//---------------------------------------------------------------------------------------------------------

	for (unsigned elems_count = 0; elems_count < m_MatrixSize; elems_count++)
	{
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

					VectBracketValue(m_ArrAnalyt_EigFunc.at(elems_count), x, y, z, num_EigFunc);

					for (unsigned d=0;d<m_Dim;d++)
					{
						m_Arr_AllNodes[elems_count*Q3+j_u*Q2+j_v*m_QuadOrder+j_w][d]=num_EigFunc[d];
						// m_Arr_RotAllNodes[elems_count*Q3+j_u*Q2+j_v*m_QuadOrder+j_w][d]=num_RotEigFunc[d];
					}//d
				}//j_w
			}//j_v
		}//j_u
	}// elems_count
}

void FiniteElementMatrix::FormFaceFunc(Bracket* bracket,unsigned* index_array,
									   unsigned* non_zero_array, std::vector<Bracket>* eig_func, unsigned ind)
{
	GainPower_t trm = {1, 0, 0, 0, 0};
	vector<GainPower_t> terms0, terms1, terms2;
	vector<double> a;
	vector<double> b;
	vector<double> c;
	if (ind==1)
	{
		a = DefVector(index_array[0]);
		b = DefVector(index_array[1]);
		c = DefVector(index_array[2]);

		LocalTermsChange(trm, index_array[1], 1);
		LocalTermsChange(trm, index_array[2], 1);
		trm.g =  1.0/(non_zero_array[1]+0.0) + 1.0/(non_zero_array[2]+0.0);            
		terms0.push_back(trm);
		terms1.push_back(trm);
		terms2.push_back(trm);
		//обнулили trm
		LocalTermsChange(trm, index_array[1], 0);
		LocalTermsChange(trm, index_array[2], 0);

		LocalTermsChange(trm, index_array[0], 1);
		LocalTermsChange(trm, index_array[2], 1);
		trm.g =  -1.0/(non_zero_array[2]+0.0);
		terms0.push_back(trm);
		terms1.push_back(trm);
		terms2.push_back(trm);
		//обнулили trm
		LocalTermsChange(trm, index_array[1], 0);
		LocalTermsChange(trm, index_array[2], 0);

		LocalTermsChange(trm, index_array[0], 1);
		LocalTermsChange(trm, index_array[1], 1);
		trm.g =  -1.0/(non_zero_array[1]+0.0);
		terms0.push_back(trm);
		terms1.push_back(trm);
		terms2.push_back(trm);
	}// 1
	if (ind==2)
	{
		a = DefVector(index_array[1]);
		b = DefVector(index_array[0]);
		c = DefVector(index_array[2]);

		LocalTermsChange(trm, index_array[0], 1);
		LocalTermsChange(trm, index_array[2], 1);
		trm.g =  1.0/(non_zero_array[0]+0.0) + 1.0/(non_zero_array[2]+0.0);            
		terms0.push_back(trm);
		terms1.push_back(trm);
		terms2.push_back(trm);
		//обнулили trm
		LocalTermsChange(trm, index_array[0], 0);
		LocalTermsChange(trm, index_array[2], 0);

		LocalTermsChange(trm, index_array[1], 1);
		LocalTermsChange(trm, index_array[2], 1);
		trm.g =  -1.0/(non_zero_array[2]+0.0);
		terms0.push_back(trm);
		terms1.push_back(trm);
		terms2.push_back(trm);
		//обнулили trm
		LocalTermsChange(trm, index_array[1], 0);
		LocalTermsChange(trm, index_array[2], 0);

		LocalTermsChange(trm, index_array[0], 1);
		LocalTermsChange(trm, index_array[1], 1);
		trm.g =  -1.0/(non_zero_array[0]+0.0);
		terms0.push_back(trm);
		terms1.push_back(trm);
		terms2.push_back(trm);
	}// 2
	terms0.at(0).g = terms0.at(0).g*a.at(0);
	terms0.at(1).g = terms0.at(1).g*b.at(0);
	terms0.at(2).g = terms0.at(2).g*c.at(0);
	Bracket br0(terms0);

	terms1.at(0).g = terms1.at(0).g*a.at(1);
	terms1.at(1).g = terms1.at(1).g*b.at(1);
	terms1.at(2).g = terms1.at(2).g*c.at(1);
	Bracket br1(terms1);

	terms2.at(0).g = terms2.at(0).g*a.at(2);
	terms2.at(1).g = terms2.at(1).g*b.at(2);
	terms2.at(2).g = terms2.at(2).g*c.at(2);
	Bracket br2(terms2);

	Mult(bracket,&br0,&eig_func->at(0));
	Mult(bracket,&br1,&eig_func->at(1));
	Mult(bracket,&br2,&eig_func->at(2));
}

void FiniteElementMatrix::FormInsideFunc(Bracket* bracket,unsigned* node_array,
										 std::vector<Bracket>* eig_func, unsigned ind)
										 // node_array - массив [i, j, k, l]										 
{
	GainPower_t trm = {1, 0, 0, 0, 0};
	vector<GainPower_t> terms0, terms1, terms2;
	// сначала идет коэффициент перед градиентом ksi_1, потом перед ksi_2, потом перед ksi_3
	// затем перед ksi_4. 
	if (ind==1)
	{
		trm.p2 = 1;
		trm.p3 = 1;
		trm.p4 = 1;
		trm.g =  1.0/(node_array[1]+0.0)/(node_array[3]+0.0) + 1.0/(node_array[2]+0.0)/(node_array[3]+0.0)
			+ 1.0/(node_array[1]+0.0)/(node_array[2]+0.0);
		terms0.push_back(trm);
		terms1.push_back(trm);
		terms2.push_back(trm);

		trm.p1 = 1;
		trm.p2 = 0;
		trm.g = - 1.0/(node_array[2]+0.0)/(node_array[3]+0.0);
		terms0.push_back(trm);
		terms1.push_back(trm);
		terms2.push_back(trm);

		trm.p3 = 0;
		trm.p2 = 1;
		trm.g = - 1.0/(node_array[1]+0.0)/(node_array[3]+0.0);
		terms0.push_back(trm);
		terms1.push_back(trm);
		terms2.push_back(trm);

		trm.p3=1;
		trm.p4=0;
		trm.g= - 1.0/(node_array[1]+0.0)/(node_array[2]+0.0);
		terms0.push_back(trm);
		terms1.push_back(trm);
		terms2.push_back(trm);
	}
	if (ind==2)
	{
		trm.p2 = 1;
		trm.p3 = 1;
		trm.p4 = 1;
		trm.g =  -1.0/(node_array[2]+0.0)/(node_array[3]+0.0);
		terms0.push_back(trm);
		terms1.push_back(trm);
		terms2.push_back(trm);

		trm.p1 = 1;
		trm.p2 = 0;
		trm.g = 1.0/(node_array[0]+0.0)/(node_array[3]+0.0) + 1.0/(node_array[2]+0.0)/(node_array[3]+0.0)+
			1.0/(node_array[0]+0.0)/(node_array[2]+0.0);
		terms0.push_back(trm);
		terms1.push_back(trm);
		terms2.push_back(trm);

		trm.p3 = 0;
		trm.p2 = 1;
		trm.g = - 1.0/(node_array[0]+0.0)/(node_array[3]+0.0);
		terms0.push_back(trm);
		terms1.push_back(trm);
		terms2.push_back(trm);

		trm.p3=1;
		trm.p4=0;
		trm.g= - 1.0/(node_array[0]+0.0)/(node_array[2]+0.0);
		terms0.push_back(trm);
		terms1.push_back(trm);
		terms2.push_back(trm);
	}
	if (ind==3)
	{
		trm.p2 = 1;
		trm.p3 = 1;
		trm.p4 = 1;
		trm.g =  -1.0/(node_array[1]+0.0)/(node_array[3]+0.0);
		terms0.push_back(trm);
		terms1.push_back(trm);
		terms2.push_back(trm);

		trm.p1 = 1;
		trm.p2 = 0;
		trm.g = - 1.0/(node_array[0]+0.0)/(node_array[3]+0.0);
		terms0.push_back(trm);
		terms1.push_back(trm);
		terms2.push_back(trm);

		trm.p3 = 0;
		trm.p2 = 1;
		trm.g = 1.0/(node_array[1]+0.0)/(node_array[3]+0.0) + 1.0/(node_array[0]+0.0)/(node_array[3]+0.0) +
			1.0/(node_array[0]+0.0)/(node_array[2]+0.0);
		terms0.push_back(trm);
		terms1.push_back(trm);
		terms2.push_back(trm);

		trm.p3=1;
		trm.p4=0;
		trm.g= - 1.0/(node_array[0]+0.0)/(node_array[1]+0.0);
		terms0.push_back(trm);
		terms1.push_back(trm);
		terms2.push_back(trm);
	}

	terms0.at(0).g = terms0.at(0).g*m_GradKsi_1.at(0);
	terms0.at(1).g = terms0.at(1).g*m_GradKsi_2.at(0);
	terms0.at(2).g = terms0.at(2).g*m_GradKsi_3.at(0);
	terms0.at(3).g = terms0.at(3).g*m_GradKsi_4.at(0);
	Bracket br0(terms0);

	terms1.at(0).g = terms1.at(0).g*m_GradKsi_1.at(1);
	terms1.at(1).g = terms1.at(1).g*m_GradKsi_2.at(1);
	terms1.at(2).g = terms1.at(2).g*m_GradKsi_3.at(1);
	terms1.at(3).g = terms1.at(3).g*m_GradKsi_4.at(1);
	Bracket br1(terms1);

	terms2.at(0).g = terms2.at(0).g*m_GradKsi_1.at(2);
	terms2.at(1).g = terms2.at(1).g*m_GradKsi_2.at(2);
	terms2.at(2).g = terms2.at(2).g*m_GradKsi_3.at(2);
	terms2.at(3).g = terms2.at(3).g*m_GradKsi_4.at(2);
	Bracket br2(terms2);

	Mult(bracket,&br0,&eig_func->at(0));
	Mult(bracket,&br1,&eig_func->at(1));
	Mult(bracket,&br2,&eig_func->at(2));
}



//----------------Массив m_Arr_AllNodes заполняется элементами---------------------------
//----------------Элементами являются значения собственных функций, вычисленные во------
// всех узлах кубатурных формул. Массив 3 (=m_Dim) строки и m_MatrixSize*m_QuadOrder. Сначала
//идет первая собственная функция, вычисленная во всех узлах, затем вторая и т. д.
// Код метода сильно повторяет код метода FormArrAnalytical.
void FiniteElementMatrix::FormArrayNum()
{
	LARGE_INTEGER Frequency, StartPerformCount, StopPerformCount;
	int bHighRes = QueryPerformanceFrequency (&Frequency);
	QueryPerformanceCounter (&StartPerformCount);
	//-------------------------------------------------------------------------------------------------
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
								elems_count++;
								nm1.clear();
								vect_bracket.clear();
							}//проверка условия (i+j+k+l==m_P+2)
						}//l
					}//k
				}//j
			}//i
		}//beta
	}//gamma
	//display(m_MatrixSize*Q3,m_Dim,m_Arr_AllNodes);
	//-------------------------------------------------------------------------------------------
	QueryPerformanceCounter (&StopPerformCount);
	double msTime = (double)(StopPerformCount.QuadPart - StartPerformCount.QuadPart) / (double)Frequency.QuadPart * 1.E3;

	cout << "FormArrayNum: ellapsed time = " << msTime << endl;
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

//В файл записываются значения базисных функций, вычисленных в точках тетраэдра. Массив
//точек задается в 3-х файлах, их имена передаются в качестве параметров s1, s2, s3. В каждом файле
//находятся массивы значений ksi1, ksi2, ksi3
void FiniteElementMatrix::ExportFuncValues(string s1, string s2, string s3)
{
	// Вычисление значений собственных функций в точках тетраэдра
	vector<double> arr_ksi_1, arr_ksi_2, arr_ksi_3;
	double value;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
	ifstream ksi_1_file(s1);
	while (!ksi_1_file.eof())
	{
		ksi_1_file>>value;
		arr_ksi_1.push_back(value);
	}

	ifstream ksi_2_file(s2);
	while (!ksi_2_file.eof())
	{
		ksi_2_file>>value;
		arr_ksi_2.push_back(value);
	}

	ifstream ksi_3_file(s3);
	while (!ksi_3_file.eof())
	{
		ksi_3_file>>value;
		arr_ksi_3.push_back(value);
	}

	ofstream func_file("F:\\TestBracket\\Array of basis functions values\\Values_1_face.txt");

	unsigned Size = arr_ksi_1.size();

	vector<double> points_vect;

	//for (unsigned i0 = 0;i0 < m_MatrixSize;i0++)
	//{
	for (unsigned i=0;i<Size - 1;i++)
	{
		VectBracketValue(m_ArrAnalyt_EigFunc[2], arr_ksi_1[i], arr_ksi_2[i], arr_ksi_3[i], points_vect);
		func_file << points_vect[0] << " ";
		func_file << points_vect[1] << " ";
		func_file << points_vect[2] << " ";
		func_file << endl;
	}
	//}

	ksi_1_file.close();
	ksi_2_file.close();
	ksi_3_file.close();

	func_file.close();
}
