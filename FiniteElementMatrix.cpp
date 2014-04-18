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

FiniteElementMatrix::FiniteElementMatrix(unsigned p,
										 double simplex_peaks[m_CountPeaks][m_Dim],
										 double Eps[m_Dim][m_Dim],
										 double Mu[m_Dim][m_Dim])
										 : m_P(p)
{
	for (unsigned i=0;i<m_CountPeaks;i++)
	{
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

	ifstream tfile("Roots.txt");
	while (!tfile.eof())
	{
		tfile>>value;
		m_Roots.push_back(value);
	}
	tfile.close();

	tfile.open("Weights.txt");
	while (!tfile.eof())
	{
		tfile>>value;
		m_Weights.push_back(value);
	}
	tfile.close();

	MatrixInit();
	ShowMatrixs();
}

FiniteElementMatrix::~FiniteElementMatrix()
{
	delete[] m_MetrMatrix;
	delete[] m_EulerMatrix;
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

	unsigned ilow, ihigh,
		jlow, jhigh,
		klow, khigh,
		llow, lhigh;

	unsigned in_ilow, in_ihigh,
		in_jlow, in_jhigh,
		in_klow, in_khigh,
		in_llow, in_lhigh;

	unsigned count_rows=0;
	unsigned count_columns=0;

	vector<Bracket> vectCur_EigFunc;
	vector<Bracket> vectIn_EigFunc;

	std::vector<Bracket> m_CurVectBracket;    //эта скобка соответсвует текущему элементу и будет храниться
	//при проходе по всем остальным элементам
	std::vector<Bracket> m_InnerVectBracket;  //эта скобка соответствует второму элементу (внутренний цикл)

	vector<Bracket> vectCur_RotEigFunc;
	vector<Bracket> vectIn_RotEigFunc;

	Bracket cur_bracket(1);
	Bracket inner_bracket(1);

	Bracket euler_br(1);
	Bracket metr_br(1);

	unsigned n1, m1, n2, m2;
	vector<unsigned> nm1, nm2;

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
								AddSilvester(gamma, beta, 1, i, m_CurVectBracket);
								AddSilvester(gamma, beta, 2, j, m_CurVectBracket);
								AddSilvester(gamma, beta, 3, k, m_CurVectBracket);
								AddSilvester(gamma, beta, 4, l, m_CurVectBracket);

								//Формирую скобку cur_bracket для текущего элемента, а именно перемножаю
								//скобки массива m_CurVectBracket

								cur_bracket=m_CurVectBracket[0];

								for (unsigned counter=1; counter<m_CurVectBracket.size();counter++)
								{
									cur_bracket=cur_bracket*m_CurVectBracket[counter];
								}

								Def_nm(gamma,beta,nm1);
								n1=nm1[0];
								m1=nm1[1];

								vectCur_EigFunc=FormVectEigFunc(cur_bracket,n1,m1);
								vectCur_RotEigFunc=RotorCalc(cur_bracket,n1,m1);

								//теперь запускаю внутрение циклы

								for (unsigned in_gamma=1;in_gamma<=4;in_gamma++)
								{
									for (unsigned in_beta = in_gamma+1;in_beta<=4;in_beta++)
									{  
										findIndex(in_gamma, in_beta, 1, in_ilow, in_ihigh);
										findIndex(in_gamma, in_beta, 2, in_jlow, in_jhigh);
										findIndex(in_gamma, in_beta, 3, in_klow, in_khigh);
										findIndex(in_gamma, in_beta, 4, in_llow, in_lhigh);

										for (unsigned in_i=in_ilow;in_i<=in_ihigh;in_i++)
										{
											for (unsigned in_j=in_jlow;in_j<=in_jhigh;in_j++)
											{
												for (unsigned in_k=in_klow;in_k<=in_khigh;in_k++)
												{
													for (unsigned in_l=in_llow;in_l<=in_lhigh;in_l++)
													{
														if ((in_i+in_j+in_k+in_l)==(m_P+2))
														{
															//В список будут добавляться 4 многочлена Сильвестра
															AddSilvester(in_gamma, in_beta, 1, in_i, m_InnerVectBracket);
															AddSilvester(in_gamma, in_beta, 2, in_j, m_InnerVectBracket);
															AddSilvester(in_gamma, in_beta, 3, in_k, m_InnerVectBracket);
															AddSilvester(in_gamma, in_beta, 4, in_l, m_InnerVectBracket);

															//теперь нужно добавить в список элемент свертки

															Def_nm(in_gamma,in_beta,nm2);
															n2=nm2[0];
															m2=nm2[1];

															nm2.clear();
															//----------------Перемножение скобок--------------------
															inner_bracket = m_InnerVectBracket[0];

															for (unsigned counter=1;counter<m_InnerVectBracket.size();counter++)
															{
																inner_bracket=inner_bracket*m_InnerVectBracket[counter];
															}
															//-------------------------------------------------------
															//----------------Матрица Эйлера----------------------------

															vectIn_EigFunc=FormVectEigFunc(inner_bracket,n2,m2);

															euler_br=GeneralVectorTensorVectorProduct(vectCur_EigFunc,vectIn_EigFunc,m_MatrixEps);

															std::cout<< " Rows = " << count_rows;
															std::cout<< " Columns = " << count_columns;
															std::cout << std::endl;

															std::cout<< " Euler_matr = " << Integrate(euler_br);
															std::cout << std::endl;

															*(m_EulerMatrix+m_MatrixSize*count_rows+count_columns)=Integrate(euler_br);

															//-------------------------------------------------------------------

															vectIn_EigFunc.clear();

															//--------------Метрическая матрица----------------------------------
															vectIn_RotEigFunc=RotorCalc(inner_bracket,n2,m2);

															metr_br=GeneralVectorTensorVectorProduct(vectCur_RotEigFunc,vectIn_RotEigFunc,m_MatrixMu);

															*(m_MetrMatrix+m_MatrixSize*count_rows+count_columns)=Integrate(metr_br);

															//-------------------------------------------------------------------

															vectIn_RotEigFunc.clear();

															count_columns++;
															//-------------------------------------------------------
															//----------------Удаление скобки------------------------
															m_InnerVectBracket.clear();
															inner_bracket.BracketCleanUp();
															//-------------------------------------------------------
														}
													}//l внутрениий
												}//k внутренний
											}//j внутренний
										}//i внутренний
									}//beta внутренний
								}//gamma внутренний

								nm1.clear();
								vectCur_EigFunc.clear();
								vectCur_RotEigFunc.clear();
								m_CurVectBracket.clear();
								cur_bracket.BracketCleanUp();
								count_rows++;
								count_columns=0;
							}//проверка условия (i+j+k+l==m_P+2)
						}//l
					}//k
				}//j
			}//i
		}//beta
	}//gamma
}//MatrixInit

void FiniteElementMatrix::AddToVectBracket(std::vector<double>& gains, std::vector<Power_t>& powers,
										   std::vector<Bracket>& vect_bracket)
{
	vect_bracket.push_back(Bracket(gains,powers));
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

//-----------------Определение индексов n, m если заданы gamma и beta---------------------------------
void FiniteElementMatrix::Def_nm(unsigned gamma, unsigned beta, vector<unsigned>& vect)
{
	unsigned n, m;
	if ((gamma==1)&&(beta==2))
	{
		n=3;
		m=4;
	}
	if ((gamma==1)&&(beta==3))
	{
		n=2;
		m=4;
	}
	if ((gamma==1)&&(beta==4))
	{
		n=2;
		m=3;
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
	if ((gamma==3)&&(beta==4))
	{
		n=1;
		m=2;
	}
	vect.push_back(n);
	vect.push_back(m);
}
//возвращает вектор, vect[0]=n, vect[1]=m; n<m;
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
	Power_t local_powers={0,0,0,0};
	std::vector<Power_t> powers;
	std::vector<double> gains;

	if ((gamma == numb)||(beta==numb)) //добавление несмещенного многочлена Сильвестра
	{
		if (ind==0)
		{
			if (vect_bracket.empty())
			{
				gains.push_back(1);
				powers.push_back(local_powers);
				AddToVectBracket(gains,powers,vect_bracket);
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
					gains.push_back(g);
					LocalPowersChange(local_powers,numb,1);
					powers.push_back(local_powers);
				}

				g=(-1.*s)/(s+1.0);
				if (g!=0)
				{
					gains.push_back(g);
					LocalPowersChange(local_powers,numb,0);
					powers.push_back(local_powers);
				}

				AddToVectBracket(gains,powers,vect_bracket);
				gains.clear();
				powers.clear();
			}
		}
	}
	else                    //добавление смещенного многочлена Сильвестра
	{
		if (ind==1)
		{
			if (vect_bracket.empty())
			{
				gains.push_back(1);
				powers.push_back(local_powers);
				AddToVectBracket(gains,powers,vect_bracket);
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
					gains.push_back(g);
					LocalPowersChange(local_powers,numb,1);
					powers.push_back(local_powers);
				}

				g=(-1.*s-1.0)/(s+1.0)/(s+1.0);
				if (g!=0)
				{
					gains.push_back(g);
					LocalPowersChange(local_powers,numb,0);
					powers.push_back(local_powers);
				}

				AddToVectBracket(gains,powers,vect_bracket);
				gains.clear();
				powers.clear();
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
	std::vector<Power_t> arr_local_powers=br.GetPowers();
	std::vector<double> gains=br.GetGains();
	Power_t local_powers;
	double i1, i2, i3, i4, i_s;
	for (unsigned i=0;i<gains.size();i++)
	{

		local_powers=arr_local_powers.at(i);
		pow1=local_powers.p1;
		pow2=local_powers.p2;
		pow3=local_powers.p3;
		pow4=local_powers.p4;
		i1=CalcFact(pow1);
		i2=CalcFact(pow2);
		i3=CalcFact(pow3);
		i4=CalcFact(pow4);
		i_s=CalcFact(pow1+pow2+pow3+pow4+3);
		I=I+gains[i]*i1*i2*i3*i4*6.0/i_s;
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

	std::vector<Power_t> powers, current_powers;
	std::vector<double> gains, current_gains;
	Power_t pw;
	double gn;

	current_powers=br.GetPowers();
	current_gains=br.GetGains();

	for (unsigned i=0;i<br.BracketSize();i++)
	{
		if ((current_powers[i]).p1>0)
		{
			//Добавление градиента ksi_1 с соответствующим коэффициентом,
			//если текущий элемент скобки зависит от переменной ksi_1
			for (unsigned j=0;j<m_Dim;j++)
			{
				pw.p1=(unsigned)(current_powers[i].p1-1);
				pw.p2=current_powers[i].p2;
				pw.p3=current_powers[i].p3;
				pw.p4=current_powers[i].p4;

				gn=m_GradKsi_1[j]*current_gains[i]*((double)(current_powers[i]).p1);

				powers.push_back(pw);
				gains.push_back(gn);

				local_br.SetPowers(powers);
				local_br.SetGains(gains);

				powers.clear();
				gains.clear();

				vect_nabla_phi[j]=vect_nabla_phi[j]+local_br;
			}
		}//ksi_1

		//Аналогично градиенты ksi_2, ksi_3, ksi_4
		if ((current_powers[i]).p2>0)
		{
			//Добавление градиента ksi_1 с соответствующим коэффициентом,
			//если текущий элемент скобки зависит от переменной ksi_1
			for (unsigned j=0;j<m_Dim;j++)
			{
				pw.p1=current_powers[i].p1;
				pw.p2=(unsigned)(current_powers[i].p2-1);;
				pw.p3=current_powers[i].p3;
				pw.p4=current_powers[i].p4;

				gn=m_GradKsi_2[j]*current_gains[i]*((double)(current_powers[i]).p2);

				powers.push_back(pw);
				gains.push_back(gn);

				local_br.SetPowers(powers);
				local_br.SetGains(gains);

				powers.clear();
				gains.clear();

				vect_nabla_phi[j]=vect_nabla_phi[j]+local_br;
			}
		}//ksi_2

		if ((current_powers[i]).p3>0)
		{
			//Добавление градиента ksi_1 с соответствующим коэффициентом,
			//если текущий элемент скобки зависит от переменной ksi_1
			for (unsigned j=0;j<m_Dim;j++)
			{
				pw.p1=current_powers[i].p1;
				pw.p2=current_powers[i].p2;
				pw.p3=(unsigned)(current_powers[i].p3-1);;
				pw.p4=current_powers[i].p4;

				gn=m_GradKsi_3[j]*current_gains[i]*((double)(current_powers[i]).p3);

				powers.push_back(pw);
				gains.push_back(gn);

				local_br.SetPowers(powers);
				local_br.SetGains(gains);

				powers.clear();
				gains.clear();

				vect_nabla_phi[j]=vect_nabla_phi[j]+local_br;
			}
		}//ksi_3

		if ((current_powers[i]).p4>0)
		{
			//Добавление градиента ksi_1 с соответствующим коэффициентом,
			//если текущий элемент скобки зависит от переменной ksi_1
			for (unsigned j=0;j<m_Dim;j++)
			{
				pw.p1=current_powers[i].p1;
				pw.p2=current_powers[i].p2;
				pw.p3=current_powers[i].p3;
				pw.p4=(unsigned)(current_powers[i].p4-1);

				gn=m_GradKsi_4[j]*current_gains[i]*((double)(current_powers[i]).p4);

				powers.push_back(pw);
				gains.push_back(gn);

				local_br.SetPowers(powers);
				local_br.SetGains(gains);

				powers.clear();
				gains.clear();

				vect_nabla_phi[j]=vect_nabla_phi[j]+local_br;
			}
		}//ksi_4
	}//vect_nabla_phi готов

	pw.p1=0;
	pw.p2=0;
	pw.p3=0;
	pw.p4=0;

	gains.push_back(1.0);
	local_br.SetGains(gains);

	//(ksi_n*nabla(ksi_m) - ksi_m*nabla(ksi_n)) = vect_ksi_nm
	std::vector<Bracket> vect_ksi_nm;

	LocalPowersChange(pw,n,1);
	powers.push_back(pw);
	local_br.SetPowers(powers);

	for (unsigned i=0;i<m_Dim;i++)
		vect_ksi_nm.push_back(local_br*v_1[i]);

	LocalPowersChange(pw,n,0);
	LocalPowersChange(pw,m,1);

	powers.clear();
	powers.push_back(pw);
	local_br.SetPowers(powers);

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

void FiniteElementMatrix::ShowMatrixs()
{
	std::cout << "============EulerMatrix==============" << std::endl;
	for (unsigned i=0;i<m_MatrixSize;i++)
	{
		for (unsigned j=0;j<m_MatrixSize;j++)
		{
			std::cout<< " " << *(m_EulerMatrix+i*m_MatrixSize+j);

		}
		std::cout << std::endl;
	}

	std::cout << "============MetrMatrix==============" << std::endl;
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

	Power_t pw={0,0,0,0};

	std::vector<Power_t> powers;
	std::vector<double> gains;

	LocalPowersChange(pw,n,1);

	powers.push_back(pw);
	gains.push_back(1.0);

	local_br.SetGains(gains);
	local_br.SetPowers(powers);

	Bracket res_br=local_br*(a[0]);
	result.push_back(res_br);

	res_br=local_br*(a[1]);
	result.push_back(res_br);

	res_br=local_br*(a[2]);
	result.push_back(res_br);

	powers.clear();

	LocalPowersChange(pw,n,0);
	LocalPowersChange(pw,m,1);

	powers.push_back(pw);

	local_br.SetPowers(powers);
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


//---------------Численное интегрирование трехмерного вектора, координатами которого являются скобки с помощью 
// квадратурных формул Гаусса 15-го порядка---------------------------------------------------------------------
void FiniteElementMatrix::NumIntegration(std::vector<Bracket>& vect_br, std::vector<double>& result)
{
	result.assign(m_Dim,0.0);
	std::vector<double> added_vector;
	for (unsigned i=0;i<m_QuadOrder;i++)
	{
		for (unsigned j=0;j<m_QuadOrder;j++)
		{
			for (unsigned k=0;k<m_QuadOrder;k++)
			{
				VectBracketValue(vect_br,added_vector,m_Roots[i], m_Roots[j], m_Roots[k]);
				MultNumber(added_vector, m_Weights[i]*m_Weights[j]*m_Weights[k]*pow(m_Roots[i],2)*m_Roots[j]);
				SumVector(result,added_vector);
			}
		}
	}
	MultNumber(result,0.5);
}

//-------------------------------Формирование матрицы Эйлера и метрической матрицы с использованием численного интегрирования-------
//-------------------------------------------------код сильно повторяет метод MatrixInit()------------------------------------------
void FiniteElementMatrix::NumMatrixInit()
{
	//Выделяю память под матрицу (метрическая и Эйлера)

	m_MetrMatrix = new double[m_MatrixSize*m_MatrixSize];
	m_EulerMatrix = new double[m_MatrixSize*m_MatrixSize];

	//--------------------------------------------------

	unsigned ilow, ihigh,
		jlow, jhigh,
		klow, khigh,
		llow, lhigh;

	unsigned in_ilow, in_ihigh,
		in_jlow, in_jhigh,
		in_klow, in_khigh,
		in_llow, in_lhigh;

	unsigned count_rows=0;
	unsigned count_columns=0;

	vector<Bracket> vectCur_EigFunc;
	vector<Bracket> vectIn_EigFunc;

	std::vector<Bracket> m_CurVectBracket;    //эта скобка соответсвует текущему элементу и будет храниться
	//при проходе по всем остальным элементам
	std::vector<Bracket> m_InnerVectBracket;  //эта скобка соответствует второму элементу (внутренний цикл)

	vector<Bracket> vectCur_RotEigFunc;
	vector<Bracket> vectIn_RotEigFunc;

	//--------------------Числовые векторы, соответсвующие собственным функциям и их роторам--------------------
	vector<double> numCur_EigFunc;
	vector<double> numIn_EigFunc;
	vector<double> numCur_RotEigFunc;
	vector<double> numIn_RotEigFunc;
	//---------------------------------------------------------------------------------------------------------

	Bracket cur_bracket(1);
	Bracket inner_bracket(1);

	Bracket euler_br(1);
	Bracket metr_br(1);

	unsigned n1, m1, n2, m2;
	vector<unsigned> nm1, nm2;

	double elem;

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
								AddSilvester(gamma, beta, 1, i, m_CurVectBracket);
								AddSilvester(gamma, beta, 2, j, m_CurVectBracket);
								AddSilvester(gamma, beta, 3, k, m_CurVectBracket);
								AddSilvester(gamma, beta, 4, l, m_CurVectBracket);

								//Формирую скобку cur_bracket для текущего элемента, а именно перемножаю
								//скобки массива m_CurVectBracket

								cur_bracket=m_CurVectBracket[0];

								for (unsigned counter=1; counter<m_CurVectBracket.size();counter++)
								{
									cur_bracket=cur_bracket*m_CurVectBracket[counter];
								}

								Def_nm(gamma,beta,nm1);
								n1=nm1[0];
								m1=nm1[1];

								vectCur_EigFunc=FormVectEigFunc(cur_bracket,n1,m1);
								vectCur_RotEigFunc=RotorCalc(cur_bracket,n1,m1);

								//-----------Выполняю численное интегрирование полученных собственных функций----------
								NumIntegration(vectCur_EigFunc,numCur_EigFunc);
								NumIntegration(vectCur_RotEigFunc,numCur_RotEigFunc);
								//-------------------------------------------------------------------------------------
								//теперь запускаю внутрение циклы

								for (unsigned in_gamma=1;in_gamma<=4;in_gamma++)
								{
									for (unsigned in_beta = in_gamma+1;in_beta<=4;in_beta++)
									{  
										findIndex(in_gamma, in_beta, 1, in_ilow, in_ihigh);
										findIndex(in_gamma, in_beta, 2, in_jlow, in_jhigh);
										findIndex(in_gamma, in_beta, 3, in_klow, in_khigh);
										findIndex(in_gamma, in_beta, 4, in_llow, in_lhigh);

										for (unsigned in_i=in_ilow;in_i<=in_ihigh;in_i++)
										{
											for (unsigned in_j=in_jlow;in_j<=in_jhigh;in_j++)
											{
												for (unsigned in_k=in_klow;in_k<=in_khigh;in_k++)
												{
													for (unsigned in_l=in_llow;in_l<=in_lhigh;in_l++)
													{
														if ((in_i+in_j+in_k+in_l)==(m_P+2))
														{
															//В список будут добавляться 4 многочлена Сильвестра
															AddSilvester(in_gamma, in_beta, 1, in_i, m_InnerVectBracket);
															AddSilvester(in_gamma, in_beta, 2, in_j, m_InnerVectBracket);
															AddSilvester(in_gamma, in_beta, 3, in_k, m_InnerVectBracket);
															AddSilvester(in_gamma, in_beta, 4, in_l, m_InnerVectBracket);

															//теперь нужно добавить в список элемент свертки

															Def_nm(in_gamma,in_beta,nm2);
															n2=nm2[0];
															m2=nm2[1];

															nm2.clear();
															//----------------Перемножение скобок--------------------
															inner_bracket = m_InnerVectBracket[0];

															for (unsigned counter=1;counter<m_InnerVectBracket.size();counter++)
															{
																inner_bracket=inner_bracket*m_InnerVectBracket[counter];
															}
															//-------------------------------------------------------
															//----------------Матрица Эйлера----------------------------

															vectIn_EigFunc=FormVectEigFunc(inner_bracket,n2,m2);

															NumIntegration(vectIn_EigFunc,numIn_EigFunc); //численное интегрирование 

															/*std::cout<< " Rows = " << count_rows;
															std::cout<< " Columns = " << count_columns;
															std::cout << std::endl;

															std::cout<< " Euler_matr = " << Integrate(euler_br);
															std::cout << std::endl;*/

															elem=NumericalVectorTensorVectorProduct(numCur_EigFunc,numIn_EigFunc,m_MatrixEps);
															*(m_EulerMatrix+m_MatrixSize*count_rows+count_columns)=elem;
															//-------------------------------------------------------------------

															vectIn_EigFunc.clear();

															//--------------Метрическая матрица----------------------------------
															vectIn_RotEigFunc=RotorCalc(inner_bracket,n2,m2);

															NumIntegration(vectIn_RotEigFunc,numIn_RotEigFunc); //численное интегрирование

															elem=NumericalVectorTensorVectorProduct(numCur_RotEigFunc,numIn_RotEigFunc,m_MatrixMu);
															*(m_MetrMatrix+m_MatrixSize*count_rows+count_columns)=elem;
															//-------------------------------------------------------------------

															vectIn_RotEigFunc.clear();

															count_columns++;
															//-------------------------------------------------------
															//----------------Удаление скобки------------------------
															m_InnerVectBracket.clear();
															inner_bracket.BracketCleanUp();
															//-------------------------------------------------------
														}
													}//l внутрениий
												}//k внутренний
											}//j внутренний
										}//i внутренний
									}//beta внутренний
								}//gamma внутренний

								nm1.clear();
								vectCur_EigFunc.clear();
								vectCur_RotEigFunc.clear();
								m_CurVectBracket.clear();
								cur_bracket.BracketCleanUp();
								count_rows++;
								count_columns=0;
							}//проверка условия (i+j+k+l==m_P+2)
						}//l
					}//k
				}//j
			}//i
		}//beta
	}//gamma
	double test=0;
}//NumMatrixInit
