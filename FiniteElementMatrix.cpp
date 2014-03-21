 #include "FiniteElementMatrix.h"
#include "bracket.h"

#include <cstdlib>
#include <ctime>



FiniteElementMatrix::FiniteElementMatrix(unsigned p,double simplex_peaks[m_CountPeaks][m_Dim],
                                         double Eps[m_Dim][m_Dim], double Mu[m_Dim][m_Dim]) : m_P(p)
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
    m_MatrixSize=(unsigned)((p+1)*(p+3)*(p+4)/2);

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
}

FiniteElementMatrix::~FiniteElementMatrix()
{

}

Power_t FiniteElementMatrix::DefPowers(unsigned n,unsigned m)
{
    Power_t powers={0,0,0,0};

    if (n==1)
    {
        powers.p1=powers.p1+1;
    }
    if (n==2)
    {
        powers.p2=powers.p2+1;
    }
    if (n==3)
    {
        powers.p3=powers.p3+1;
    }
    if (n==4)
    {
        powers.p4=powers.p4+1;
    }

    if (m==1)
    {
        powers.p1=powers.p1+1;
    }
    if (m==2)
    {
        powers.p2=powers.p2+1;
    }
    if (m==3)
    {
        powers.p3=powers.p3+1;
    }
    if (m==4)
    {
        powers.p4=powers.p4+1;
    }
    return powers;
}

void FiniteElementMatrix::AddVTVProduct(std::vector<double> a, std::vector<double> b, std::vector<double> c,
                                        std::vector<double> d, double M[][m_Dim], unsigned n1,
                                        unsigned m1, unsigned n2, unsigned m2,
                                        std::vector<Bracket>& vect_bracket)
{
    double s_a, s_b;
    std::vector<double> M_a;
    std::vector<double> M_b;
    for (unsigned i=0;i<m_Dim;i++)
    {
        s_a=0; s_b=0;
        for (unsigned j=0;j<m_Dim;j++)
        {
            s_a=s_a+a[j]*M[j][i];
            s_b=s_b+b[j]*M[j][i];
        }

        M_a.push_back(s_a);
        M_b.push_back(s_b);
    }

    //сначала формируем вектор коэффициентов

    double s;
    std::vector<double> gains;

    //формируем вектор степеней

    Power_t local_powers={0,0,0,0};
    std::vector<Power_t> powers;

    s=0;
    for (unsigned i=0;i<m_Dim;i++)
    {
        s=s+c[i]*M_a[i];
    }
    if (s!=0)
    {
        gains.push_back(s);
        local_powers=DefPowers(n1,n2);
        powers.push_back(local_powers);
    }

    //аналогично для чисел n1,m2; m1,n2; m2,m2
    s=0;
    for (unsigned i=0;i<m_Dim;i++)
    {
        s=s-d[i]*M_a[i];
    }
    if (s!=0)
    {
        gains.push_back(s);
        local_powers=DefPowers(n1,m2);
        powers.push_back(local_powers);
    }

    s=0;
    for (unsigned i=0;i<m_Dim;i++)
    {
        s=s-c[i]*M_b[i];
    }
    if (s!=0)
    {
        gains.push_back(s);
        local_powers=DefPowers(n2,m1);
        powers.push_back(local_powers);
    }

    s=0;
    for (unsigned i=0;i<m_Dim;i++)
    {
        s=s+d[i]*M_b[i];
    }
    if (s!=0)
    {
        gains.push_back(s);
        local_powers=DefPowers(m1,m2);
        powers.push_back(local_powers);
    }

    //полученную скобку добавляем в список
    AddToVectBracket(gains, powers,vect_bracket);
}


void FiniteElementMatrix::MatrixInit()

{
    //Выделяю память под матрицу (метрическая и Эйлера)

    m_MetrMatrix=new double*[m_MatrixSize];
    m_EulerMatrix=new double*[m_MatrixSize];

    for (unsigned i=0;i<m_MatrixSize;i++)
    {
        m_MetrMatrix[i]=new double [m_MatrixSize];
        m_EulerMatrix[i]=new double [m_MatrixSize];
    }
    //--------------------------------------------------

    unsigned ilow, ihigh,
    jlow, jhigh,
    klow, khigh,
    llow, lhigh;

    unsigned in_ilow, in_ihigh,
    in_jlow, in_jhigh,
    in_klow, in_khigh,
    in_llow, in_lhigh;

    for (unsigned gamma=1;gamma<=4;gamma++)
    {
        for (unsigned beta = gamma;beta<=4;beta++)
        {
            if ((gamma==1)||(beta==1))
            {
                ilow=0;
                ihigh=m_P;
            }
            else
            {
                ilow=1;
                ihigh=m_P+1;
            }
            if ((gamma==2)||(beta==2))
            {
                jlow=0;
                jhigh=m_P;
            }
            else
            {
                jlow=1;
                jhigh=m_P+1;
            }
            if ((gamma==3)||(beta==3))
            {
                klow=0;
                khigh=m_P;
            }
            else
            {
                klow=1;
                khigh=m_P+1;
            }
            if ((gamma==4)||(beta==4))
            {
                llow=0;
                lhigh=m_P;
            }
            else
            {
                llow=1;
                lhigh=m_P+1;
            }

            for (unsigned i=ilow;i<=ihigh;i++)
            {
                for (unsigned j=jlow;j<=jhigh;j++)
                {
                    for (unsigned k=klow;k<=khigh;k++)
                    {
                        for (unsigned l=llow;l<lhigh;l++)
                        {
                            if (i+j+k+l==m_P+2)
                            {
                                //В список добавляются 4 многочлена Сильвестра
                                //В список будут добавляться 4 многочлена Сильвестра
                                AddSilvester(gamma, beta, 1, i, m_CurVectBracket);
                                AddSilvester(gamma, beta, 2, j, m_CurVectBracket);
                                AddSilvester(gamma, beta, 3, k, m_CurVectBracket);
                                AddSilvester(gamma, beta, 4, l, m_CurVectBracket);

                                //Формирую скобку cur_bracket для текущего элемента, а именно перемножаю
                                //скобки массива m_CurVectBracket

                                Bracket cur_bracket=m_CurVectBracket[0];
                                for (unsigned counter=1; counter<m_CurVectBracket.size();counter++)
                                {
                                    cur_bracket=cur_bracket*m_CurVectBracket[counter];
                                }

                                //теперь запускаю внутрение циклы

                                for (unsigned in_gamma=1;in_gamma<=4;in_gamma++)
                                {
                                    for (unsigned in_beta = in_gamma;in_beta<=4;in_beta++)
                                    {
                                        if ((in_gamma==1)||(in_beta==1))
                                        {
                                            in_ilow=0;
                                            in_ihigh=m_P;
                                        }
                                        else
                                        {
                                            in_ilow=1;
                                            in_ihigh=m_P+1;
                                        }
                                        if ((in_gamma==2)||(in_beta==2))
                                        {
                                            in_jlow=0;
                                            in_jhigh=m_P;
                                        }
                                        else
                                        {
                                            in_jlow=1;
                                            in_jhigh=m_P+1;
                                        }
                                        if ((in_gamma==3)||(in_beta==3))
                                        {
                                            in_klow=0;
                                            in_khigh=m_P;
                                        }
                                        else
                                        {
                                            in_klow=1;
                                            in_khigh=m_P+1;
                                        }
                                        if ((in_gamma==4)||(in_beta==4))
                                        {
                                            in_llow=0;
                                            in_lhigh=m_P;
                                        }
                                        else
                                        {
                                            in_llow=1;
                                            in_lhigh=m_P+1;
                                        }

                                        for (unsigned in_i=in_ilow;in_i<=in_ihigh;in_i++)
                                        {
                                            for (unsigned in_j=in_jlow;in_j<=in_jhigh;in_j++)
                                            {
                                                for (unsigned in_k=in_klow;in_k<=in_khigh;in_k++)
                                                {
                                                    for (unsigned in_l=in_llow;in_l<in_lhigh;in_l++)
                                                    {
                                                        if (in_i+in_j+in_k+in_l==m_P+2)
                                                        {
                                                            AddToVectBracket(cur_bracket.GetGains(),
                                                                             cur_bracket.GetPowers(),
                                                                             m_VectBracket);
                                                            //В список будут добавляться 4 многочлена Сильвестра
                                                            AddSilvester(in_gamma, in_beta, 1, in_i, m_VectBracket);
                                                            AddSilvester(in_gamma, in_beta, 2, in_j, m_VectBracket);
                                                            AddSilvester(in_gamma, in_beta, 3, in_k, m_VectBracket);
                                                            AddSilvester(in_gamma, in_beta, 4, in_l, m_VectBracket);

                                                            //теперь нужно добавить в список элемент свертки
                                                            unsigned n1, m1, n2, m2;
                                                            std::vector<double> a, b, c, d;
                                                            std::vector<unsigned> nm1, nm2;
                                                            nm1=Def_nm(gamma,beta);
                                                            nm2=Def_nm(in_gamma,in_beta);
                                                            n1=nm1[0];
                                                            m1=nm1[1];
                                                            n2=nm2[0];
                                                            m2=nm2[1];
                                                            a=DefVector(n1);
                                                            b=DefVector(m1);
                                                            c=DefVector(n2);
                                                            d=DefVector(m2);
                                                            AddVTVProduct(a, b, c, d, m_MatrixEps, n1, m1, n2, m2,
                                                                          m_VectBracket);
                                                            //----------------Перемножение скобок--------------------
                                                            Bracket finite_bracket = m_VectBracket[0];
                                                            for (unsigned counter=1;counter<m_VectBracket.size();counter++)
                                                            {
                                                                finite_bracket=finite_bracket*m_VectBracket[counter];
                                                            }
                                                            //-------------------------------------------------------
                                                            //----------------Интегрирование-------------------------

                                                            //-------------------------------------------------------
                                                            //----------------Удаление скобки------------------------
                                                            m_VectBracket.clear();
                                                            finite_bracket.BracketCleanUp();
                                                            //-------------------------------------------------------
                                                        }
                                                    }//l внутрениий
                                                }//k внутренний
                                            }//j внутренний
                                        }//i внутренний
                                    }//beta внутренний
                                }//gamma внутренний

                                m_CurVectBracket.clear();
                                cur_bracket.BracketCleanUp();

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

void FiniteElementMatrix::ShowVectBracket(std::vector<Bracket>& vect_bracket)
{
    for (unsigned i=0;i<m_VectBracket.size();i++)
    {
        ((Bracket)(vect_bracket.at(i))).ShowElements();
    }
}
//-----------------------Определение вектора----------------------------------------------------------
std::vector<double> FiniteElementMatrix::DefVector(unsigned ind)
{
    std::vector<double> result;
    if (ind==1)
        result=m_GradKsi_1;
    if (ind==2)
        result=m_GradKsi_2;
    if (ind==3)
        result=m_GradKsi_3;
    if (ind==4)
        result=m_GradKsi_4;
    return result;
}
//----------------------------------------------------------------------------------------------------

//-----------------Определение индексов n, m если заданы gamma и beta---------------------------------
std::vector<unsigned> FiniteElementMatrix::Def_nm(unsigned gamma, unsigned beta)
{
    unsigned n, m;
    std::vector<unsigned> vect;
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
    if ((gamma==3)&&(beta==4))
    {
        n=1;
        m=2;
    }
    vect.push_back(n);
    vect.push_back(m);
    return vect;
}
//возвращает вектор, vect[0]=n, vect[1]=m; n<m;
//----------------------------------------------------------------------------------------------------
//----------------Векторное произведение двух векторов------------------------------------------------
std::vector<double> FiniteElementMatrix::VectProduct(std::vector<double> a, std::vector<double> b)
{
    std::vector<double> result;
    result.push_back(a[1]*b[2]-a[2]*b[1]);
    result.push_back(a[2]*b[0]-a[0]*b[2]);
    result.push_back(a[0]*b[1]-a[1]*b[0]);
    return result;
}
//--------------------------------------------------------------------------------------------------------
//---------------Векторное произведение двух векторов, элементами кот. являются скобки--------------------
std::vector<Bracket> FiniteElementMatrix::VectBracketProduct(std::vector<Bracket> a, std::vector<Bracket> b)
{
    std::vector<Bracket> result;

    Bracket b1=a[1]*b[2];
    Bracket b2=a[2]*b[1];

    result.push_back(b1-b2);

    b1.BracketCleanUp();
    b2.BracketCleanUp();

    b1=a[2]*b[0];
    b2=a[0]*b[2];

    result.push_back(b1-b2);

    b1.BracketCleanUp();
    b2.BracketCleanUp();

    b1=a[0]*b[1];
    b2=a[1]*b[0];

    result.push_back(b1-b2);

    b1.BracketCleanUp();
    b2.BracketCleanUp();

    return result;
}
//----------------------------------------------------------------------------------------------------
//-------------Скалярное произведение двух векторов---------------------------------------------------
double FiniteElementMatrix::ScalarProduct(std::vector<double> a, std::vector<double> b)
{
    double result;
    result=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    return result;
}
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
                g=(m_P+2.0)/(s+1.0);
                if (g!=0)
                {
                    gains.push_back(g);
                    LocalPowersChange(local_powers,numb,1);
                    powers.push_back(local_powers);
                }

                g=(-s+0.0)/(s+1.0);
                if (g!=0)
                {
                    gains.push_back(g);
                    LocalPowersChange(local_powers,numb,0);
                    powers.push_back(local_powers);
                }

                AddToVectBracket(gains,powers,vect_bracket);
                gains.clear();
            }
        }
    }
    else                    //добавление смещенного многочлена Сильвестра
    {
        if (ind==1)
        {
            if (m_VectBracket.empty())
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
                g=(m_P+2.0)/(s+1.0);
                if (g!=0)
                {
                    gains.push_back(g);
                    LocalPowersChange(local_powers,numb,1);
                    powers.push_back(local_powers);
                }
                g=(-s-1.0)/(s+1.0);
                if (g!=0)
                {
                    gains.push_back(g);
                    LocalPowersChange(local_powers,numb,0);
                    powers.push_back(local_powers);
                }

                AddToVectBracket(gains,powers,vect_bracket);
                gains.clear();
            }
        }
    }
}
//---------------------------------------------------------------------------------------------
//----------В переменной local_powers присвоить полю с номером ind значение value--------------
void FiniteElementMatrix::LocalPowersChange(Power_t& local_powers, unsigned ind, unsigned value)
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
double FiniteElementMatrix::Integrate (Bracket_t& br)
{
    double I=0;
    double pow1, pow2, pow3, pow4;
    std::vector<Power_t> arr_local_powers=br.GetPowers();
    Power_t local_powers;
    double i1, i2, i3, i4, i_s;
    for (unsigned i=0;i<(br.GetGains()).size();i++)
    {

        local_powers=arr_local_powers.at(i);
        pow1=local_powers.p1;
        pow2=local_powers.p2;
        pow3=local_powers.p3;
        pow4=local_powers.p4;
        i1=CalcFact(unsigned(pow1));
        i2=CalcFact(unsigned(pow2));
        i3=CalcFact(unsigned(pow3));
        i4=CalcFact(unsigned(pow4));
        i_s=CalcFact(unsigned(pow1+pow2+pow3+pow4+3));
        I=I+(br.GetGains())[i]*i1*i2*i3*i4/i_s;
    }
    return I;
}
//--------------------------Произведение вектора тензора и вектора в общем виде-------------------------
Bracket FiniteElementMatrix::GeneralVectorTensorVectorProduct(std::vector<Bracket>& vect1,
                                                              std::vector<Bracket>& vect2, double M[][m_Dim])
{
    //создаю скобку, которая ничего не содержит
    std::vector<double> zero_gains;
    zero_gains.push_back(0.0);
    std::vector<Power_t> zero_powers;
    Power_t p={0, 0, 0, 0};
    zero_powers.push_back(p);

    Bracket result(1);
    result.SetGains(zero_gains);
    result.SetPowers(zero_powers);

    Bracket br_product=vect1[0]*(M[0][0]);

    Bracket br_sum(1);

    for (unsigned i=0;i<m_Dim;i++)
    {
        br_sum.SetGains(zero_gains);
        br_sum.SetPowers(zero_powers);

        for (unsigned k=0;k<m_Dim;k++)
        {
            if ((i!=0)||(k!=0))
            {
                br_product=vect1[k]*(M[k][i]);
            }
            br_sum=br_sum+br_product;
        }
        br_sum=br_sum*vect2[i];

        result=result+br_sum;

        br_sum.BracketCleanUp();
    }
    //AddToVectBracket(result.GetGains(),result.GetPowers(),m_VectBracket);
    return result;
}

//----------Вычисление ротора от вектора, который задан в виде: скобка * (ksi_n*nabla(ksi_m) - ksi_m*nabla(ksi_n))
// n < m
std::vector<Bracket> FiniteElementMatrix::RotorCalc(Bracket br, unsigned n, unsigned m)
{
    std::vector<Bracket> result;
    //Скалярная функция * градиент вектора
    std::vector<double> vect=ProductGradKsi(n,m);
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
    std::vector<double> v_1=DefVector(m);
    std::vector<double> v_2=DefVector(n);

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

    std::vector<Bracket> res_1=VectBracketProduct(vect_nabla_phi,vect_ksi_nm);

    for (unsigned i=0;i<m_Dim;i++)
        result[i]=result[i]+res_1[i];

    /*for (unsigned i=0;i<m_Dim;i++)
    {
        ((Bracket)(result[i])).ShowElements();
    }*/

    return result;
}
//------------------------------Вычисление nabla(ksi_1)*nabla(ksi_2)----------------------------------------------
std::vector<double> FiniteElementMatrix::ProductGradKsi(unsigned ind1, unsigned ind2)
{
    std::vector<double> result;

    if((ind1==1)&&(ind2==2))
        result = m_UnVect_3;

    if((ind1==2)&&(ind2==1))
    {
        for (unsigned i=0;i<m_Dim;i++)
        {
            result.push_back((-1.0)*m_UnVect_3[i]);
        }
    }

    if ((ind1==1)&&(ind2==3))
    {
        for (unsigned i=0;i<m_Dim;i++)
        {
            result.push_back((-1.0)*m_UnVect_2[i]);
        }
    }

    if((ind1==3)&&(ind2==1))
        result = m_UnVect_2;

    if ((ind1==2)&&(ind2==3))
        result = m_UnVect_1;

    if ((ind1==3)&&(ind2==2))
    {
        for (unsigned i=0;i<m_Dim;i++)
        {
            result.push_back((-1.0)*m_UnVect_1[i]);
        }
    }

    if ((ind1==1)&&(ind2==4))
    {
        for (unsigned i=0;i<m_Dim;i++)
        {
            result.push_back((-1.0)*m_UnVect_3[i]+m_UnVect_2[i]);
        }
    }

    if ((ind1==4)&&(ind2==1))
    {
        for (unsigned i=0;i<m_Dim;i++)
        {
            result.push_back(m_UnVect_3[i]+(-1.0)*m_UnVect_2[i]);
        }
    }

    if ((ind1==2)&&(ind2==4))
    {
        for (unsigned i=0;i<m_Dim;i++)
        {
            result.push_back((-1.0)*m_UnVect_1[i]+m_UnVect_3[i]);
        }
    }

    if ((ind1==4)&&(ind2==2))
    {
        for (unsigned i=0;i<m_Dim;i++)
        {
            result.push_back(m_UnVect_1[i]+(-1.0)*m_UnVect_3[i]);
        }
    }

    if ((ind1==3)&&(ind2==4))
    {
        for (unsigned i=0;i<m_Dim;i++)
        {
            result.push_back(m_UnVect_1[i]+(-1.0)*m_UnVect_2[i]);
        }
    }

    if ((ind1==4)&&(ind2==3))
    {
        for (unsigned i=0;i<m_Dim;i++)
        {
            result.push_back((-1.0)*m_UnVect_1[i]+m_UnVect_2[i]);
        }
    }

    for (unsigned i=0;i<m_Dim;i++)
    {
        result[i]=result[i]/m_Icob;
    }

    return result;
}
