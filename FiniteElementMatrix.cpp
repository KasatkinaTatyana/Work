#include "FiniteElementMatrix.h"
#include "bracket.h"

FiniteElementMatrix::FiniteElementMatrix(unsigned p) : m_P(p)
{

}

FiniteElementMatrix::~FiniteElementMatrix()
{

}

Power_t FiniteElementMatrix::FormPowers(unsigned n,unsigned m)
{
    Power_t powers;

    if ((n==1)||(m==2))
    {
        powers.p1=1;
        powers.p2=1;
    }
    if ((n==1)||(m==3))
    {
        powers.p1=1;
        powers.p3=1;
    }
    if ((n==1)||(m==4))
    {
        powers.p1=1;
        powers.p4=1;
    }
    if ((n==2)||(m==3))
    {
        powers.p2=1;
        powers.p3=1;
    }
    if ((n==3)||(m==4))
    {
        powers.p3=1;
        powers.p4=1;
    }

    return powers;
}

void FiniteElementMatrix::AddVTVProduct(std::vector<float> a, std::vector<float> b, std::vector<float> c,
                                        std::vector<float> d, double M[3][3], unsigned n1,
                                        unsigned m1, unsigned n2, unsigned m2)
{
       double s_a, s_b;
       std::vector<double> M_a;
       std::vector<double> M_b;
       for (unsigned i=0;i<2;i++)
       {
           s_a=0; s_b=0;
           for (unsigned j=0;j<2;j++)
           {
               s_a=s_a+a[j]*M[j][i];
               s_b=s_b+b[j]*M[j][i];
           }

           M_a[i]=s_a;
           M_b[i]=s_b;
       }

       //сначала формируем вектор коэффициентов

       double s;
       std::vector<float> gains;

       s=0;
       for (unsigned i=0;i<2;i++)
       {
           s=s+c[i]*M_a[i];
       }
       gains[0]=s;

       s=0;
       for (unsigned i=0;i<2;i++)
       {
           s=s-d[i]*M_a[i];
       }
       gains[1]=s;

       s=0;
       for (unsigned i=0;i<2;i++)
       {
           s=s-c[i]*M_b[i];
       }
       gains[2]=s;

       s=0;
       for (unsigned i=0;i<2;i++)
       {
           s=s+d[i]*M_b[i];
       }
       gains[3]=s;

       //формируем вектор степеней

       Power_t local_powers={0,0,0,0};
       std::vector<Power_t> powers;

       local_powers=FormPowers(n1,m1); //можно ли здесь так сделать?

       powers.push_back(local_powers);

       //аналогично нужно сделать еще 3 раза
       //для чисел n1,m2; m1,n2; n2,m2

       //полученную скобку нужно добавить в список

}


void FiniteElementMatrix::MatrixInit(double eps[3][3], double mu[3][3])
{
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
                                //здесь будут формироваться многочлены Сильвестра;
                                //наверно, мне надо будет создать под это отдельный метод

                                //В список будут добавляться 4 многочлена Сильвестра

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
                                                            //здесь будут формироваться многочлены Сильвестра;
                                                            //наверно, мне надо будет создать под это отдельный метод

                                                            //В список будут добавляться 4 многочлена Сильвестра

                                                            //теперь нужно добавить в список элемент свертки

                                                            //пусть процедура добавления будет addVTVProduct
                                                            //AddVTVProduct(std::vector a, std::vector b, std::vector c, std::vector d, double M[3][3], unsigned n1, unsigned m1, unsigned n2, unsigned m2);

                                                        }
                                                    }//l внутрениий
                                                }//k внутренний
                                            }//j внутренний
                                        }//i внутренний
                                    }//beta внутренний
                                }//gamma внутренний
                            }//проверка условия (i+j+k+l==p+2)

                        }//l
                    }//k
                }//j
            }//i
        }//beta
    }//gamma
}//MatrixInit
