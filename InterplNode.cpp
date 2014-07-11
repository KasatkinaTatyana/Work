
//-------по номеру интерполяционного узла функция определяет, где этот узел------------
//-------находится. Если среди индексов узла 2 нуля, то этот узел находится------------
//-------на ребре, если 1 - то на грани, если вообще нет нулей, то данный--------------
//----------------------узел принадлежит внутренней области---------------------------
void WhereIsNode(unsigned i, unsigned j, unsigned k, unsigned l,unsigned& flag)
{
	if ((i==0)&&(j==0)||(i==0)&&(k==0)||(i==0)&&(l==0)||
		(j==0)&&(k==0)||(j==0)&&(l==0)||
		(k==0)&&(l==0))
		flag = 1;
	else if ((i==0)||(j==0)||(k==0)||(l==0))
		 {
			 flag = 2;
		 }
		 else
		 {
			 flag = 3;
		 }
}

//----------------------по двум нулевым индексам определяется значение gamma, beta-------------------
// и массив index_array - номера барицентрических координат, которые не равны 0 для этого узла
void EdgeIndDefine(unsigned i, unsigned j, unsigned k, unsigned l, unsigned& gamma, unsigned& beta,
				   unsigned* index_array)
{
	if ((i==0)&&(j==0))
	{
		gamma = 1;
		beta = 2;
		index_array[0]=3;
		index_array[1]=4;
	}
	if ((i==0)&&(k==0))
	{
		gamma = 1;
		beta = 3;
		index_array[0]=2;
		index_array[1]=4;
	}
	if ((i==0)&&(l==0))
	{
		gamma = 1;
		beta = 4;
		index_array[0]=2;
		index_array[1]=3;
	}
	if ((j==0)&&(k==0))
	{
		gamma = 2;
		beta = 3;
		index_array[0]=1;
		index_array[1]=4;
	}
	if ((j==0)&&(l==0))
	{
		gamma = 2;
		beta = 4;
		index_array[0]=1;
		index_array[1]=3;
	}
	if ((k==0)&&(l==0))
	{
		gamma = 3;
		beta = 4;
		index_array[0]=1;
		index_array[1]=2;
	}
}

//-------------Определяются индексы вершин, которые образауют грань, в которой------------
//-------------лежит этот интерполяционный узел, и записываются в index_array------------
//-------------в массив non_zero_arr записываются ненулевые индексы i, j, k, l----------
void DefFaceInd(unsigned i, unsigned j, unsigned k, unsigned l, 
				unsigned* index_array, unsigned* non_zero_arr)
{
	if (i==0)
	{
		index_array[0] = 2;
		index_array[1] = 3;
		index_array[2] = 4;
		non_zero_arr[0] = j;
		non_zero_arr[1] = k;
		non_zero_arr[2] = l;
	}
	if (j==0)
	{
		index_array[0] = 1;
		index_array[1] = 3;
		index_array[2] = 4;
		non_zero_arr[0] = i;
		non_zero_arr[1] = k;
		non_zero_arr[2] = l;
	}
	if (k==0)
	{
		index_array[0] = 1;
		index_array[1] = 2;
		index_array[2] = 4;
		non_zero_arr[0] = i;
		non_zero_arr[1] = j;
		non_zero_arr[2] = l;
	}
	if (l==0)
	{
		index_array[0] = 1;
		index_array[1] = 2;
		index_array[2] = 3;
		non_zero_arr[0] = i;
		non_zero_arr[1] = j;
		non_zero_arr[2] = k;
	}
}