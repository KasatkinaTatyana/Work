#ifndef INTERPLNODE_H
#define INTERPLNODE_H

void WhereIsNode(unsigned i, unsigned j, unsigned k, unsigned l, unsigned& flag);
void EdgeIndDefine(unsigned i, unsigned j, unsigned k, unsigned l, unsigned& gamma, unsigned& beta,
				   unsigned* index_array);
void DefFaceInd(unsigned i, unsigned j, unsigned k, unsigned l, 
				unsigned* index_array, unsigned* non_zero_arr);

#endif // INTERPLNODE_H