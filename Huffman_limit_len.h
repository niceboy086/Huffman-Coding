
// Huffman.h : 头文件
//


#pragma once

#include <deque>
#include <algorithm>
#include <cmath>
#include <map>
#include <vector>
using namespace std;


/*哈夫曼编码*/
class CHuffmanCode
{
public:
	CHuffmanCode()
		:m_CodeLen(0), m_pCode(nullptr)
	{
	}
	CHuffmanCode(int codeLen)
		:m_CodeLen(codeLen) 
	{
		m_pCode = new char[m_CodeLen];
		for(int i=0; i<m_CodeLen; i++) 
			m_pCode[i] = 0;
	}
	CHuffmanCode(int codeLen, CHuffmanCode& initVal)
		:m_CodeLen(codeLen) 
	{
		m_pCode = new char[m_CodeLen];
		for(int i=0; i<m_CodeLen; i++) 
			m_pCode[i] = 0;

		int start;
		char * ptr = m_pCode + m_CodeLen - 1;

		if(m_CodeLen >= initVal.m_CodeLen)
		{
			start = 0;
		}
		else
		{
			start = initVal.m_CodeLen - m_CodeLen;
		}

		for(int i=initVal.m_CodeLen-1; i>=start; i--)
		{
			*ptr = initVal.m_pCode[i];
			ptr--;
		}
	}

	~CHuffmanCode() 
	{
		if(m_pCode != nullptr) 
		{
			delete[] m_pCode;
			m_pCode = nullptr;
		}
		m_CodeLen = 0;
	}

	CHuffmanCode& operator+=(char bit) 
	{
		for(int i=m_CodeLen-1; i>=0; i--)
		{
			m_pCode[i] += bit;

			if(m_pCode[i] < 2)
			{
				break;
			}

			m_pCode[i] = 0;
		}

        return *this;
    }

	CHuffmanCode& operator<<=(int k) 
	{
		if(k >= m_CodeLen)
		{
			for(int i=m_CodeLen-1; i>=0; i--)
			{
				m_pCode[i] = 0;
			}

		}
		else
		{
			for(int i=0; i<m_CodeLen-k; i++)
			{
				m_pCode[i] = m_pCode[i+k];
			}
			for(int i=m_CodeLen-k; i<m_CodeLen; i++)
			{
				m_pCode[i] = 0;
			}
		}

        return *this;
    }

	void GetCode(char ** ppCode, int * pCodeLen)
	{
		if(m_pCode != nullptr)
		{
			*ppCode = new char[m_CodeLen];
			memcpy(*ppCode, m_pCode, m_CodeLen);
			*pCodeLen = m_CodeLen;
		}
		else
		{
			*ppCode = nullptr;
			*pCodeLen = 0;
		}
	}

	void SimpReset()
	{
		m_pCode = nullptr;
		m_CodeLen = 0;
	}

	char * GetCodePtr()
	{
		return m_pCode;
	}
	int GetCodeLen()
	{
		return m_CodeLen;
	}

	void PaintSkin(CHuffmanCode& hfmcode)
	{
		if(this != &hfmcode)
		{
			if(m_pCode != nullptr)
			{
				delete[] m_pCode;
			}
			m_pCode = hfmcode.GetCodePtr();
			m_CodeLen = hfmcode.GetCodeLen();
			hfmcode.SimpReset();
		}
	}


private:
	char * m_pCode;
	int m_CodeLen;
};

/*哈夫曼树的节点定义*/
template <typename T>
struct HuffmanNode
{
    HuffmanNode(T k,HuffmanNode<T>*l=nullptr,HuffmanNode<T>* r=nullptr,HuffmanNode<T>* p=nullptr)
        :key(k), lchild(l), rchild(r), parent(p), code(-1), idx(-1){}
    ~HuffmanNode(){};
	int idx;					// 创建时节点索引, nodes和create()参数a[]中的索引, 解码时输出这个索引
	char code;					// 节点哈夫曼编码
    T key;                      //节点的权值
    HuffmanNode<T>* parent;        //父节点
    HuffmanNode<T>* lchild;        //节点左孩
    HuffmanNode<T>* rchild;        //节点右孩
};

typedef char*	_CodePtr;

template <typename T>
class CHuffman
{
public:
    void preOrder();              //前序遍历哈夫曼树
    void inOrder();                  //中序遍历哈夫曼树
    void postOrder();              //后序遍历哈夫曼树
 
	// iLimit: 0 - 不限长, 大于0限长
	bool CanonicCreat(T w[],int size, int iLimit = 0);
    void creat(T a[], int size);		//创建哈夫曼树
    void recreat();						//根据范式编码重建哈夫曼树
    void destroy();						//销毁哈夫曼树
    void print();						//打印哈夫曼树
    void printCode();						//打印哈夫曼树
	bool getCode(int iIndex, char ** ppCode, int * plen);	// 指定索引的哈夫编码, 返回编码长度
	bool getCode(int * plen);	// 返回编码长度
	void Reset(){ destroy(); }
	HuffmanNode<T>* GetRoot(){return root;}
	void ClearCodePtr();

    CHuffman();
    virtual ~CHuffman(){ destroy(); ClearCodePtr(); TRACE("called destructor of class CHuffman!\r\n"); };
 
private:
    void preOrder(HuffmanNode<T>* pnode);
    void inOrder(HuffmanNode<T>* pnode);
    void postOrder(HuffmanNode<T>*pnode);
    void print(HuffmanNode<T>*pnode);
	// 返回编码长度
	bool getCodeLen();
	// 限定编码长度
	bool CodeLenLimit(int iLimitLen);
	void CanonicCodeByLens(T w[],int sz);

protected:
	vector<int>	m_vecCodeLens;		// 编码长度
	_CodePtr * m_pCodePtr;
 
private:
    HuffmanNode<T>* root;			//哈夫曼树根节点
    deque<HuffmanNode<T>*> nodes;	//叶子节点
    deque<HuffmanNode<T>*> forest;	//森林
};

template<typename T>
void CHuffman<T>::ClearCodePtr()
{
	if(m_pCodePtr != nullptr)
	{
		int size = m_vecCodeLens.size();
		for(int i=0; i<size; i++)
		{
			if(m_pCodePtr[i] != nullptr)
			{
				delete[] m_pCodePtr[i];
				m_pCodePtr[i] = nullptr;
			}
		}

		delete[] m_pCodePtr;
		m_pCodePtr = nullptr;
	}
	m_vecCodeLens.clear();
}

template<typename T>
void CHuffman<T>::printCode()						//打印哈夫曼树
{
	int size = m_vecCodeLens.size();

	for(int i=0; i<size; i++)
	{
		int len = m_vecCodeLens[i];
		_CodePtr codePtr = m_pCodePtr[i];

		for(int j=0; j<len; j++)
		{
			TRACE("%d", codePtr[j]);
		}
		TRACE(" -- %d\r\n", i);
	}
}

/*创建哈夫曼树*/
template<typename T>
CHuffman<T>::CHuffman()
{
	root = nullptr;
	m_pCodePtr = nullptr;
}

// 返回编码长度
template<typename T>
bool CHuffman<T>::getCodeLen()
{
	int size = nodes.size();

	if(size == 1)
	{
		m_vecCodeLens.push_back(1);
	}
	else
	{
		for(int j=0; j<size; j++)
		{
			int codeLen = 0;
			HuffmanNode<T>* node = nodes[j];
			while(node->parent != nullptr)
			{
				codeLen += 1;
				node = node->parent;
			}

			m_vecCodeLens.push_back(codeLen);
		}
	}

	return true;
}

template<typename T>
bool CHuffman<T>::CodeLenLimit(int iLimitLen)
{
	if(m_vecCodeLens.empty())
		return false;

	sort(m_vecCodeLens.begin(), m_vecCodeLens.end());

	if(m_vecCodeLens.size() <= 3)
	{
		if(iLimitLen < m_vecCodeLens.back())
		{
			// 该情况不能减码长
			return false;
		}
		else
		{
			// 符合码长要求, 返回true
			return true;
		}
	}

	bool brst = true;

	while(m_vecCodeLens.back() > iLimitLen)
	{
		// 取出最后两个
		int iLast1 = m_vecCodeLens.back();
		m_vecCodeLens.pop_back();
		int iLast2 = m_vecCodeLens.back();
		m_vecCodeLens.pop_back();

		if(iLast1 != iLast2)
		{
			// 正常这两个应该相等
			brst = false;
			break;
			// return false;
		}

		bool bfind = false;
		vector<int>::iterator iter = m_vecCodeLens.end();
		while(iter != m_vecCodeLens.begin())
		{
			iter--;
			if(*iter <= iLimitLen - 1)
			{
				bfind = true;
				break;
			}
		}

		if(bfind)
		{
			int icodelen = *iter;
			m_vecCodeLens.erase(iter);
			m_vecCodeLens.push_back(icodelen + 1);
			m_vecCodeLens.push_back(icodelen + 1);
			m_vecCodeLens.push_back(iLast1 - 1);
			sort(m_vecCodeLens.begin(), m_vecCodeLens.end());
		}
		else
		{
			brst = false;
			break;
		}
	}

	return brst;
}

template<typename T>
void CHuffman<T>::CanonicCodeByLens(T w[],int sz)
{
	// 编码长度与索引对应列表, 并按编码长度排序
	deque<pair<int, int>>	deqCodeLen;

	int size = m_vecCodeLens.size();

	for(int i=0; i<size; i++)
	{
		deqCodeLen.push_back(pair<int,int>(m_vecCodeLens[i], i));
	}
	// 排序
	sort(deqCodeLen.begin(), deqCodeLen.end(), [](pair<int,int> & it1, pair<int,int> & it2){return it1.first < it2.first;});

	m_pCodePtr = new _CodePtr[size];

	int codeLen = deqCodeLen[0].first;
	int idx = deqCodeLen[0].second;

	CHuffmanCode hfmCode(codeLen);
	int icodelen = 0;
	hfmCode.GetCode(m_pCodePtr + idx, &icodelen);

	for(int ii=0; ii<icodelen; ii++)
	{
		TRACE("%d",m_pCodePtr[idx][ii]);
	}
	TRACE(" -- idx:%d, cnt:%d, code len:%d, id:%d - huffman code\r\n", idx, w[idx], codeLen, 0);

	for(int i=1; i<size; i++)
	{
		int curcodeLen = deqCodeLen[i].first;
		int curidx = deqCodeLen[i].second;

		CHuffmanCode curhfmCode(curcodeLen, hfmCode);
		curhfmCode += 1;

		if(curcodeLen > codeLen)
		{
			curhfmCode <<= (curcodeLen - codeLen);
		}

		int icodelen = 0;
		curhfmCode.GetCode(m_pCodePtr + curidx, &icodelen);
		codeLen = curcodeLen;
		hfmCode.PaintSkin(curhfmCode);

		for(int ii=0; ii<icodelen; ii++)
		{
			TRACE("%d",m_pCodePtr[curidx][ii]);
		}
		TRACE(" -- idx:%d, cnt:%d, code len:%d, id:%d - huffman code\r\n", curidx, w[curidx], curcodeLen, i);
	}
}


template<typename T>
bool CHuffman<T>::CanonicCreat(T w[],int size, int iLimit)
{
	creat(w,size);

	// 返回编码长度
	getCodeLen();
	if(iLimit > 0)
	{
		// 限长时要先对权重和原字符进行排序, 权重一级, 字符二级, 两级都是由大到小
		bool bok = CodeLenLimit(iLimit);
		if(!bok)
			return false;
	}

	//范式编码
	CanonicCodeByLens(w,size);
	// 范式编码重建哈夫曼树
	destroy();
	recreat();
	return true;
}

//根据编码重建哈夫曼树
template<typename T>
void CHuffman<T>::recreat()
{
	int size = m_vecCodeLens.size();

	if(size == 1)
	{
		root = new HuffmanNode<T>(0,nullptr,nullptr,nullptr); 
		root->idx = 0;
		root->code = 0;
	}
	else
	{
		root = new HuffmanNode<T>(0,nullptr,nullptr,nullptr); 

		for(int i=0; i<size; i++)
		{
			_CodePtr codePtr = m_pCodePtr[i];
			int codeLen = m_vecCodeLens[i];

			HuffmanNode<T>* node = root;

			for(int j=0; j<codeLen; j++)
			{
				if(codePtr[j] == 0)
				{
					if(node->lchild == nullptr)
					{
						HuffmanNode<T> * newnode = new HuffmanNode<T>(0,nullptr,nullptr,node); 
						newnode->code = 0;
						node->lchild = newnode;
						if(j == (codeLen - 1))
						{
							newnode->idx = i;
							nodes.push_back(newnode);
						}
					}

					node = node->lchild;
				}

				if(codePtr[j] == 1)
				{
					if(node->rchild == nullptr)
					{
						HuffmanNode<T> * newnode = new HuffmanNode<T>(0,nullptr,nullptr,node); 

						newnode->code = 1;
						node->rchild = newnode;

						if(j == (codeLen - 1))
						{
							newnode->idx = i;
							nodes.push_back(newnode);
						}
					}

					node = node->rchild;
				}
			}

		}
	}
}

/*创建哈夫曼树*/
template<typename T>
void CHuffman<T>::creat(T a[],int size)
{
    for (int i = 0; i < size; i++) //每个节点都作为一个森林
    {
        //为初始序列的元素构建节点。每个节点作为一棵树加入森林中。
        HuffmanNode<T>* ptr = new HuffmanNode<T>(a[i],nullptr,nullptr,nullptr);  
		ptr->idx = i;
        nodes.push_back(ptr);
        forest.push_back(ptr);
		
    }

	if(size == 1)
	{
		forest[0]->code = 0;
	}

    for (int i = 0; i < size - 1; i++)
    {
        //排序，以选出根节点权值最小两棵树
		//for(int i=0; i<forest.size(); i++)
		//	TRACE("%.1f(l%d, r%d), ", (float)forest[i]->key, (short)forest[i]->lchild, (short)forest[i]->rchild);
		//TRACE("org\r\n");
        sort(forest.begin(), forest.end(), [](HuffmanNode<T>* a, HuffmanNode<T>*b){return a->key < b->key; });    
		//for(int i=0; i<forest.size(); i++)
		//	TRACE("%.1f(l%d, r%d), ", (float)forest[i]->key, (short)forest[i]->lchild, (short)forest[i]->rchild);
		//TRACE("sort\r\n");
        HuffmanNode<T>*node = new HuffmanNode<T>(forest[0]->key + forest[1]->key, forest[0], forest[1]); //构建新节点
        forest[0]->parent = node;
        forest[1]->parent = node;
		forest[0]->code = 0;
		forest[1]->code = 1;
		forest.push_back(node);  //新节点加入森林中
        forest.pop_front(); //删除两棵权值最小的树
        forest.pop_front();
    }
    root = forest.front();
    forest.clear();
}

/*打印哈夫曼树*/
template<typename T>
void CHuffman<T>::print()
{
	print(root);
}

/*打印哈夫曼树*/
template<typename T>
void CHuffman<T>::print(HuffmanNode<T>* pnode)
{
    if (pnode != nullptr)
    {
		TRACE("当前结点: %.2f. ", (float)pnode->key);
        if (pnode->lchild != nullptr)
            TRACE("它的左孩子节点为： %.2f.", (float)pnode->lchild->key);
        else 
			TRACE("它没有左孩子.");
        if (pnode->rchild != nullptr)
            TRACE("它的右孩子节点为： %.2f.", (float)pnode->rchild->key);
        else 
			TRACE("它没有右孩子.");
		TRACE("\r\n");

        print(pnode->lchild);
        print(pnode->rchild);
    }
}

//销毁哈夫曼树
template<typename T>
void CHuffman<T>::destroy()
{
	//结点集
    deque<HuffmanNode<T>*> nodes; 
	HuffmanNode<T>*pnode = root;

    while(pnode != nullptr)
    {
		if(pnode->lchild != nullptr)
		{
			nodes.push_back(pnode->lchild);
		}
		if(pnode->rchild != nullptr)
		{
			nodes.push_back(pnode->rchild);
		}

		delete pnode;
		pnode = nullptr;

		if(nodes.empty())
		{
			break;
		}

		pnode = nodes[0];
		nodes.pop_front();
	}

	this->nodes.clear();
	root = nullptr;
}

// ppCode - 指定索引的哈夫编码, len - 返回编码长度
template<typename T>
bool CHuffman<T>::getCode(int iIndex, char ** ppCode, int * plen)
{
	if(root == nullptr)
	{	
		return false;
	}

	int size = nodes.size();
	if(iIndex >= size)
	{
		return false;
	}

	if(size == 1)
	{
		*ppCode = new char[1];
		(*ppCode)[0] = nodes[0]->code;
		*plen = 1;
	}
	else
	{
		int codeLen = 0;
		HuffmanNode<T>* node = nodes[iIndex];
		while(node->parent != nullptr)
		{
			codeLen += 1;
			node = node->parent;
		}

		char * pcode = new char[codeLen];
		node = nodes[iIndex];
		for(int i=codeLen-1; i>=0; i--)
		{
			pcode[i] = node->code;
			node = node->parent;
		}

		*ppCode = pcode;
		*plen = codeLen;
	}

	return true;
}


template<typename T>
class CElemStat
{
public:
	typedef pair <T, int> Elem_Pair;

	int Stat(T * pText, int size);
	int GetElemNum();
	int GetStat(T * pElems, int * pCnts, int size);
	void Clear();

	CElemStat(){}
    virtual ~CElemStat(){TRACE("called destructor of class CElemStat!\r\n");}
 
private:
	map<T, int>		m_mapStat;
};

template<typename T>
void CElemStat<T>::Clear()
{
	if(!m_mapStat.empty())
		m_mapStat.clear();
}

template<typename T>
int CElemStat<T>::Stat(T * pText, int size)
{
	Clear();

	for(int i=0; i<size; i++)
	{
		map<T, int>::iterator iter = m_mapStat.find(pText[i]);
		if(iter == m_mapStat.end())
		{
			m_mapStat.insert(Elem_Pair(pText[i], 1));
		}
		else
		{
			iter->second++;
		}
	}

	return m_mapStat.size();
}

template<typename T>
int CElemStat<T>::GetStat(T * pElems, int * pCnts, int size)
{
	int rsize = m_mapStat.size();
	if(rsize <= size)
	{
		T * ptrElem = pElems;
		int * ptrCnt = pCnts;

		map<T, int>::iterator iter = m_mapStat.begin();
		while(iter != m_mapStat.end())
		{
			*ptrElem = iter->first;
			*ptrCnt = iter->second;
			ptrElem++;
			ptrCnt++;
			iter++;
		}

		return rsize;
	}
	else
	{
		return 0;
	}
}

template<typename T>
int CElemStat<T>::GetElemNum()
{
	int rsize = m_mapStat.size();
	return rsize;
}

template<typename _EL, typename _WT>
class CHuffmanCodec: 
	public CElemStat<_EL>, public CHuffman<_WT>
{
public:
	typedef CElemStat<_EL> _ElemStat;
	typedef CHuffman<_WT> _Huffman;
	typedef char*	_CodePtr;

	CHuffmanCodec():_ElemStat(), _Huffman()
	{
		m_pElems = nullptr;
		m_pWeights = nullptr;
		m_pCodePtr = nullptr;
		m_pCodeLen = nullptr;
		m_iElemNum = 0;
	}
	virtual ~CHuffmanCodec(){Reset(); TRACE("called destructor of class CHuffmanCodec!\r\n");}

	int Encode(_EL * pText, int iTextLen, char ** ppOutput, int * pOutputLen);
	int Decode(_EL * pText, int iTextLen, char ** ppOutput, int * pOutputLen);

public:
	void Reset();
	static bool cmp(pair<_EL, _WT> & it1, pair<_EL, _WT> & it2);

private:
	_EL * m_pElems;
	_WT * m_pWeights;
	_CodePtr *	m_pCodePtr;
	int * m_pCodeLen;
	map<_EL, int>	m_mapElemIdx;
	int	  m_iElemNum;
	vector<char>	m_vecEnText;
	vector<char>	m_vecDeText;
};


template<typename _EL, typename _WT>
void CHuffmanCodec<_EL, _WT>::Reset()
{
	if(m_pElems != nullptr)
	{
		delete[] m_pElems;
		m_pElems = nullptr;
	}
	if(m_pWeights != nullptr)
	{
		delete[] m_pWeights;
		m_pWeights = nullptr;
	}
	if(m_pCodePtr != nullptr)
	{
		for(int i=0; i<m_iElemNum; i++)
		{
			delete[] m_pCodePtr[i];
			m_pCodePtr[i] = nullptr;
		}

		delete[] m_pCodePtr;
		m_pCodePtr = nullptr;
	}

	if(m_pCodeLen != nullptr)
	{
		delete[] m_pCodeLen;
		m_pCodeLen = nullptr;
	}

	m_mapElemIdx.clear();
	m_vecEnText.clear();

	_ElemStat::Clear();
	_Huffman::Reset();
}

template<typename _EL, typename _WT>
bool CHuffmanCodec<_EL, _WT>::cmp(pair<_EL, _WT> & it1, pair<_EL, _WT> & it2)
{
	if(it1.second > it2.second)
	{return true;}
	else if(it1.second == it2.second)
	{return it1.first > it2.first;}
	else
	{return false;}
}

template<typename _EL, typename _WT>
int CHuffmanCodec<_EL, _WT>::Encode(_EL * pText, int iTextLen, char ** ppOutput, int * pOutputLen)
{
	Reset();

	int elemnum = Stat(pText, iTextLen);

	m_pElems = new _EL[elemnum];
	m_pWeights = new _WT[elemnum];
	m_iElemNum = GetStat(m_pElems, m_pWeights, elemnum);

	// int iLimit = 0;
	int iLimit = 3;
	if(iLimit > 0)
	{
		// 限长时要先对权重和原字符进行排序, 权重一级, 字符二级, 两级都是由大到小
		deque<pair<_EL, _WT>>	deq;
		for(int i=0; i<m_iElemNum; i++)
		{
			deq.push_back(pair<_EL, _WT>(m_pElems[i], m_pWeights[i]));
		}
		sort(deq.begin(), deq.end(), cmp);
		for(int i=0; i<m_iElemNum; i++)
		{
			pair<_EL, _WT> & it = deq[i];
			m_pElems[i] = it.first;
			m_pWeights[i] = it.second;
		}
	}

	TRACE("Elem: ");
	for(int i=0; i<elemnum; i++)
	{
		TRACE("%c, ", m_pElems[i]);
	}
	TRACE("\r\n");

	TRACE("cnt: ");
	for(int i=0; i<elemnum; i++)
	{
		TRACE("%d, ", m_pWeights[i]);
	}
	TRACE("\r\n");

	TRACE("idx: ");
	for(int i=0; i<elemnum; i++)
	{
		TRACE("%d, ", i);
	}
	TRACE("\r\n");

	bool bok = CanonicCreat(m_pWeights, m_iElemNum, iLimit);
	if(!bok)
	{
		return -1;
	}

	for(int i=0; i<elemnum; i++)
	{
		m_mapElemIdx[m_pElems[i]] = i;
	}

	m_pCodePtr = new _CodePtr[elemnum];
	m_pCodeLen = new int[elemnum];

	for(int i=0; i<elemnum; i++)
	{
		m_pCodeLen[i] = 0;
		m_pCodePtr[i] = nullptr;
		bool bret = getCode(i, m_pCodePtr + i, m_pCodeLen + i);
		for(int j=0; j<m_pCodeLen[i]; j++)
		{
			//TRACE("%d", m_pCodePtr[i][j]);
		}
		//TRACE(" -- char: %c (%d) - cnt: %d 哈夫曼编码.\r\n", m_pElems[i], m_pElems[i], m_pWeights[i]);
	}
	
	for(int i=0; i<iTextLen; i++)
	{
		int idx = m_mapElemIdx[pText[i]];
		_CodePtr codeptr = CHuffman<_WT>::m_pCodePtr[idx];

		for(int j=0; j<CHuffman<_WT>::m_vecCodeLens[idx]; j++)
		{
			m_vecEnText.push_back(codeptr[j]);
		}
	}

	int iEnTextLen = m_vecEnText.size();
	char * pEnText = new char[iEnTextLen];

	for(int i=0; i<iEnTextLen; i++)
	{
		pEnText[i] = m_vecEnText[i];
	}

	*ppOutput = pEnText;
	*pOutputLen = iEnTextLen;
	
	return iEnTextLen;
}

template<typename _EL, typename _WT>
int CHuffmanCodec<_EL, _WT>::Decode(_EL * pText, int iTextLen, char ** ppOutput, int * pOutputLen)
{
	HuffmanNode<_WT>*pnode = GetRoot();

	if((pnode->lchild == nullptr) && (pnode->rchild == nullptr))
	{
		for(int i=0; i<iTextLen; i++)
		{
			m_vecDeText.push_back(m_pElems[0]);
		}
	}
	else
	{
		for(int i=0; i<iTextLen; i++)
		{
			if(pText[i] == 0)
			{
				pnode = pnode->lchild;
			}
			else
			{
				pnode = pnode->rchild;
			}

			if((pnode->lchild == nullptr) && (pnode->rchild == nullptr))
			{
				m_vecDeText.push_back(m_pElems[pnode->idx]);
				pnode = GetRoot();
			}
		}
	}

	int iDeTextLen = m_vecDeText.size();
	char * pDeText = new char[iDeTextLen+1];

	for(int i=0; i<iDeTextLen; i++)
	{
		pDeText[i] = m_vecDeText[i];
	}

	pDeText[iDeTextLen] = '\0';

	*ppOutput = pDeText;
	*pOutputLen = iDeTextLen;
	
	return iDeTextLen;
}

/** // test code 1
char g_text[] = "ADFHFAAAAHFGKKKKJJJJJJJJJJEEvkwwuuuuu";

void CHuffmanDemoDlg::OnBnClickedButton1()
{
	// TODO: 在此添加控件通知处理程序代码
	TRACE("text:%s\r\n", g_text); 
	CHuffmanCodec<char, int> hfmCodec;
	char * pOutput;
	int iOutputLen;
	int textlen = strlen(g_text);
	int len = hfmCodec.Encode(g_text, textlen, &pOutput, &iOutputLen);
	TRACE("encode:");
	for(int i=0; i<iOutputLen; i++)
	{
		TRACE("%d", pOutput[i]);
		if((1+i)%130==0)
			TRACE("\r\n");
	}

	TRACE("\r\ntext bit count: %d, after encoding:%d \r\n", textlen * 8, iOutputLen);

	char * pText;// = new char[1];
	int iTextLen;
	int len2 = hfmCodec.Decode(pOutput, iOutputLen, &pText, &iTextLen);
	TRACE("decode:%s\r\n", pText);

	delete[] pOutput;
	delete[] pText;
	
}

**/
/** // test code 2
char g_text[] = 
"BOOL CHuffmanDemoDlg::OnInitDialog()\r\n"\
"{\r\n"\
"	CDialogEx::OnInitDialog();\r\n"\
"\r\n"\
"	// 将“关于...”菜单项添加到系统菜单中。\r\n"\
"\r\n"\
"	// IDM_ABOUTBOX 必须在系统命令范围内。\r\n"\
"	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);\r\n"\
"	ASSERT(IDM_ABOUTBOX < 0xF000);\r\n"\
"\r\n"\
"	CMenu* pSysMenu = GetSystemMenu(FALSE);\r\n"\
"	if (pSysMenu != NULL)\r\n"\
"	{\r\n"\
"		BOOL bNameValid;\r\n"\
"		CString strAboutMenu;\r\n"\
"		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);\r\n"\
"		ASSERT(bNameValid);\r\n"\
"		if (!strAboutMenu.IsEmpty())\r\n"\
"		{\r\n"\
"			pSysMenu->AppendMenu(MF_SEPARATOR);\r\n"\
"			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);\r\n"\
"		}\r\n"\
"	}\r\n"\
"\r\n"\
"	// 设置此对话框的图标。当应用程序主窗口不是对话框时，框架将自动\r\n"\
"	//  执行此操作\r\n"\
"	SetIcon(m_hIcon, TRUE);			// 设置大图标\r\n"\
"	SetIcon(m_hIcon, FALSE);		// 设置小图标\r\n"\
"\r\n"\
"	// TODO: 在此添加额外的初始化代码\r\n"\
"\r\n"\
"	return TRUE;  // 除非将焦点设置到控件，否则返回 TRUE\r\n"\
"}\r\n";




void CHuffmanDemoDlg::OnBnClickedButton1()
{
	// TODO: 在此添加控件通知处理程序代码
/**	CElemStat<char> elemstat;
	int textlen = strlen(g_text);
	int elemnum = elemstat.Stat(g_text, strlen(g_text));

	char pElems[130];
	int  pCnts[130];
	int elemnum_ = elemstat.GetStat(pElems, pCnts, 130);

	for(int i=0; i<elemnum_; i++)
	{
		TRACE("%d: %c (%u) - %d\r\n", i+1, pElems[i], (BYTE)pElems[i], pCnts[i]);
	}
	
/**    CHuffman<int> huff;
    // int a[] = { 10,20,30,40,5,15 };
    // huff.CanonicCreat(a, 6);    //构建一棵哈夫曼树
    //int a[] = { 10,20,30,40,5,15 };
    //huff.creat(pCnts, elemnum_);    //构建一棵哈夫曼树
    huff.CanonicCreat(pCnts, elemnum_);    //构建一棵哈夫曼树
    //huff.print();        //打印节点间关系
	huff.printCode();
	char * pcode = nullptr;
	int codelen = 0;
	for(int j=0;j<elemnum_;j++)
	{
		bool bret = huff.getCode(j, &pcode, &codelen);
		for(int i=0; i<codelen; i++)
		{
			TRACE("%d", pcode[i]);
		}
		TRACE(" -- char: %c (%d) - cnt: %d 哈夫曼编码.\r\n", pElems[j], pElems[j], pCnts[j]);

		if(pcode != nullptr)
		{
			delete[] pcode;
			pcode = nullptr;
		}
	}
	huff.destroy();
	huff.ClearCodePtr();
	**
	
	CHuffmanCodec<char, int> hfmCodec;
	char * pOutput;
	int iOutputLen;
	int textlen = strlen(g_text);
	int len = hfmCodec.Encode(g_text, textlen, &pOutput, &iOutputLen);
	for(int i=0; i<iOutputLen; i++)
	{
		TRACE("%d", pOutput[i]);
		if((1+i)%130==0)
	TRACE("\r\n");
	}

	TRACE("\r\n");
	pOutput[14] = '\0';

	char * pText;// = new char[1];
	int iTextLen;
	int len2 = hfmCodec.Decode(pOutput, iOutputLen, &pText, &iTextLen);
	TRACE("%s\r\n", pText);

	delete[] pOutput;
	delete[] pText;
	
}
**/