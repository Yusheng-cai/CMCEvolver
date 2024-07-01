#ifndef SPARSE_ENTRY_H
#define SPARSE_ENTRY_H

//! Sparse Entry class \ingroup MathSuite
class Sparse_Entry
{
public:

	//! Index ID
	int index;

	//! Real value
	double value;

public:
	//! constructor
	Sparse_Entry(int ind, double v = 0)
		:index(ind), value(v)
	{}

	//! destructor
	~Sparse_Entry()
	{}

	//! The compare function for sorting
	inline bool operator< (const Sparse_Entry& m_r) const
	{
		return index < m_r.index;
	}
};

#endif //SPARSE_ENTRY_H
