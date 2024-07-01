#ifndef LITE_SPARSE_MATRIX_H
#define LITE_SPARSE_MATRIX_H

#include <list>
#include <vector>
#include <algorithm>
#include <cassert>
#include <ostream>
#include "Sparse_Entry.h"

// \addtogroup MathSuite 
//@{

//! Symmetric status 
enum SYMMETRIC_STATE
{
	NOSYM, /*!< general case */ 
	SYM_UPPER, /*!< symmetric (store upper triangular part) */
	SYM_LOWER, /*!< symmetric (store lower triangular part) */
	SYM_BOTH /*!< symmetric (store both upper and lower triangular part) */	
};

//!	Storage way
enum SPARSE_STORAGE
{
	CCS, /*!< compress column format */
	CRS, /*!< compress row format */
	TRIPLE /*!< row-wise coordinate format */
};

//! C or Fortran
enum ARRAYTYPE
{
	FORTRAN_TYPE, C_TYPE
};

//////////////////////////////////////////////////////////////////////////
//! Lite Sparse Matrix Class
class Lite_Sparse_Matrix
{
private:

	//!	Status for creating sparse solver
	enum STORE_STATE
	{
		ENABLE, 
		DISABLE, 
		LOCK	 
	};

	STORE_STATE state_fill_entry;
	SYMMETRIC_STATE sym_state; 
	SPARSE_STORAGE s_store;
	ARRAYTYPE arraytype;

	int nrows; //!< number of rows
	int ncols; //!< number of columns
	int nonzero; //!< number of nonzeros
	//! pointers to where columns begin in rowind and values 0-based, length is (col+1) 
	/*!
	* When s_store is CRS, colptr stores column indices;
	*/
	int* colptr; 
	//! row indices, 0-based
	/*!
	* When s_store is CRS, rowind stores row-pointers
	*/
	int* rowind; 
	double* values; //!< nonzero values of the sparse matrix

	std::vector< std::vector<Sparse_Entry> > entryset; //!< store temporary sparse entries

public:

	//! Sparse matrix constructor
	/*!
	* \param m row dimension
	* \param n column dimension
	* \param symmetric_state 
	* \param m_store the storage format
	* \param atype Fortran or C type of array
	*/
	Lite_Sparse_Matrix(int m, int n, SYMMETRIC_STATE symmetric_state = NOSYM, SPARSE_STORAGE m_store = CCS, ARRAYTYPE atype = C_TYPE)
		:nrows(m), ncols(n), sym_state(symmetric_state), s_store(m_store), arraytype(atype),
		nonzero(0), colptr(NULL), rowind(NULL), values(NULL), state_fill_entry(DISABLE)
	{
		if (m != n)
		{
			symmetric_state = NOSYM;
		}

		int nn = (m_store == CCS? ncols:nrows);
		entryset.resize(nn);
	}

	//! Sparse matrix destructor
	~Lite_Sparse_Matrix()
	{
		clear_mem();
	}
	//! Start to build sparse matrix pattern 
	inline void begin_fill_entry()
	{
		state_fill_entry = ENABLE;
	}
	//! Construct sparse pattern
	void end_fill_entry()
	{
		assert(state_fill_entry==ENABLE);

		clear_mem();

		state_fill_entry = LOCK;

		int inc =  (arraytype==FORTRAN_TYPE?1:0);

		if (s_store == CCS)
		{
			//construct map and ccs matrix
			int i, j, k = 0;
			colptr = new int[ncols+1];
			colptr[0] = inc;
			for (j=1; j<ncols+1; j++) {
				colptr[j] = (int)entryset[j-1].size() + colptr[j-1];
			}

			nonzero = colptr[ncols];

			if (nonzero > 0)
			{
				rowind = new int[nonzero];
				values = new double[nonzero];

				for ( j=0; j<ncols; j++) {
					for ( i=0; i<colptr[j+1]-colptr[j]; i++) {
						rowind[k] = entryset[j][i].index+inc;
						values[k] = entryset[j][i].value;
						k++;
					}
				}

			}

		}
		else if (s_store == CRS)
		{
			//construct map and crs matrix
			int i, j, k = 0;
			rowind = new int[nrows+1];
			rowind[0] = inc;
			for (j=1; j<nrows+1; j++) {
				rowind[j] = (int)entryset[j-1].size() + rowind[j-1];
			}
			nonzero = rowind[nrows];
			if (nonzero > 0)
			{
				colptr = new int[nonzero];
				values = new double[nonzero];

				for ( j = 0; j<nrows; j++) {
					for ( i = 0; i<rowind[j+1]-rowind[j]; i++) {
						colptr[k] = entryset[j][i].index+inc;
						values[k] = entryset[j][i].value;
						k++;
					}
				}
			}
		}
		else if (s_store == TRIPLE)
		{
			int i, j, k = 0;
			nonzero = 0;
			for (i = 0; i < nrows; i++)
			{
				nonzero += (int)entryset[i].size();
			}

			if (nonzero > 0)
			{
				rowind = new int[nonzero];
				colptr = new int[nonzero];
				values = new double[nonzero];

				for (i = 0; i < nrows; i++)
				{
					int jsize = (int)entryset[i].size();
					for (j = 0; j < jsize; j++)
					{
						rowind[k] = i+inc;
						colptr[k] = entryset[i][j].index+inc;
						values[k] = entryset[i][j].value;
						k++;
					}
				}
			}
		}
		entryset.clear();
	}

	//! Fill matrix entry \f$ 	Mat[rowind, colind] += val \f$
	void fill_entry(int row_index, int col_index, double val = 0 )
	{
		if(sym_state == NOSYM)
		{
			fill_entry_internal(row_index, col_index, val);
		}
		else if (sym_state == SYM_UPPER)
		{
			if (row_index <= col_index)
			{
				fill_entry_internal(row_index, col_index, val);
			}
			else
			{
				fill_entry_internal(col_index, row_index, val);
			}
		}
		else if (sym_state == SYM_LOWER)
		{
			if (row_index <= col_index)
			{
				fill_entry_internal(col_index, row_index, val);
			}
			else
			{
				fill_entry_internal(row_index, col_index, val);
			}
		}
		else if (sym_state == SYM_BOTH)
		{
			fill_entry_internal(row_index, col_index, val);

			if (row_index != col_index)
			{
				fill_entry_internal(col_index, row_index, val);
			}
		}
	}

	//! get the number of nonzeros
	inline int get_nonzero() const
	{
		return nonzero;
	}
	//! get the row dimension
	inline int rows() const
	{
		return nrows;
	}
	//! get the column dimension
	inline int cols() const
	{
		return ncols;
	}
	//! return the symmetric state
	inline bool issymmetric() const
	{
		return sym_state != NOSYM;
	}

	//! tell whether the matrix is upper or lower symmetric
	inline bool issym_store_upper_or_lower() const
	{
		return (sym_state == SYM_LOWER) || (sym_state == SYM_UPPER);
	}

	//! return symmetric state
	inline SYMMETRIC_STATE symmetric_state() const
	{
		return sym_state;
	}

	//! tell whether the matrix is square
	inline bool issquare() const
	{
		return nrows == ncols;
	}

	//! return the storage format
	inline SPARSE_STORAGE storage() const
	{
		return s_store;
	}

	//! return array type
	inline ARRAYTYPE get_arraytype() const
	{
		return arraytype;
	}

	//! get rowind
	inline int* get_rowind() const
	{
		return rowind;
	}

	//! get colptr
	inline int* get_colptr() const
	{
		return colptr;
	}

	//! get the values array
	inline double* get_values() const
	{
		return values;
	}

	//////////////////////////////////////////////////////////////////////////
private:
	//! Clear memory
	void clear_mem()
	{
		if (colptr != NULL) 
		{
			delete[] colptr;
			colptr = NULL;
		}

		if (rowind !=NULL) 
		{
			delete[] rowind;
			rowind = NULL;
		}

		if (values !=NULL) 
		{
			delete[] values;
			values = NULL;
		}
	}

	//! fill matrix entry (internal) \f$ Mat[rowid][colid] += val \f$
	bool fill_entry_internal(int row_index, int col_index, double val = 0)
	{
		assert(state_fill_entry==ENABLE);

		int search_index = (s_store==CCS?row_index:col_index);
		int pos_index = (s_store==CCS?col_index:row_index);

		Sparse_Entry forcompare(search_index);

		std::vector<Sparse_Entry>::iterator iter =
			std::lower_bound(entryset[pos_index].begin(), entryset[pos_index].end(), forcompare);
		if (iter != entryset[pos_index].end())
		{
			if (iter->index == search_index)
			{
				iter->value += val;
			}
			else
				entryset[pos_index].insert(iter, Sparse_Entry(search_index, val));
		}
		else
		{
			entryset[pos_index].push_back(Sparse_Entry(search_index, val));
		}
		return true;
	}
	//////////////////////////////////////////////////////////////////////////
};
//! print sparse matrix
std::ostream& operator<<(std::ostream &s, const Lite_Sparse_Matrix *A);

//@}


#endif //Lite_Sparse_Matrix_H

