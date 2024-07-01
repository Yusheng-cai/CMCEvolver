#include "Lite_Sparse_Matrix.h"

std::ostream& operator<<(std::ostream &s, const Lite_Sparse_Matrix *A)
{
	s.precision(16);
	if (A == NULL) {
		s<<"the matrix is null !\n ";
	}

	int row = A->rows();
	int col = A->cols();
	int nonzero = A->get_nonzero();
	int *rowind = A->get_rowind();
	int *colptr = A->get_colptr();
	double *values = A->get_values();

	s<<"row :"<<row<<" col :"<<col<<" Nonzero: "<<nonzero<<"\n\n";

	s<<"matrix --- (i, j, value)\n\n";

	SPARSE_STORAGE s_store = A->storage();
	int inc = (A->get_arraytype()==FORTRAN_TYPE?-1:0);
	if (s_store == CCS)
	{
		int k = 0;
		for (int i=1; i<col+1; i++) {
			for (int j=0; j<colptr[i]-colptr[i-1]; j++) {
				s<<rowind[k]+inc<<" "<<i-1<<" "<< std::scientific<<values[k]<<"\n";
				k++;
			}
		}
	}
	else if (s_store == CRS)
	{
		int k = 0;
		for (int i=1; i<row+1; i++) {
			for (int j=0; j<rowind[i]-rowind[i-1]; j++) {
				s<<i-1<<" "<<colptr[k]+inc<<" "<< std::scientific<<values[k]<<"\n";
				k++;
			}
		}
	}
	else if (s_store == TRIPLE)
	{
		for (int k = 0; k < nonzero; k++)
		{
			s << rowind[k]+inc <<" " << colptr[k]+inc << " " << std::scientific<<values[k]<<"\n";
		}
	}

	return s;
}