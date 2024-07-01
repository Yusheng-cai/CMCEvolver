
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// CMC Surface Evolver                                                        //
// http://i.cs.hku.hk/~hpan/code/cmc_source.zip        						  //
//                                                                            //
// CMC Surface Evolver is a demo implementation of the surface modeling method//
// proposed in the paper "Robust Modeling of Constant Mean Curvature Surfaces"//
//                                                                            //
// Author: Hao Pan (hpan@cs.hku.hk)                                           //
// Last Update: 2014-3-30                                                     //
//																			  //
// The software is freely available for non-commercial purposes.              //
//																			  //
////////////////////////////////////////////////////////////////////////////////

//The CMC_Evolver class is the interface of CMC surface computation. 
//Fast restricted Voronoi diagram is done by simply finding the Voronoi regions 
//of vertices over their 1-ring neighborhood facets.

#ifndef _FAST_RDT_
#define _FAST_RDT_

#include "edgemesh.h"

using namespace Geex;

namespace CMC {

	class FRDT_RestrictedVoronoiCell{
	public:
		double f_;
		vec3 g_;
		double mass_;
		vec3 vol_grad_;
		vec3 cvt_grad_;

		void reset(){
			f_ = 0;
			g_ = vec3(0,0,0);
			mass_ = 0;
			vol_grad_ = vec3(0,0,0);
			cvt_grad_ = vec3(0,0,0);
		}
	};

	class FRDT_RestrictedVoronoiDiagram{
	public:
		std::vector<FRDT_RestrictedVoronoiCell> rvcells_;
		EdgeMesh* p_mesh_;
		
	public:
		FRDT_RestrictedVoronoiDiagram():p_mesh_(NULL){}
		void compute_centroids(std::vector<vec3> *centroids, double vol_t);
		void get_fg(double &f, double *g, double vol_t = 0);
		
		//estimate t_i for each vertex.
		void estimate_ti();
	};

	class CMC_Evolver{
	public:
		EdgeMesh* edge_mesh_;
		FRDT_RestrictedVoronoiDiagram fast_rvd_;
		double vol_t_;
	protected:
		FRDT_RestrictedVoronoiDiagram* fast_rvd(){return &fast_rvd_;}
	public:
		CMC_Evolver();
		static CMC_Evolver* instance();

		void set_mesh(EdgeMesh* e_mesh);

		void estimate_t_i(bool delaunay = true);
		void set_volume_weight(double vol_t){ vol_t_ = vol_t;}

		void cmc_lloyd(int nb_iter, int out_itr_count);
		void cmc_qnewton(int nb_iter, int out_itr_count);
		void print_energy_cell_size();

		static CMC_Evolver* instance_ ;
	};

	//Edge flipping based Delaunay triangulation of surface mesh. This procedure handles topological changes.
	//For details, please refer to the paper "Robust Modeling of Constant Mean Curvature Surfaces".
	void frdt_delaunaize_mesh(EdgeMesh* mesh);
	//Compute energy and gradient data given a surface mesh.
	void frdt_RVD(EdgeMesh * mesh, FRDT_RestrictedVoronoiDiagram * frdt_rvd, double vol_t = 0);
}
#endif