
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

#include "fast_rdt.h"
#include "geex_utils.h"
#include <queue>
#include <set>
#include <vector>
#include "HLBFGS/HLBFGS.h"

using namespace Geex;

namespace CMC {

	//global variables
	
	CMC_Evolver * CMC_Evolver::instance_ = NULL;
	static std::vector<int> v_edges;
	static std::vector<bool> edge_marks ;
	static EdgeMesh temp_mesh;

	typedef struct {
		small_set<int,2> edge_vts;
		int edge_idx;
	} search_edge_t;

	struct sedge_comp {
		bool operator() (const search_edge_t& lhs, const search_edge_t& rhs) const
		{
			return lhs.edge_vts < rhs.edge_vts;	
		}
	};
	
	std::set<search_edge_t,sedge_comp> set_edges;

	//static variables
	static double time_for_DT = 0;
	static int frdt_curriter;
	static int frdt_call_nb_iter;
	

	int intersect(std::vector<int> & a, std::vector<int> & b){
		for (int i = 0; i < a.size(); i ++)
		{
			for (int j = 0; j < b.size(); j ++)
			{
				if (a[i] == b[j])
				{
					return a[i];
				}
			}
		}
		return -1;
	}


	void frdt_RVD(EdgeMesh * mesh, FRDT_RestrictedVoronoiDiagram * frdt_rvd, double vol_t /*= 0*/){
		
		frdt_rvd->p_mesh_ = mesh;
		
		frdt_rvd->rvcells_.resize(mesh->vertices_.size());

		for(int i=0; i< mesh->vertices_.size(); i++) {
			frdt_rvd->rvcells_[i].reset();
		}

		//compute CMC energy, gradient in two passes. First pass iterate through all facets. Second pass collect the values for all vertices.
		for (int i = 0; i < mesh->facets_.size(); i++)
		{
			
			/*compute the mid points, the centroid of the facet, that defines the voronoi region on this facet*/
			double facet_area, facet_energy;
			vec3 pts[3], mid_pt[3], fcentroid, n_edge[3];
			int pt_idx [3];
			double varea = 0;
			vec3 gradi (0,0,0);
			double bary_coeff[3];
			double pyramid_h;
			vec3 facet_area_normal;

			//condition check
			if (mesh->facets_[i].e_[0] == -1)
			{
				continue;
			}

			for (int j = 0; j < 3; j ++)
			{
				pt_idx[j] = mesh->edges_[mesh->facets_[i].e_[j]].intersect(mesh->edges_[mesh->facets_[i].e_[(j+1)%3]]);
				pts[j] = mesh->vertices_[pt_idx[j]].pos_;
			}
			
			for (int j = 0; j < 3; j ++)
			{
				mid_pt[j] = 0.5*(pts[(j+1)%3]+pts[(j+2)%3]);
				n_edge[j] = normalize( dot(pts[j] - pts[(j+1)%3], pts[(j+2)%3] - pts[(j+1)%3]) / (pts[(j+2)%3] - pts[(j+1)%3]).length2() 
					* (pts[(j+2)%3] - pts[(j+1)%3]) - (pts[j] - pts[(j+1)%3])) ;
			}

			facet_area = tri_area(pts[0],pts[1],pts[2]);
			if (facet_area < 1e-8)
			{
				continue;//degeneracy
			}
			facet_area_normal = 0.5 * cross(pts[2]-pts[0],pts[1]-pts[0]);
			pyramid_h = dot(normalize(facet_area_normal), pts[0]);

			for (int j = 0; j < 3; j ++)
			{
				bary_coeff[j] = (pts[(j+1)%3]-pts[(j+2)%3]).length2()*dot(pts[(j+1)%3] - pts[j], pts[(j+2)%3] - pts[j]);
			}

			fcentroid = 1.0/(8*facet_area*facet_area)*(bary_coeff[0]*pts[0]+bary_coeff[1]*pts[1]+bary_coeff[2]*pts[2]);

			/*compute the energy and gradient*/

			//assign t.
			pyramid_h *= vol_t;

			int obtuse_angle = -1;
			for (int j = 0; j < 3; j ++)
			{
				if (bary_coeff[j] < 0)
				{
					obtuse_angle = j;
					break;
				}
			}

			//compute
			if (obtuse_angle < 0)
			{
				//no obtuse angle
				for (int j = 0; j < 3; j ++)
				{
					facet_energy = Lloyd_energy(pts[j], pts[j], mid_pt[(j+1)%3], fcentroid, gradi, varea) ;
					frdt_rvd->rvcells_[pt_idx[j]].f_ += facet_energy + pyramid_h * varea ;
					frdt_rvd->rvcells_[pt_idx[j]].mass_ += varea;
					frdt_rvd->rvcells_[pt_idx[j]].g_ += gradi;
					frdt_rvd->rvcells_[pt_idx[j]].cvt_grad_ += gradi;

					facet_energy = Lloyd_energy(pts[j], pts[j], mid_pt[(j+2)%3], fcentroid, gradi, varea) ;
					frdt_rvd->rvcells_[pt_idx[j]].f_ += facet_energy + pyramid_h * varea ;
					frdt_rvd->rvcells_[pt_idx[j]].mass_ += varea;
					frdt_rvd->rvcells_[pt_idx[j]].g_ += gradi;
					frdt_rvd->rvcells_[pt_idx[j]].cvt_grad_ += gradi;
				}
			}
			else
			{
				//obtuse at obtuse_angle
				vec3 p_b_1, p_b_2;
				double tan_b_1 = 2 * (pts[obtuse_angle] - pts[(obtuse_angle+2)%3]).length2()*facet_area/bary_coeff[(obtuse_angle+1)%3],
					tan_b_2 = 2 * (pts[obtuse_angle] - pts[(obtuse_angle+1)%3]).length2()*facet_area/bary_coeff[(obtuse_angle+2)%3];
				p_b_1 = mid_pt[(obtuse_angle+2)%3] + tan_b_1 * 0.5 * 
					(pts[(obtuse_angle+1)%3] - pts[obtuse_angle]).length() * (- n_edge[(obtuse_angle+2)%3]);
				p_b_2 = mid_pt[(obtuse_angle+1)%3] + tan_b_2 * 0.5 * 
					(pts[(obtuse_angle+2)%3] - pts[obtuse_angle]).length() * (- n_edge[(obtuse_angle+1)%3]);

				facet_energy = Lloyd_energy(pts[obtuse_angle], pts[obtuse_angle], mid_pt[(obtuse_angle+2)%3], p_b_1, gradi, varea);
				frdt_rvd->rvcells_[pt_idx[obtuse_angle]].f_ += facet_energy + pyramid_h * varea ;
				frdt_rvd->rvcells_[pt_idx[obtuse_angle]].mass_ += varea;
				frdt_rvd->rvcells_[pt_idx[obtuse_angle]].g_ += gradi;
				frdt_rvd->rvcells_[pt_idx[obtuse_angle]].cvt_grad_ += gradi;

				facet_energy = Lloyd_energy(pts[obtuse_angle], pts[obtuse_angle], mid_pt[(obtuse_angle+1)%3], p_b_2, gradi, varea);
				frdt_rvd->rvcells_[pt_idx[obtuse_angle]].f_ += facet_energy + pyramid_h * varea ;
				frdt_rvd->rvcells_[pt_idx[obtuse_angle]].mass_ += varea;
				frdt_rvd->rvcells_[pt_idx[obtuse_angle]].g_ += gradi;
				frdt_rvd->rvcells_[pt_idx[obtuse_angle]].cvt_grad_ += gradi;

				facet_energy = Lloyd_energy(pts[obtuse_angle], pts[obtuse_angle], p_b_1, p_b_2, gradi, varea);
				frdt_rvd->rvcells_[pt_idx[obtuse_angle]].f_ += facet_energy + pyramid_h * varea ;
				frdt_rvd->rvcells_[pt_idx[obtuse_angle]].mass_ += varea;
				frdt_rvd->rvcells_[pt_idx[obtuse_angle]].g_ += gradi;
				frdt_rvd->rvcells_[pt_idx[obtuse_angle]].cvt_grad_ += gradi;

				facet_energy = Lloyd_energy(pts[(obtuse_angle+1)%3], pts[(obtuse_angle+1)%3], p_b_1, mid_pt[(obtuse_angle+2)%3], gradi, varea);
				frdt_rvd->rvcells_[pt_idx[(obtuse_angle+1)%3]].f_ += facet_energy + pyramid_h * varea ;
				frdt_rvd->rvcells_[pt_idx[(obtuse_angle+1)%3]].mass_ += varea;
				frdt_rvd->rvcells_[pt_idx[(obtuse_angle+1)%3]].g_ += gradi;
				frdt_rvd->rvcells_[pt_idx[(obtuse_angle+1)%3]].cvt_grad_ += gradi;

				facet_energy = Lloyd_energy(pts[(obtuse_angle+2)%3], pts[(obtuse_angle+2)%3], p_b_2, mid_pt[(obtuse_angle+1)%3], gradi, varea);
				frdt_rvd->rvcells_[pt_idx[(obtuse_angle+2)%3]].f_ += facet_energy + pyramid_h * varea ;
				frdt_rvd->rvcells_[pt_idx[(obtuse_angle+2)%3]].mass_ += varea;
				frdt_rvd->rvcells_[pt_idx[(obtuse_angle+2)%3]].g_ += gradi;
				frdt_rvd->rvcells_[pt_idx[(obtuse_angle+2)%3]].cvt_grad_ += gradi;
			}

			for (int j = 0; j < 3; j ++)
			{
				vec3 cvt_2nd_term = 1.0 / 24.0 * ((pts[(j+1)%3]-pts[j]).length()*(pts[(j+1)%3]-pts[j]).length2()*n_edge[(j+2)%3] +
					(pts[(j+2)%3]-pts[j]).length()*(pts[(j+2)%3]-pts[j]).length2()*n_edge[(j+1)%3] );
				frdt_rvd->rvcells_[pt_idx[j]].g_ += cvt_2nd_term;
				frdt_rvd->rvcells_[pt_idx[j]].cvt_grad_ += cvt_2nd_term;
				frdt_rvd->rvcells_[pt_idx[j]].g_ += vol_t * facet_area_normal;
				frdt_rvd->rvcells_[pt_idx[j]].vol_grad_ += facet_area_normal;
			}
		}
		
	}



	void frdt_delaunaize_mesh(EdgeMesh* mesh){

		int four_edges[4];
		double edge_length[4];

		set_edges.clear();

		for (int i = 0;i < mesh->edges_.size(); i ++)
		{
			v_edges.push_back(i);
			edge_marks[i] = true;

			search_edge_t nedge;
			nedge.edge_vts.insert(mesh->edges_[i].v_[0]);
			nedge.edge_vts.insert(mesh->edges_[i].v_[1]);
			nedge.edge_idx = i;
			set_edges.insert(nedge);
		}

		while(v_edges.size()!=0){
			int edge_idx = v_edges.back();
			v_edges.pop_back();
			edge_marks[edge_idx] = false;
			if (mesh->edges_[edge_idx].adj_facets_.size() == 2)
			{
				//get the other four edges
				int nb_edge_counter = 0;
				for (int i = 0; i < mesh->edges_[edge_idx].adj_facets_.size(); i ++)
				{
					int facet_idx = mesh->edges_[edge_idx].adj_facets_[i];
					for (int j = 0; j < 3; j ++)
					{
						if (mesh->facets_[facet_idx].e_[j] == edge_idx)
						{
							four_edges[nb_edge_counter++] = mesh->facets_[facet_idx].e_[(j+1)%3];
							four_edges[nb_edge_counter++] = mesh->facets_[facet_idx].e_[(j+2)%3];
							break;
						}
					}
				}
				if (mesh->edges_[four_edges[0]].intersect(mesh->edges_[four_edges[3]]) < 0)
				{
					//swap
					std::swap(four_edges[2],four_edges[3]);
				}
				//check delaunay
				for (int i = 0; i < 4; i ++)
				{
					EdgeMesh_Edge& edge = mesh->edges_[four_edges[i]];
					edge_length[i] = (mesh->vertices_[edge.v_[0]].pos_ - mesh->vertices_[edge.v_[1]].pos_).length();
				}
				double cur_edge_lgth = (mesh->vertices_[mesh->edges_[edge_idx].v_[0]].pos_ - 
					mesh->vertices_[mesh->edges_[edge_idx].v_[1]].pos_).length();
				double alpha = (edge_length[0]*edge_length[0] + edge_length[1]*edge_length[1] - cur_edge_lgth*cur_edge_lgth)
					/ (2 * edge_length[0] * edge_length[1]);
				double beta = (edge_length[2]*edge_length[2] + edge_length[3]*edge_length[3] - cur_edge_lgth*cur_edge_lgth)
					/ (2 * edge_length[2] * edge_length[3]);
				if (!is_nan(alpha + beta) && alpha + beta < 0)
				{
					//erase the current edge from set_edges.
					search_edge_t nedge_cur;
					nedge_cur.edge_vts.insert(mesh->edges_[edge_idx].v_[0]);
					nedge_cur.edge_vts.insert(mesh->edges_[edge_idx].v_[1]);
					set_edges.erase(set_edges.find(nedge_cur));


					int nv_1 = mesh->edges_[four_edges[0]].intersect(mesh->edges_[four_edges[1]]);
					int nv_2 = mesh->edges_[four_edges[2]].intersect(mesh->edges_[four_edges[3]]);

					int cur_edge_adj_facet_idx_0 = mesh->edges_[edge_idx].adj_facets_[0];
					int cur_edge_adj_facet_idx_1 = mesh->edges_[edge_idx].adj_facets_[1];

					//check if the flipped edge already exists
					search_edge_t nedge_flip;
					nedge_flip.edge_vts.insert(nv_1);
					nedge_flip.edge_vts.insert(nv_2);
					nedge_flip.edge_idx = edge_idx;
					std::pair<std::set<search_edge_t,sedge_comp>::iterator,bool> set_ret = set_edges.insert(nedge_flip);

					if (!set_ret.second)
					{ //this is an unflippable edge.
						int overlap_edge_idx = (*set_ret.first).edge_idx;


						int common_facet_idx = intersect(mesh->edges_[four_edges[0]].adj_facets_, mesh->edges_[four_edges[3]].adj_facets_);
						if (common_facet_idx > -1)
						{//the flipped facet already exists, clear duplicate facet.

							EdgeMesh_Facet& adj0_facet = mesh->facets_[cur_edge_adj_facet_idx_0];
							
							for (int i = 0; i < 3; i ++)
							{
								adj0_facet.e_[i] = -1;
							}

							//clear pointers in edges
							mesh->edges_[four_edges[0]].adj_facets_.erase(std::find(mesh->edges_[four_edges[0]].adj_facets_.begin(),mesh->edges_[four_edges[0]].adj_facets_.end(),cur_edge_adj_facet_idx_0));
							mesh->edges_[four_edges[3]].adj_facets_.erase(std::find(mesh->edges_[four_edges[3]].adj_facets_.begin(),mesh->edges_[four_edges[3]].adj_facets_.end(),cur_edge_adj_facet_idx_1));
						}
						else{

							mesh->edges_[overlap_edge_idx].adj_facets_.push_back(cur_edge_adj_facet_idx_0);

							(*std::find(mesh->edges_[four_edges[3]].adj_facets_.begin(),mesh->edges_[four_edges[3]].adj_facets_.end(),cur_edge_adj_facet_idx_1)) = cur_edge_adj_facet_idx_0;

							EdgeMesh_Facet& nf1 = mesh->facets_[cur_edge_adj_facet_idx_0];
							nf1.e_[0] = four_edges[0];nf1.e_[1] = overlap_edge_idx; nf1.e_[2] = four_edges[3];
						}

						common_facet_idx = intersect(mesh->edges_[four_edges[1]].adj_facets_,mesh->edges_[four_edges[2]].adj_facets_);
						if (common_facet_idx > -1)
						{
							//the flipped facet already exists, clear duplicate facet.
							EdgeMesh_Facet& adj1_facet = mesh->facets_[cur_edge_adj_facet_idx_1];

							for (int i = 0; i < 3; i ++)
							{
								adj1_facet.e_[i] = -1;
							}

							//clear pointers in edges
							mesh->edges_[four_edges[1]].adj_facets_.erase(std::find(mesh->edges_[four_edges[1]].adj_facets_.begin(),mesh->edges_[four_edges[1]].adj_facets_.end(),cur_edge_adj_facet_idx_0));
							mesh->edges_[four_edges[2]].adj_facets_.erase(std::find(mesh->edges_[four_edges[2]].adj_facets_.begin(),mesh->edges_[four_edges[2]].adj_facets_.end(),cur_edge_adj_facet_idx_1));
						}
						else{

							mesh->edges_[overlap_edge_idx].adj_facets_.push_back(cur_edge_adj_facet_idx_1);

							(*std::find(mesh->edges_[four_edges[1]].adj_facets_.begin(),mesh->edges_[four_edges[1]].adj_facets_.end(),cur_edge_adj_facet_idx_0)) = cur_edge_adj_facet_idx_1;

							EdgeMesh_Facet& nf2 = mesh->facets_[cur_edge_adj_facet_idx_1];
							nf2.e_[0] = overlap_edge_idx; nf2.e_[1] = four_edges[1]; nf2.e_[2] = four_edges[2];
						}
						
						//clear the edge.
						mesh->edges_[edge_idx].v_[0] = mesh->edges_[edge_idx].v_[1] = -1;
						mesh->edges_[edge_idx].adj_facets_.clear();

						//std::cout <<"unflippable edge.."<<std::endl;
					}
					else
					{
						//flip
						mesh->edges_[edge_idx].v_[0] = nv_1;
						mesh->edges_[edge_idx].v_[1] = nv_2;

						std::vector<int> & adj_facets = mesh->edges_[four_edges[3]].adj_facets_;
						for (int i = 0; i < adj_facets.size(); i ++)
						{
							if (adj_facets[i] == cur_edge_adj_facet_idx_1)
							{
								adj_facets[i] = cur_edge_adj_facet_idx_0;
								break;
							}
						}
						std::vector<int> & adj_facets_an = mesh->edges_[four_edges[1]].adj_facets_;
						for (int i = 0; i < adj_facets_an.size(); i ++)
						{
							if (adj_facets_an[i] == cur_edge_adj_facet_idx_0)
							{
								adj_facets_an[i] = cur_edge_adj_facet_idx_1;
								break;
							}
						}

						EdgeMesh_Facet& nf1 = mesh->facets_[cur_edge_adj_facet_idx_0];
						nf1.e_[0] = four_edges[0];nf1.e_[1] = edge_idx; nf1.e_[2] = four_edges[3];
						
						EdgeMesh_Facet& nf2 = mesh->facets_[cur_edge_adj_facet_idx_1];
						nf2.e_[0] = edge_idx; nf2.e_[1] = four_edges[1]; nf2.e_[2] = four_edges[2];

					}

					//push edges
					for (int i = 0; i < 4; i ++)
					{
						if (!edge_marks[four_edges[i]])
						{
							edge_marks[four_edges[i]] = true;
							v_edges.push_back(four_edges[i]);
						}
					}
				}
			}
		}

		//marking and cleaning.
		for (int i = 0; i < mesh->edges_.size(); i ++)
		{
			int v0,v1;
			v0 = mesh->edges_[i].v_[0];
			v1 = mesh->edges_[i].v_[1];
			if (mesh->edges_[i].adj_facets_.size() == 1 && mesh->vertices_[v0].type_ == COMMON_VTX && mesh->vertices_[v1].type_ == COMMON_VTX)
			{
				int adj_facet_idx = mesh->edges_[i].adj_facets_[0];
				EdgeMesh_Facet& adj_facet = mesh->facets_[adj_facet_idx];
				for (int j = 0; j < 3; j ++)
				{
					std::vector<int>::iterator v_itr = std::find(mesh->edges_[adj_facet.e_[j]].adj_facets_.begin(),mesh->edges_[adj_facet.e_[j]].adj_facets_.end(),adj_facet_idx);
					if (v_itr != mesh->edges_[adj_facet.e_[j]].adj_facets_.end())
					{
						mesh->edges_[adj_facet.e_[j]].adj_facets_.erase(v_itr);
					}

					adj_facet.e_[j] = -1;
				}

				mesh->edges_[i].adj_facets_.clear();
				mesh->edges_[i].v_[0] = mesh->edges_[i].v_[1] = -1;

				//std::cout << "clean 1 boundary edge." << std::endl;
			}
		}

	}

	void FRDT_RestrictedVoronoiDiagram::compute_centroids(std::vector<vec3> *centroids, double mc)
	{
		centroids->resize(rvcells_.size());
		int rvcell_itr;
#pragma omp parallel private(rvcell_itr) shared(centroids)
		{
#pragma omp for
			for (rvcell_itr = 0; rvcell_itr < rvcells_.size(); rvcell_itr ++)
			{
				FRDT_RestrictedVoronoiCell & cell = rvcells_[rvcell_itr];

				(*centroids)[rvcell_itr] = - cell.g_ / (2*cell.mass_);

				vec3 ctrd = (*centroids)[rvcell_itr];

				if (is_nan(ctrd.x) || is_nan(ctrd.y) || is_nan(ctrd.z) || p_mesh_->vertices_[rvcell_itr].type_ != COMMON_VTX)
				{
					(*centroids)[rvcell_itr] = vec3(0,0,0) ;
				}
			}
		}
	}

	void FRDT_RestrictedVoronoiDiagram::get_fg( double &f, double *g, double vol_t) 
	{
		f = 0;
		for (int i = 0; i <  rvcells_.size(); i ++)
		{
			f += rvcells_[i].f_;
			if (p_mesh_->vertices_[i].type_ != COMMON_VTX)
			{
				g[i*3] = g[i*3+1] = g[i*3+2] = 0;
			}
			else{
				g[i*3] = rvcells_[i].g_.x;
				g[i*3+1] = rvcells_[i].g_.y;
				g[i*3+2] = rvcells_[i].g_.z;
			}
		}
	}

	void FRDT_RestrictedVoronoiDiagram::estimate_ti()
	{

		int rvcell_itr;
		
		vec3 cvt_grad, vol_grad;
		int t_counter = 0;
		double avg_t = 0, max_t = -1e10, min_t = 1e10, fb_t, avg_cvt_grad = 0, avg_vol_grad = 0;

		for (rvcell_itr = 0; rvcell_itr < rvcells_.size(); rvcell_itr ++)
		{
			if (p_mesh_->vertices_[rvcell_itr].type_ != COMMON_VTX)
			{
				continue;
			}

			cvt_grad = rvcells_[rvcell_itr].cvt_grad_;
			vol_grad = rvcells_[rvcell_itr].vol_grad_;

			if (vol_grad.length() > 1e-8)
			{
				fb_t = -dot(cvt_grad,vol_grad)/vol_grad.length2();	
			}
			else{
				fb_t = 0;
			}

			if (fb_t < min_t)
			{
				min_t = fb_t;
			}
			if (fb_t > max_t)
			{
				max_t = fb_t;
			}
			avg_t += fb_t;
			avg_cvt_grad += cvt_grad.length();
			avg_vol_grad += vol_grad.length();
			t_counter ++;
			p_mesh_->vertices_[rvcell_itr].t_i_ = fb_t;
		}

		avg_t /= t_counter;
		std::cout<<"--- Feedback avg t is: " << avg_t << " min t is: " << min_t << " max t is: "<< max_t <<std::endl;

		std::cout.precision(16);
		avg_cvt_grad /= t_counter;
		avg_vol_grad /= t_counter;
		std::cout<<"Avg. ECVT gradient and avg. Vol gradient: \n"<< avg_cvt_grad << "\t" << avg_vol_grad << std::endl;

	}

	
	void CMC_Evolver::cmc_lloyd( int nb_iter, int out_itr_count)
	{
		std::vector<vec3> centroids;
		std::cout.precision(16);

		double *g, f0, f, step_lgth;

		centroids.resize(this->edge_mesh_->vertices_.size());
		g = new double[this->edge_mesh_->vertices_.size()*3];

		for(int i = 0; i < out_itr_count; i ++)
		{

			frdt_delaunaize_mesh(edge_mesh_);
			frdt_RVD(edge_mesh_, &fast_rvd_, vol_t_);
			f0 = 1; f = 0;
			for (int j = 0; j < nb_iter; j ++) // && (f0-f)/f0 > 0.0005
			{
				//evaluate gradient.
				fast_rvd_.get_fg(f0,g,vol_t_);
				fast_rvd_.compute_centroids(&centroids,vol_t_);

				//line search.
				step_lgth = 1;
				do 
				{
					temp_mesh = *edge_mesh_;
					for (int k = 0; k < temp_mesh.vertices_.size(); k ++)
					{
						temp_mesh.vertices_[k].pos_ += step_lgth * centroids[k];
					}
					//compute RVD
					frdt_delaunaize_mesh(&temp_mesh);
					frdt_RVD(&temp_mesh, &fast_rvd_, vol_t_);
					//get f & g
					fast_rvd_.get_fg(f,g,vol_t_);
					//update step length
					step_lgth *= 0.8;
				} while (f > f0);

				//update
				*edge_mesh_ = temp_mesh;
				std::cout << j <<": " << f << std::endl;
			}

		}

		delete g;
		
	}



	void frdt_funcgrad_cmc(int N, double* x, double* prev_x, double* f, double* g)
	{
		CMC_Evolver * frdt_ds = CMC_Evolver::instance();
		int counter = 0;
		
		temp_mesh = *(frdt_ds->edge_mesh_);
		for (int i = 0; i < temp_mesh.vertices_.size(); i ++)
		{
			temp_mesh.vertices_[i].pos_.x = x[3*counter];
			temp_mesh.vertices_[i].pos_.y = x[3*counter+1];
			temp_mesh.vertices_[i].pos_.z = x[3*counter+2];
			counter ++;
		}
		frdt_delaunaize_mesh(&temp_mesh);
		
		//compute RVD
		frdt_RVD(&temp_mesh, &(frdt_ds->fast_rvd_), frdt_ds->vol_t_);
		//get f & g
		frdt_ds->fast_rvd_.get_fg(*f,g,frdt_ds->vol_t_);

		std::cout << "f: " << *f << std::endl;
		
		frdt_call_nb_iter++;
	}

	void frdt_newiteration_cmc(int iter, int call_iter, double* x, double *f, double* g, double* gnorm)
	{

		CMC_Evolver * frdt_ds = CMC_Evolver::instance();

		int counter = 0;
		for (int i = 0; i < frdt_ds->edge_mesh_->vertices_.size(); i ++)
		{
			frdt_ds->edge_mesh_->vertices_[i].pos_.x = x[3*counter];
			frdt_ds->edge_mesh_->vertices_[i].pos_.y = x[3*counter+1];
			frdt_ds->edge_mesh_->vertices_[i].pos_.z = x[3*counter+2];
			counter ++;
		}

		frdt_delaunaize_mesh(frdt_ds->edge_mesh_);

		std::cout.precision(16);
		std::cout << frdt_curriter <<": " << frdt_call_nb_iter <<" " << std::scientific << *f  << " " << *gnorm  << std::endl;

		frdt_curriter ++;
		
	}

	void CMC_Evolver::cmc_qnewton( int nb_iter , int out_itr_count)
	{
		
		int nv = 0;
		for (int i = 0; i < edge_mesh_->vertices_.size(); i ++)
		{
			nv ++;
		}
		int n = nv * 3 ;
		int m = 7;

		//alloc and init the variables with mesh vertices
		double *x = new double[n];

		v_edges.clear();
		edge_marks.resize(edge_mesh_->edges_.size());

		//start the variational optimization process
		for (int opt_c = 0; opt_c < out_itr_count; opt_c ++)
		{

			int counter = 0;
			for(unsigned int i=0; i< edge_mesh_->vertices_.size(); i++) {
				x[3*counter ] = edge_mesh_->vertices_[i].pos_.x;
				x[3*counter+1] = edge_mesh_->vertices_[i].pos_.y;
				x[3*counter+2] = edge_mesh_->vertices_[i].pos_.z;
				counter ++;
			}

			//init the optimizer object

			frdt_curriter = 0;
			frdt_call_nb_iter = 0;

			double epsg = 1e-10;

			double parameter[20] ;
			int hlbfgs_info[20] ;

			//initialize parameters and infos
			INIT_HLBFGS(parameter, hlbfgs_info) ; 
			hlbfgs_info[3] = 0; //determines whether we use m1qn3
			hlbfgs_info[4] = nb_iter ; // max iterations
			hlbfgs_info[7] = 0; //without hessian
			hlbfgs_info[10] = 0; //determines whether we use cg
			parameter[5] = 0; // disabled
			parameter[6] = epsg;
			HLBFGS(
				n,
				m,
				x, 
				frdt_funcgrad_cmc,
				0, 
				HLBFGS_UPDATE_Hessian,
				frdt_newiteration_cmc,
				parameter,
				hlbfgs_info
				);

		}
		
		//clean allocated variables.
		delete [] x;

	}

	CMC_Evolver::CMC_Evolver()
	{
		instance_ = this;
		vol_t_ = 0;
	}

	CMC_Evolver* CMC_Evolver::instance()
	{
		return instance_;
	}

	void CMC_Evolver::estimate_t_i(bool delaunay /*=true*/)
	{
		if (delaunay)
		{
			frdt_delaunaize_mesh(edge_mesh_);
		}
		
		frdt_RVD(edge_mesh_, &fast_rvd_);
		fast_rvd_.estimate_ti();
		std::cout << "Assigning t_i to vertices...End!" << std::endl;
	}

	void CMC_Evolver::set_mesh(EdgeMesh* edge_mesh )
	{
		edge_mesh_ = edge_mesh;
		//delaunay vector initialization.
		v_edges.clear();
		edge_marks.resize(edge_mesh_->edges_.size());
	}

	void CMC_Evolver::print_energy_cell_size()
	{
		frdt_delaunaize_mesh(edge_mesh_);
		frdt_RVD(edge_mesh_, &fast_rvd_);

		std::ofstream ofile ("cell_energy_size.txt");
		for (int i = 0; i < fast_rvd_.rvcells_.size(); i ++)
		{
			ofile << fast_rvd_.rvcells_[i].f_ << "\t" << fast_rvd_.rvcells_[i].mass_ << std::endl;
		}
		ofile.close();
		std::cout<< "Energy and cell size saved." << std::endl;
	}
}
