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

//A simplified wing-edge datastructure for tri-mesh is defined here. Edge based structure enables
//efficient edge-flipping for Delaunay mesh construction.

#ifndef EDGE_BASED_MESH_DS
#define EDGE_BASED_MESH_DS

#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include "geex_utils.h"

// BEGIN YC_EDITS
#include "tools/CommonTypes.h"
// END YC_EDITS

using namespace Geex;

namespace CMC {
	using Real3 = CommonTypes::Real3;
	using INT3  = CommonTypes::index3;

	enum VertexType {COMMON_VTX, FIXED_VTX};

	class EdgeMesh_Vertex{
	public:
		vec3 pos_;
		int type_;
		double t_i_;

		EdgeMesh_Vertex(vec3 point){
			pos_ = point;
			type_ = VertexType::COMMON_VTX;
			t_i_ = 0;
		}
	};

	class EdgeMesh_Edge{
	public:
		int v_[2];
		std::vector<int> adj_facets_;
		EdgeMesh_Edge(int v1,int v2){
			v_[0] = v1;
			v_[1] = v2;
		}
		int intersect(EdgeMesh_Edge& rhs){
			if (v_[0] == rhs.v_[0] || v_[0] == rhs.v_[1])
			{
				return v_[0];
			}
			if (v_[1] == rhs.v_[0] || v_[1] == rhs.v_[1])
			{
				return v_[1];
			}
			return -1;
		}
		std::vector<int> intersect_facets(EdgeMesh_Edge& rhs){
			std::vector<int> ret_v;
			for (int i = 0; i < adj_facets_.size(); i ++)
			{
				for (int j = 0; j < rhs.adj_facets_.size(); j ++)
				{
					if (adj_facets_[i] == rhs.adj_facets_[j])
					{
						ret_v.push_back(adj_facets_[i]);
					}
				}
			}
			return ret_v;
		}
	};

	class EdgeMesh_Facet{
	public:
		int e_[3]; 

		EdgeMesh_Facet(){
			e_[0] = e_[1] = e_[2] = -1;
		}
	};

	class EdgeMesh
	{
	public:
		std::vector<EdgeMesh_Vertex> vertices_;
		std::vector<EdgeMesh_Edge> edges_;
		std::vector<EdgeMesh_Facet> facets_;

		EdgeMesh (){}

		EdgeMesh (const std::vector<Real3>& verts, const std::vector<INT3>& faces){
			typedef small_set<int, 2> edge_t;
			typedef small_set<int, 3> facet_t;
			std::vector<edge_t> v_edges; //vector of edges to examine
			std::vector<edge_t>::iterator e_idx;

			for (int i=0;i<verts.size();i++){
				vec3 p;
				p.x = verts[i][0]; p.y = verts[i][1]; p.z = verts[i][2];

				vertices_.push_back(p);
			}

			for (int i=0;i<faces.size();i++){
				edge_t e1,e2,e3;
				e1.insert(faces[i][0]); e1.insert(faces[i][1]);
				e2.insert(faces[i][1]); e2.insert(faces[i][2]);
				e3.insert(faces[i][2]); e3.insert(faces[i][0]);

				EdgeMesh_Facet n_facet;

				if ( (e_idx = std::find(v_edges.begin(), v_edges.end(), e1)) == v_edges.end() )
				{
					v_edges.push_back(e1);
					EdgeMesh_Edge n_edge (e1[0],e1[1]);
					n_edge.adj_facets_.push_back(facets_.size());
					edges_.push_back(n_edge);
					n_facet.e_[0] = edges_.size() - 1;
				}
				else{
					n_facet.e_[0] = e_idx - v_edges.begin();
					edges_[e_idx - v_edges.begin()].adj_facets_.push_back(facets_.size());
				}

				if ( (e_idx = std::find(v_edges.begin(), v_edges.end(), e2)) == v_edges.end() )
				{
					v_edges.push_back(e2);
					EdgeMesh_Edge n_edge (e2[0],e2[1]);
					n_edge.adj_facets_.push_back(facets_.size());
					edges_.push_back(n_edge);
					n_facet.e_[2] = edges_.size() - 1;
				}
				else{
					n_facet.e_[2] = e_idx - v_edges.begin();
					edges_[e_idx - v_edges.begin()].adj_facets_.push_back(facets_.size());
				}

				if ( (e_idx = std::find(v_edges.begin(), v_edges.end(), e3)) == v_edges.end() )
				{
					v_edges.push_back(e3);
					EdgeMesh_Edge n_edge (e3[0],e3[1]);
					n_edge.adj_facets_.push_back(facets_.size());
					edges_.push_back(n_edge);
					n_facet.e_[1] = edges_.size() - 1;
				}
				else{
					n_facet.e_[1] = e_idx - v_edges.begin();
					edges_[e_idx - v_edges.begin()].adj_facets_.push_back(facets_.size());
				}
				//n_facet.body_[0] = f.body_[0];
				//n_facet.body_[1] = f.body_[1];
				facets_.push_back(n_facet);
			}
		}

		//read obj file.
		EdgeMesh (std::string filename){

			typedef small_set<int, 2> edge_t;
			typedef small_set<int, 3> facet_t;

			std::vector<edge_t> v_edges; //vector of edges to examine
			std::vector<edge_t>::iterator e_idx;

			std::cout << "Load from file "<<filename<<std::endl;

			std::ifstream file (filename) ;
			std::vector<vec3> vertices ;
			
#define LINE_BUF_SIZE 512

			char line[LINE_BUF_SIZE];
			file.getline(line, LINE_BUF_SIZE);

			while (!file.eof()) {
				std::istringstream iss (line);
				
				std::string keyword ;
				iss >> keyword;
            
				if(keyword == "v") {
					vec3 p ;
					iss >> p ;
					vertices_.push_back(p) ;
			
				} else if(keyword == "f") {
					std::vector<int> cur_facet ;
					while(!iss.eof()) {
						std::string s ;
						iss >> s ;
						if(s.length() > 0) {
							std::istringstream v_input(s) ;
							int index ;
							v_input >> index ;
							if(index < 1 || index > (int)vertices_.size()) {
								std::cerr << "Out of bounds vertex index" << std::endl ;
							} else {
								cur_facet.push_back(index - 1) ;
							}
						}
					}

					if(cur_facet.size() != 3) {
						std::cerr << "facet with vertices other than 3, ignoring" << std::endl ;
					} else {
						//build edge and facet.
						edge_t e1,e2,e3;
						e1.insert(cur_facet[0]); e1.insert(cur_facet[1]);
						e2.insert(cur_facet[1]); e2.insert(cur_facet[2]);
						e3.insert(cur_facet[2]); e3.insert(cur_facet[0]);

						EdgeMesh_Facet n_facet;

						if ( (e_idx = std::find(v_edges.begin(), v_edges.end(), e1)) == v_edges.end() )
						{
							v_edges.push_back(e1);
							EdgeMesh_Edge n_edge (e1[0],e1[1]);
							n_edge.adj_facets_.push_back(facets_.size());
							edges_.push_back(n_edge);
							n_facet.e_[0] = edges_.size() - 1;
						}
						else{
							n_facet.e_[0] = e_idx - v_edges.begin();
							edges_[e_idx - v_edges.begin()].adj_facets_.push_back(facets_.size());
						}

						if ( (e_idx = std::find(v_edges.begin(), v_edges.end(), e2)) == v_edges.end() )
						{
							v_edges.push_back(e2);
							EdgeMesh_Edge n_edge (e2[0],e2[1]);
							n_edge.adj_facets_.push_back(facets_.size());
							edges_.push_back(n_edge);
							n_facet.e_[2] = edges_.size() - 1;
						}
						else{
							n_facet.e_[2] = e_idx - v_edges.begin();
							edges_[e_idx - v_edges.begin()].adj_facets_.push_back(facets_.size());
						}

						if ( (e_idx = std::find(v_edges.begin(), v_edges.end(), e3)) == v_edges.end() )
						{
							v_edges.push_back(e3);
							EdgeMesh_Edge n_edge (e3[0],e3[1]);
							n_edge.adj_facets_.push_back(facets_.size());
							edges_.push_back(n_edge);
							n_facet.e_[1] = edges_.size() - 1;
						}
						else{
							n_facet.e_[1] = e_idx - v_edges.begin();
							edges_[e_idx - v_edges.begin()].adj_facets_.push_back(facets_.size());
						}
						//n_facet.body_[0] = f.body_[0];
						//n_facet.body_[1] = f.body_[1];
						facets_.push_back(n_facet);

					} 
				}

				file.getline(line, LINE_BUF_SIZE);
			}
		}

		//save to .obj file.
		void save_obj(std::string filename){
			std::ofstream ofs (filename.c_str());

			for (int i = 0; i < vertices_.size(); i ++)
			{
				ofs << "v " << vertices_[i].pos_.x << " " << vertices_[i].pos_.y << " " << vertices_[i].pos_.z << std::endl;
			}

			for (int i = 0; i < facets_.size(); i ++)
			{
				if (facets_[i].e_[0] == -1)
				{
					continue;
				}

				ofs << "f "
					<< edges_[facets_[i].e_[0]].intersect(edges_[facets_[i].e_[1]]) + 1<< " " 
					<< edges_[facets_[i].e_[1]].intersect(edges_[facets_[i].e_[2]]) + 1<< " " 
					<< edges_[facets_[i].e_[2]].intersect(edges_[facets_[i].e_[0]]) + 1<< std::endl;
			}
			ofs.close();
			std::cout << filename <<" mesh saved. " << std::endl;
		}

		// BEGIN YC EDITS
		void get_verts_faces(std::vector<Real3>& verts, std::vector<INT3>& faces){
			verts.clear();
			faces.clear();

			for (int i=0;i<vertices_.size();i++){
				Real3 pos;
				pos[0] = vertices_[i].pos_.x; pos[1] = vertices_[i].pos_.y; pos[2] = vertices_[i].pos_.z;
				verts.push_back(pos);
			}

			for (int i=0;i<facets_.size();i++){
				if (facets_[i].e_[0] == -1){continue;}

				INT3 face;
				face[0] = edges_[facets_[i].e_[0]].intersect(edges_[facets_[i].e_[1]]);
				face[1] = edges_[facets_[i].e_[1]].intersect(edges_[facets_[i].e_[2]]);
				face[2] = edges_[facets_[i].e_[2]].intersect(edges_[facets_[i].e_[0]]);

				faces.push_back(face);
			}
		}

		// END YC EDITS



		//save to a self-defined format with edge infor. For faster loading.
		void save(std::string filename){
			std::ofstream ofs (filename.c_str());
			ofs << vertices_.size() << " " << edges_.size() << " " << facets_.size() << std::endl;
			for (int i = 0; i < vertices_.size(); i ++)
			{
				ofs << vertices_[i].pos_.x << " " << vertices_[i].pos_.y << " " << vertices_[i].pos_.z << std::endl;
			}
			for (int i = 0; i < edges_.size(); i ++)
			{
				ofs << edges_[i].v_[0] << " " << edges_[i].v_[1];
				ofs << " " << edges_[i].adj_facets_.size() ;
				for (int j = 0; j < edges_[i].adj_facets_.size(); j ++)
				{
					ofs << " " << edges_[i].adj_facets_[j];
				}
				ofs << std::endl;
			}
			for (int i = 0; i < facets_.size(); i ++)
			{
				ofs << facets_[i].e_[0] << " " << facets_[i].e_[1] << " " << facets_[i].e_[2] << std::endl;
			}
			ofs.close();
			std::cout << filename <<" tri mesh saved. " << std::endl;
		}

		//load from file of the self-defined format.
		void load(std::string filename){
			std::ifstream ifs(filename.c_str());
			int v_num, e_num, f_num;
			ifs >> v_num >> e_num >> f_num;
			
			double vx,vy,vz;
			for (int i = 0; i < v_num; i ++)
			{
				ifs >> vx >> vy >> vz;
				vertices_.push_back(EdgeMesh_Vertex(vec3(vx,vy,vz)));
			}

			int ev0,ev1,nf,nf_idx;
			for (int i = 0; i < e_num; i ++)
			{
				ifs >> ev0 >> ev1 >> nf;
				EdgeMesh_Edge edge (ev0,ev1);
				for (int j = 0; j < nf; j ++)
				{
					ifs >> nf_idx;
					edge.adj_facets_.push_back(nf_idx);
				}
				edges_.push_back(edge);
			}

			for (int i = 0; i < f_num; i ++)
			{
				EdgeMesh_Facet facet;
				ifs >> facet.e_[0] >> facet.e_[1] >> facet.e_[2];
				facets_.push_back(facet);
			}
			ifs.close();
			std::cout << filename << " tri mesh loaded." << std::endl;
		}

		void get_area_volume(double & area, double & volume){
			double total_volume = 0;
			double total_area = 0;

			for (int i = 0; i < facets_.size(); i ++)
			{
				int vtx_idx[3];
				vec3 vtx_pos[3];
				EdgeMesh_Facet& efacet = facets_[i];
				if (efacet.e_[0] < 0)
				{
					continue;
				}
				for (int j = 0; j < 3; j ++)
				{
					vtx_idx[j] = edges_[efacet.e_[j]].intersect(edges_[efacet.e_[(j+1)%3]]);
					vtx_pos[j] = vertices_[vtx_idx[j]].pos_;
				}

				double facet_volume = 0.5 * dot(cross(vtx_pos[2],vtx_pos[1]),vtx_pos[0]);
				double facet_area = 0.5*cross(vtx_pos[2]-vtx_pos[0],vtx_pos[1]-vtx_pos[0]).length();
				total_volume += facet_volume;
				total_area += facet_area;
			}

			area = total_area;
			volume = total_volume;
		}

		void clear(){
			vertices_.clear();
			edges_.clear();
			facets_.clear();
		}

		//mark boundary vertices as fixed.
		void mark_boundary(){
			for (int eitr = 0; eitr < edges_.size(); eitr++)
			{
				EdgeMesh_Edge& edge = edges_[eitr];
				if (edge.adj_facets_.size()<2)
				{
					vertices_[edge.v_[0]].type_ = VertexType::FIXED_VTX;
					vertices_[edge.v_[1]].type_ = VertexType::FIXED_VTX;
				}
			}
		}
	};

}

#endif