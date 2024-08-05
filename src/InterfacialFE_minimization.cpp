#include "InterfacialFE_minimization.h"

namespace MeshRefineStrategyFactory
{
    registry_<InterfacialFE_minimization> registerInterfacialFE_minimization("InterfacialFE_minimization");
}

// InterfacialFE_Solver::InterfacialFE_Solver(std::vector<INT3>& faces) : face_(faces){

// }

// double InterfacialFE_Solver::operator()(const Eigen::VectorXd& xd, Eigen::VectorXd& grad){
//     // convert the information to a mesh
//     std::vector<Real3> verts;
//     for (int i=0;i<3;i++){
//         Real3 p;
//         p[0] = xd[i*3 + 0];
//         p[1] = xd[i*3 + 1];
//         p[2] = xd[i*3 + 2];

//         verts.push_back(p);
//     }

//     Mesh m(verts, face_);

//     MeshTools::CalculateAreaDerivatives(Mesh &m, int &dAdpi)
// }

InterfacialFE_minimization::InterfacialFE_minimization(MeshRefineStrategyInput& input)
    : MeshRefineStrategy(input)
{
    pack_.ReadNumber("maxstep", ParameterPack::KeyType::Optional, maxStep_);
    pack_.ReadNumber("stepsize", ParameterPack::KeyType::Optional, stepsize_);
    pack_.ReadNumber("L", ParameterPack::KeyType::Optional, L_);
    pack_.ReadNumber("temperature", ParameterPack::KeyType::Optional, temperature_);
    pack_.ReadNumber("print_every", ParameterPack::KeyType::Optional, print_every);
    pack_.ReadNumber("tolerance", ParameterPack::KeyType::Optional, tol_);
    pack_.ReadNumber("optimize_every", ParameterPack::KeyType::Optional, optimize_every);
    pack_.Readbool("MaxStepCriteria", ParameterPack::KeyType::Optional, MaxStepCriteria);
    pack_.Readbool("boundaryMaxStepCriteria", ParameterPack::KeyType::Optional, boundaryMaxStepCriteria_);

    // boundary terms
    pack_.ReadNumber("boundarymaxstep", ParameterPack::KeyType::Optional, maxBoundaryStep_);
    pack_.ReadNumber("boundarystepsize", ParameterPack::KeyType::Optional,boundarystepsize_);
    pack_.ReadNumber("boundarytolerance", ParameterPack::KeyType::Optional, boundarytol_);
    pack_.ReadNumber("L2maxstep", ParameterPack::KeyType::Optional, maxL2step_);
    pack_.ReadNumber("boundary_optimize_every", ParameterPack::KeyType::Optional, boundary_optimize_every);
    pack_.ReadNumber("dgamma_gamma", ParameterPack::KeyType::Optional, dgamma_gamma_);
    pack_.ReadNumber("L2_guess", ParameterPack::KeyType::Optional, L2_);
    pack_.ReadNumber("zstar", ParameterPack::KeyType::Optional, zstar_);
    pack_.ReadNumber("L2_step_size", ParameterPack::KeyType::Optional, L2_stepsize_);
    pack_.ReadNumber("L2tolerance", ParameterPack::KeyType::Optional, L2tol_);
    pack_.ReadNumber("zstar_deviation", ParameterPack::KeyType::Optional, zstar_deviation_);
    pack_.Readbool("useNumerical", ParameterPack::KeyType::Optional, useNumerical_);
    pack_.Readbool("debug", ParameterPack::KeyType::Optional, debug_);
    pack_.Readbool("use_better_L2_update", ParameterPack::KeyType::Optional, use_better_L2_update_);

    // calculate mu and gamma based on temperature
    mu_ = CalculateMu(temperature_);
    gamma_ = CalculateGamma(temperature_);
}

void InterfacialFE_minimization::update_Mesh(Mesh& m){
    // obtain map from edge index {minIndex, maxIndex} to the face index 
    // obtain map from vertex index {index} to the Edge Index {minIndex, maxIndex}
    MeshTools::MapEdgeToFace(m, MapEdgeToFace_, MapVertexToEdge_);

    // Calculate the boundary vertices --> vertices which that has an edge shared by only 1 face
    MeshTools::CalculateBoundaryVertices(m, MapEdgeToFace_, boundaryIndicator_);

    // Calculate the vertex neighbors
    MeshTools::CalculateVertexNeighbors(m, neighborIndices_);

    // Map from vertex indices to face indices 
    MeshTools::MapVerticesToFaces(m, MapVertexToFace_);

    // Map from edges to their opposing vertex indices 
    MeshTools::MapEdgeToOpposingVertices(m, MapEdgeToFace_, MapEdgeToOpposingVerts_);

    // find number of vertices 
    const auto& v = m.getvertices();
    numVerts_ = v.size();
    newVertices_.clear();newVertices_.resize(numVerts_);
    TotalArea_.clear();TotalArea_.resize(numVerts_);
}

Real InterfacialFE_minimization::CalculateGamma(Real temperature){
    return (30.04 - 0.27477 * (270 - temperature)) * 1e-6 * 1e-18; //kJ/nm2
}

Real InterfacialFE_minimization::CalculateMu(Real temperature){
    Real mu = 0.021 * temperature - 5.695;

    // water ice free energy difference 
    return -mu;
}

void InterfacialFE_minimization::setSquaredGradients(Mesh& m){
    squared_gradients_.clear();
    const auto& vertices = m.getvertices();
    squared_gradients_.resize(vertices.size(), {0,0,0});
}


void InterfacialFE_minimization::refine(Mesh& mesh){
    // define mesh
    mesh_ = &mesh;
    mesh_->CalcVertexNormals();

    FE_.clear();

    // first generate necessary things for the mesh
    update_Mesh(*mesh_);

    Real3 Volume_shift={0,0,0};

    if (mesh_->isPeriodic()){
        Volume_shift = -0.5 * mesh_->getBoxLength();
    }

    // calculate cotangent weights 
    for (int i=0;i<maxStep_;i++){
        // calculate the derivatives dAdpi and dVdpi
        MeshTools::CalculateCotangentWeights(*mesh_, neighborIndices_, MapEdgeToFace_, MapEdgeToOpposingVerts_, dAdpi_);
        MeshTools::CalculateVolumeDerivatives(*mesh_, MapVertexToFace_, dVdpi_, Volume_shift);
    
        // obtain the vertices
        auto& verts = mesh_->accessvertices();

        // define some variables
        Real max=std::numeric_limits<Real>::lowest();
        Real avg_step = 0;
        int total_verts=0;

        // start updating
        #pragma omp parallel
        {
            Real local_max = std::numeric_limits<Real>::lowest();
            Real sum = 0;
            int n_verts= 0 ;
            #pragma omp for
            for (int j=0;j<numVerts_;j++){
                if (! MeshTools::IsBoundary(j, boundaryIndicator_)){
                    // calculate gradient 
                    Real3 gradient = dAdpi_[j] - dVdpi_[j] * rho_*(mu_ + L_) / gamma_;

                    // update position
                    verts[j].position_ = verts[j].position_ - stepsize_ * gradient;

                    // calculate size of step
                    Real step = std::sqrt(LinAlg3x3::DotProduct(gradient, gradient));
                    sum += step;
                    n_verts+=1;

                    if (step > local_max){
                        local_max = step;
                    }
                }
            }
           
            #pragma omp critical
            if (local_max > max){
                max = local_max;
            }

            #pragma omp critical
            avg_step += sum;

            #pragma omp critical
            total_verts += n_verts;
        }

        avg_step /= total_verts;

        // calculate the vertex normals 
        mesh_->CalcVertexNormals();

        if (MaxStepCriteria){if (max < tol_){break;}}
        else{if (avg_step < tol_){break;}}

        // print if necessary
        if ((i+1) % print_every == 0){
            std::vector<Real3> Normal;
            std::vector<Real> vecArea;
            Real a = MeshTools::CalculateArea(*mesh_, vecArea, Normal);
            Real V = MeshTools::CalculateVolumeDivergenceTheorem(*mesh_, vecArea, Normal, Volume_shift);
            Real E = a - rho_ * (mu_ + L_) / gamma_ * V;
            FE_.push_back(a - rho_ * (mu_ + L_) / gamma_ * V);
            std::cout << "At iteration " << i+1 << " Area = " << a << " " << " Volume = " << V << " Energy = " << E << std::endl;
            std::cout << "max step = " << max << std::endl;
            std::cout << "Avg step = " << avg_step << std::endl;
        }

        if ((i+1) % optimize_every == 0){
            MeshTools::CGAL_optimize_Mesh(*mesh_, 10, 60);
            MeshTools::ChangeWindingOrderPosZ(*mesh_);

            update_Mesh(*mesh_);
        }
    }

    mesh_->CalcVertexNormals();
}

void InterfacialFE_minimization::refineBoundary(Mesh& m, AFP_shape* shape){
    update_Mesh(m);
    need_update_ = false;

    area_list_.clear(); volume_list_.clear();
    Vnbs_list_.clear(); Anbs_list_.clear();
    Vunderneath_ = 0.0;
    Vnbs_underneath_ = 0.0;
    std::vector<Real> contact_angle_list, vecArea;
    std::vector<Real3> Normal;

    // first do we refine --> refine updates the mesh
    m_flatContact_ = m;
    this->refine(m_flatContact_);

    Real3 Volume_shift={0,0,0};

    if (m.isPeriodic()){
        Volume_shift = -0.5 * m.getBoxLength();
    }

    int num_L2_step=0;
    Real last_L2   =L2_;
    // outer loop for Lagrange refinement
    while (true) {
        Real mean_z;
        Real mean_ca;

        // solve for pi*
        int iteration = 0;
        int cont_ind  = 0;

        // set curr m to the orignal flat mesh --> also update all the corresponding things!
        Mesh curr_m = m_flatContact_;
        update_Mesh(curr_m);

        Real L2_g   = 0.0f;
        bool exceed_zstar_deviation=false;

        if (debug_){
            MeshTools::writePLY("debugBoundaryStart_" + std::to_string(num_L2_step) + ".ply" , curr_m);
        }

        // inner loop for pi, pib, L2 refinement
        while (true){
            // check if we are optimizing mesh
            if ((cont_ind+1) % boundary_optimize_every == 0){
                // then optimize mesh
                MeshTools::CGAL_optimize_Mesh(curr_m,10,60);
                MeshTools::ChangeWindingOrderPosZ(curr_m);
                update_Mesh(curr_m);
            }

                                        // boundary update //
            // calculate area derivative --> dAdr at the boundary and volume derivative dVdr 
            std::vector<Real3> dAdr ,dVdr;
            if (! MeshTools::CalculateCotangentWeights(curr_m, neighborIndices_, MapEdgeToFace_, MapEdgeToOpposingVerts_, dAdr)){
                MeshTools::CGAL_optimize_Mesh(curr_m, 10, 60);
                MeshTools::ChangeWindingOrderPosZ(curr_m);
                update_Mesh(curr_m);
                ASSERT(MeshTools::CalculateCotangentWeights(curr_m, neighborIndices_, MapEdgeToFace_, MapEdgeToOpposingVerts_, dAdr), "Something is seriously wrong");
            }
            MeshTools::CalculateVolumeDerivatives(curr_m, MapVertexToFace_, dVdr, Volume_shift);

            // calculate drdu, drdv, boundaryindices, dAnbsdv, dAnbsdu, dVnbsdv, dVnbsdu
            std::vector<int> BoundaryIndices;
            std::vector<Real2> dAnbsduv, dVnbsduv;
            std::vector<Real3> drdu, drdv;
            std::vector<Real> ulist, vlist;
            MeshTools::CalculatedAVnbsdUV(curr_m, shape, BoundaryIndices, ulist, \
                                    vlist, drdu, drdv, dAnbsduv, dVnbsduv, useNumerical_, Volume_shift);

            // access to the vertices in m
            contact_angle_list.clear();

            // set max boundary step to be lower initially = -3.40282e38
            std::vector<Real> boundarysteps;
            std::vector<Real> zlist;
            Real max_boundary_step  = std::numeric_limits<Real>::lowest();
            Real avg_boundary_step  = 0.0f;
            auto& verts             = curr_m.accessvertices();
            int N                   = BoundaryIndices.size();
            Real kk                 = rho_ * (L_ + mu_) / (2*gamma_);
            mean_z                  = 0.0;

            // pib refinement
            for (int j=0;j<BoundaryIndices.size();j++){
                // get the actual index of boundary
                int ind       = BoundaryIndices[j];

                Real dAdu     = LinAlg3x3::DotProduct(drdu[j], dAdr[ind]);
                Real dAdv     = LinAlg3x3::DotProduct(drdv[j], dAdr[ind]);

                Real dVdu     = LinAlg3x3::DotProduct(drdu[j], dVdr[ind]);
                Real dVdv     = LinAlg3x3::DotProduct(drdv[j], dVdr[ind]);

                // calculate dAnbsdu and dAnbsdv --> keep drdv the same 
                Real dAnbsdu  = dAnbsduv[j][0];
                Real dAnbsdv  = dAnbsduv[j][1];
                Real dVnbsdu  = dVnbsduv[j][0];
                Real dVnbsdv  = dVnbsduv[j][1];

                // calculate dEdv
                Real dEdv     = dAdv - rho_ * (L_ + mu_) / gamma_* (dVdv + dVnbsdv) + dgamma_gamma_ * dAnbsdv + L2_ * drdv[j][2] / (Real)N;
                if (cont_ind == 0){
                    L2_g      += (-dAdv + rho_ * (L_ + mu_) / gamma_ * (dVdv + dVnbsdv) - dgamma_gamma_ * dAnbsdv) / (drdv[j][2]);
                }

                // calculate the inverse jacobian
                auto invjac   = shape->InvJacobian(ulist[j], vlist[j],useNumerical_);
                Eigen::MatrixXd dEduv(2,1);
                dEduv << 0, dEdv;
                auto dEdr     = invjac.transpose() * dEduv;
                Real3 step; step[0]=dEdr(0,0); step[1]=dEdr(1,0); step[2]=dEdr(2,0);

                // we can calculate the contact angle by finding the dgamma_gamma where dEdv is 0
                Real ca = 1.0 / dAnbsdv * (-dAdv + rho_ * (L_ + mu_) / gamma_ * (dVdv + dVnbsdv));
                contact_angle_list.push_back(ca);

                // // update vlist 
                //Real3 new_p = shape->calculatePos(ulist[j], vlist[j] - boundarystepsize_ * dEdv);
                //Real3 shift_vec;
                //curr_m.CalculateShift(new_p, verts[ind].position_, shift_vec);
                //new_p = new_p + shift_vec;
                //verts[ind].position_ = new_p;
                verts[ind].position_ = verts[ind].position_ - boundarystepsize_ * step;

                // update the mean z
                mean_z += verts[ind].position_[2];
                zlist.push_back(verts[ind].position_[2]);

                // norm step 
                Real norm_step = LinAlg3x3::norm(step);
                boundarysteps.push_back(norm_step);

                if (norm_step > max_boundary_step){
                    max_boundary_step = norm_step;
                }

                avg_boundary_step += norm_step;
            }


            // calculate mean z
            mean_z    = mean_z / (Real)N;
            avg_boundary_step = avg_boundary_step / (Real)N;

            // check if we deviated too much from zstar
            if (std::abs(mean_z - zstar_) > zstar_deviation_){
                exceed_zstar_deviation = true;
                break;
            }

            // once we updated the boundary, let's check its perimeter and area again
            Real max_k = MeshTools::CalculateMaxCurvature(curr_m, BoundaryIndices);
            if (std::abs(kk / max_k) > 1.0){
                std::cout << "Not converged." << std::endl;
                break;
            }

            Real var  = Algorithm::calculateVariance(contact_angle_list);
            mean_ca   = Algorithm::calculateMean(contact_angle_list);
            Real zvar = Algorithm::calculateVariance(zlist);
            std::cout << "Mean z = " << mean_z << std::endl;
            std::cout << "z std = "  << std::sqrt(zvar) << std::endl;
            std::cout << "maxmimum boundary step = " << max_boundary_step << std::endl;
            std::cout << "Average boundary step = " << avg_boundary_step << std::endl;
            std::cout << "std of contact angle is " << std::sqrt(var) << std::endl;
            std::cout << "mean of contact angle is " << mean_ca << std::endl; 
            std::cout << "Using L2 = " << L2_ << std::endl;
            std::cout << "Zstar deviation = " << zstar_deviation_ << std::endl;

            if (debug_){
                MeshTools::writePLY("debugBoundary_" + std::to_string(num_L2_step) + "_" + std::to_string(cont_ind) + ".ply" , curr_m);
                std::cout << "we are now going to optimize debugBoundary " << std::to_string(num_L2_step) << "_" << std::to_string(cont_ind) << "\n"; 
            }

            // then do pi refinement --> again
            this->refine(curr_m);

            if (iteration > maxBoundaryStep_){
                std::cout << "Exceed max boundary step." << std::endl;
                break;
            }

            if (boundaryMaxStepCriteria_){
                if (max_boundary_step < boundarytol_){
                    break;
                }
            }
            else{
                if (avg_boundary_step < boundarytol_){
                    break;
                }
            }

            // update iterations
            cont_ind++;
            iteration++;
        }

        Real L2_step = (mean_z - zstar_);
        L2_list_.push_back(L2_);

        // break the while loop
        if (std::abs(L2_step) < L2tol_ || num_L2_step > maxL2step_){
            m = curr_m;
            break;
        }

        // if we are on the first step --> then we make a better guess
        if (num_L2_step == 0 && use_better_L2_update_){
            if (exceed_zstar_deviation){
                L2_ = L2_g;
            }
        }
        else{
            Real updated_L2 = L2_ + L2_stepsize_ * L2_step;
            if (std::abs(updated_L2 - last_L2) < L2tol_ * L2_stepsize_){
                updated_L2 = L2_ + L2_stepsize_ * 0.5 * L2_step;
            }

            last_L2 = L2_;
            L2_     = updated_L2;
        }
        
        num_L2_step++;
    }

    // get the information about the boundary indices 
    std::vector<int> BoundaryIndices;
    std::vector<Real> ulist, vlist;
    MeshTools::CalculateBoundaryVerticesIndex(m, BoundaryIndices,true);
    MeshTools::FindBoundaryUV(m, ulist, vlist, BoundaryIndices, shape);

    // calculate area and volume after boundary steps
    vecArea.clear(); Normal.clear();
    a_ = MeshTools::CalculateArea(m, vecArea, Normal);
    V_ = MeshTools::CalculateVolumeDivergenceTheorem(m, vecArea, Normal);
    MeshTools::CalculateAVnbs(m, shape, BoundaryIndices,
                                ulist, vlist, Anbs_, Vnbs_,10000, useNumerical_, Volume_shift);

    area_list_.push_back(a_);
    volume_list_.push_back(V_);
    Vnbs_list_.push_back(Vnbs_);
    Anbs_list_.push_back(Anbs_);

    Vunderneath_     = MeshTools::CalculateVolumeUnderneath(m, 2);
    Vnbs_underneath_ = MeshTools::CalculateVnbsUnderneath(m, shape, 2, 10000, useNumerical_); 

    need_update_ = true;
}