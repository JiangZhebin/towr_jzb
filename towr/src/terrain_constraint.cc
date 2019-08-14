/******************************************************************************
Copyright (c) 2018, Alexander W. Winkler. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************/

#include <towr/constraints/terrain_constraint.h>
#include <iostream>

namespace towr {


TerrainConstraint::TerrainConstraint (std::string ee_motion,const SplineHolder& s,
                                      EE ee,
                                      const VecTimes& phase_duration, const double T)
    :ConstraintSet(kSpecifyLater, "terrain-" )
{
  ee_=ee;
  base_linear_sp=s.base_linear_;
  ee_motion_sp = s.ee_motion_.at(ee_);
  ee_motion_id_ = ee_motion;
  phase_duration_=phase_duration;
  total_time = T;
}

void
TerrainConstraint::InitVariableDependedQuantities (const VariablesPtr& x)
{
  ee_motion_ = x->GetComponent<NodesVariablesPhaseBased>(ee_motion_id_);

  // skip first node, b/c already constrained by initial stance
  for (int id=1; id<ee_motion_->GetNodes().size(); ++id)
    node_ids_.push_back(id);

  int constraint_count = node_ids_.size();
  SetRows(2*constraint_count);
std::cout<<constraint_count<<std::endl;
  double global = 0.0;
  node_time[0]=0;
  node_time[ee_motion_->GetNodes().size()-1]= total_time;
  for(int i=0;i<phase_duration_.size()-1;i++){
    global += phase_duration_.at(i);
    node_time[ee_motion_->GetNodeIDAtStartOfPhase(i+1)]=global;
  }
  for(int i=0; i<ee_motion_->GetNodes().size()/3;i++){
    node_time[3*i+2]=(node_time[3*i+1]+node_time[3*i+3])/2;

  }
  for(int i=0; i<ee_motion_->GetNodes().size();++i)
     std::cout<<"node: "<<i << " time: "<< node_time.at(i)<<std::endl;
}

Eigen::VectorXd
TerrainConstraint::GetValues () const
{
  VectorXd g(GetRows());

  auto nodes = ee_motion_->GetNodes();
  int row = 0;
  for (int id : node_ids_) {
   // double time_now = ToGlobalTime(ee_motion_,id);
    Vector3d p = nodes.at(id).p();
    g(row++) =  base_linear_sp->GetPoint(node_time.at(id)).p().z()-p.z()-0.1;
    g(row++) = base_linear_sp->GetPoint(node_time.at(id)).p().x()-p.x();
  }

  return g;
}

TerrainConstraint::VecBound
TerrainConstraint::GetBounds () const
{
  VecBound bounds(GetRows());
  double max_distance_above_terrain = 1; // [m]

  int row = 0;
  for (int id : node_ids_) {
    if (ee_motion_->IsConstantNode(id)){
      bounds.at(2*row) = ifopt::BoundZero;
      bounds.at(2*row+1) =ifopt::BoundZero;
    }
    else{
      bounds.at(2*row) = ifopt::Bounds(0.0, max_distance_above_terrain);
      bounds.at(2*row+1) = ifopt::NoBound;
    }
    row++;
  }

  return bounds;
}

void
TerrainConstraint::FillJacobianBlock (std::string var_set, Jacobian& jac) const
{
  if (var_set == ee_motion_->GetName()) {
    auto nodes = ee_motion_->GetNodes();
    int row = 0;
    for (int id : node_ids_) {
      int idx = ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(id, kPos, Z));
      jac.coeffRef(2*row, idx) = -1.0;
      int idz = ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(id,kPos,X));
      jac.coeffRef(2*row+1,idz) = -1.0;

//      Vector3d p = nodes.at(id).p();
//      for (auto dim : {X,Y}) {
//        int idx = ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(id, kPos, dim));
//        jac.coeffRef(row, idx) = -terrain_->GetDerivativeOfHeightWrt(To2D(dim), p.x(), p.y());
//      }
      row++;
    }
  }

//  if (var_set == "base-lin") {
//    auto nodes = ee_motion_->GetNodes();
//    int row = 0;
//    for (int id : node_ids_) {

//      int idx = ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(id, kPos, Z));
//      jac.coeffRef(row, idx) = 1.0;
//    }
//  }
}

double TerrainConstraint::ToGlobalTime(const NodesVariablesPhaseBased::Ptr ee_phasebased,
                                     const int node_id) const
{
  //get global time of phase nodes
  double t_global =0.0;
  int node_id_ = node_id;
  std::vector<double> num_of_nodes_in_a_phase;
  for(int i; i<phase_duration_.size();i++){
    if(i==phase_duration_.size())
      num_of_nodes_in_a_phase.push_back(ee_phasebased->GetNodes().size()-1-ee_phasebased->GetNodeIDAtStartOfPhase(i));
    else{
      num_of_nodes_in_a_phase.push_back(ee_phasebased->GetNodeIDAtStartOfPhase(i+1)-ee_phasebased->GetNodeIDAtStartOfPhase(i));

    }
  }

  std::vector<double> each_node_time;
  for(auto i=0;i<num_of_nodes_in_a_phase.size();i++){
    each_node_time.push_back(phase_duration_.at(i)/num_of_nodes_in_a_phase.at(i));
  }
  //node in which phase
  int in_which_phase;

  for(auto i=0; i<num_of_nodes_in_a_phase.size();i++)
  {
    node_id_ -= num_of_nodes_in_a_phase.at(i);
    if(node_id_<=0){
      in_which_phase=i;
      break;
    }
  }
 // std::cout<<node_id<<" in which phase: "<<in_which_phase<<std::endl;

  for(auto i=0; i<in_which_phase;i++)
    t_global += phase_duration_.at(i);

  //std::cout<<"t_global before"<<t_global<<std::endl;

  t_global += each_node_time.at(in_which_phase)*(node_id-ee_phasebased->GetNodeIDAtStartOfPhase(in_which_phase));
  return t_global;
}
} /* namespace towr */
