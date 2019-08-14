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

#include <cmath>
#include <iostream>

#include <towr/terrain/examples/height_map_examples.h>
#include <towr/nlp_formulation.h>
#include <ifopt/ipopt_solver.h>
#include <towr/Robot_arm_model.h>

#include <towr/matplotlibcpp.h>
namespace plt = matplotlibcpp;

using namespace towr;

// A minimal example how to build a trajectory optimization problem using TOWR.
//
// The more advanced example that includes ROS integration, GUI, rviz
// visualization and plotting can be found here:
// towr_ros/src/towr_ros_app.cc
int main()
{
  NlpFormulation formulation;

  formulation.model_.kinematic_model_ = std::make_shared<RobotArmKinematicModel>();
  formulation.model_.dynamic_model_ = std::make_shared<ObjectDynamicModel>();
  // set the initial position of the hopper
  formulation.initial_base_.lin.at(kPos).z() = 1.1;
  formulation.initial_ee_W_.push_back(Eigen::Vector3d(0,0,1));

  // define the desired goal state of the hopper
  formulation.final_base_.lin.at(towr::kPos) << 1.0, 0.0, 1.1;
  formulation.final_ee_W.lin.at(kPos) << 1.0,0.0,1.0;
  // Parameters that define the motion. See c'tor for default values or
  // other values that can be modified.
  // First we define the initial phase durations, that can however be changed
  // by the optimizer. The number of swing and stance phases however is fixed.
  // alternating stance and swing:     ____-----_____-----_____-----_____
  std::vector<double> time_schedule {0.2,0.4,0.2,0.4,0.2};

  formulation.params_.ee_phase_durations_.push_back(time_schedule);
  formulation.params_.ee_in_contact_at_start_.push_back(true);

  // Initialize the nonlinear-programming problem with the variables,
  // constraints and costs.
  ifopt::Problem nlp;
  SplineHolder solution;
  for (auto c : formulation.GetVariableSets(solution))
    nlp.AddVariableSet(c);
  for (auto c : formulation.GetConstraints(solution))
    nlp.AddConstraintSet(c);
  for (auto c : formulation.GetCosts())
    nlp.AddCostSet(c);

  // You can add your own elements to the nlp as well, simply by calling:
  // nlp.AddVariablesSet(your_custom_variables);
  // nlp.AddConstraintSet(your_custom_constraints);

  // Choose ifopt solver (IPOPT or SNOPT), set some parameters and solve.
  // solver->SetOption("derivative_test", "first-order");
  auto solver = std::make_shared<ifopt::IpoptSolver>();
  solver->SetOption("jacobian_approximation", "exact"); // "finite difference-values"
  solver->SetOption("max_cpu_time", 60.0);
  std::cout<<"before solve"<<std::endl;
  solver->Solve(nlp);

  double total_time=0;
  for(int i=0; i< time_schedule.size();++i)
    total_time += time_schedule.at(i);


  std::vector<double> tp(10000),f(10000),b_a(10000),b_l(10000),b_v(10000);
  tp.at(0)=0.0;
  f.at(0)=solution.ee_force_.at(0)->GetPoint(tp.at(0)).p().x();
  b_a.at(0)=solution.base_linear_->GetPoint(tp.at(0)).a().x();
  b_l.at(0)=solution.base_linear_->GetPoint(tp.at(0)).p().x();
  b_v.at(0)=solution.base_linear_->GetPoint(tp.at(0)).v().x();

  std::vector<double> f_z(10000), b_a_z(10000), b_l_z(10000),b_v_z(10000);
  f_z.at(0)=solution.ee_force_.at(0)->GetPoint(tp.at(0)).p().z();
  b_a_z.at(0)=solution.base_linear_->GetPoint(tp.at(0)).a().z();
  b_l_z.at(0)=solution.base_linear_->GetPoint(tp.at(0)).p().z();
  b_v_z.at(0)=solution.base_linear_->GetPoint(tp.at(0)).v().z();

  std::vector<double> c_p_z(10000), c_p_x(10000),c_p_y(10000);
  c_p_x.at(0) = solution.ee_motion_.at(0)->GetPoint(tp.at(0)).p().x();
  c_p_y.at(0) = solution.ee_motion_.at(0)->GetPoint(tp.at(0)).p().y();
  c_p_z.at(0) = solution.ee_motion_.at(0)->GetPoint(tp.at(0)).p().z();

  for(auto i=1;i<total_time*1000;i++){
    tp.at(i)=tp.at(i-1)+0.001;
    f.at(i)= solution.ee_force_.at(0)->GetPoint(tp.at(i)).p().x();
    b_a.at(i)= solution.base_linear_->GetPoint(tp.at(i)).a().x();
    b_l.at(i)= solution.base_linear_->GetPoint(tp.at(i)).p().x();
    b_v.at(i)= solution.base_linear_->GetPoint(tp.at(i)).v().x();

    f_z.at(i)= solution.ee_force_.at(0)->GetPoint(tp.at(i)).p().z();
    b_a_z.at(i)= solution.base_linear_->GetPoint(tp.at(i)).a().z();
    b_l_z.at(i)= solution.base_linear_->GetPoint(tp.at(i)).p().z();
    b_v_z.at(i)= solution.base_linear_->GetPoint(tp.at(i)).v().z();
       //std::cout<<"xxx"<<std::endl;
    c_p_x.at(i) = solution.ee_motion_.at(0)->GetPoint(tp.at(i)).p().x();
    c_p_y.at(i) = solution.ee_motion_.at(0)->GetPoint(tp.at(i)).p().y();
    c_p_z.at(i) = solution.ee_motion_.at(0)->GetPoint(tp.at(i)).p().z();

  }
//  std::cout<<"before plotting "<< std::endl;

  for(auto i=1; i<3;i++){
  plt::figure(i);
  plt::subplot(2,2,1);
  plt::title("EE-force");
  i==1? plt::plot(tp,f,"-r"): plt::plot(tp,f_z,"-r");
  plt::subplot(2,2,2);
  plt::title("Base Acceleration");
  i==1? plt::plot(tp,b_a,"-b"):plt::plot(tp,b_a_z,"-b");
  plt::subplot(2,2,3);
  plt::title("Base Velocity");
  i==1? plt::plot(tp,b_v,"-y"):plt::plot(tp,b_v_z,"-y");
  plt::subplot(2,2,4);
  plt::title("Base Position");
  i==1? plt::plot(tp,b_l,"-g"):plt::plot(tp,b_l_z,"-g");
  }

  plt::figure(3);
  plt::subplot(2,2,1);
  plt::title("ee_position_x");
  plt::plot(tp,c_p_x,"-r");
  plt::subplot(2,2,2);
  plt::title("ee_position_y");
  plt::plot(tp,c_p_y,"-b");
  plt::subplot(2,2,3);
  plt::title("ee_position_z");
  plt::plot(tp,c_p_z,"-g");
  // Can directly view the optimization variables through:
  // Eigen::VectorXd x = nlp.GetVariableValues()
  // However, it's more convenient to access the splines constructed from these
  // variables and query their values at specific times:
  using namespace std;
  cout.precision(2);
  nlp.PrintCurrent(); // view variable-set, constraint violations, indices,...
  cout << fixed;
  cout << "\n====================\nMonoped trajectory:\n====================\n";

  double t = 0.0;
  while (t<=solution.base_linear_->GetTotalTime() + 1e-5) {
    cout << "t=" << t << "\n";
    cout << "Base linear position x,y,z:   \t";
    cout << solution.base_linear_->GetPoint(t).p().transpose() << "\t[m]" << endl;

    cout << "Base Euler roll, pitch, yaw:  \t";
    Eigen::Vector3d rad = solution.base_angular_->GetPoint(t).p();
    cout << (rad/M_PI*180).transpose() << "\t[deg]" << endl;

    cout << "Foot position x,y,z:          \t";
    cout << solution.ee_motion_.at(0)->GetPoint(t).p().transpose() << "\t[m]" << endl;

    cout << "Contact force x,y,z:          \t";
    cout << solution.ee_force_.at(0)->GetPoint(t).p().transpose() << "\t[N]" << endl;

    bool contact = solution.phase_durations_.at(0)->IsContactPhase(t);
    std::string foot_in_contact = contact? "yes" : "no";
    cout << "Foot in contact:              \t" + foot_in_contact << endl;

    cout << endl;

    t += 0.2;
  }
  plt::show();
}
