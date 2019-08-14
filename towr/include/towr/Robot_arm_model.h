#include <towr/models/kinematic_model.h>
#include <towr/models/single_rigid_body_dynamics.h>

namespace towr{

class RobotArmKinematicModel : public KinematicModel{
public:
  RobotArmKinematicModel () : KinematicModel(1)
  {
    nominal_stance_.at(0) = Eigen::Vector3d(0,0,-0.1);
    max_dev_from_nominal_ << 1, 1, 1;

  }
  ~RobotArmKinematicModel(){};
};

class ObjectDynamicModel : public SingleRigidBodyDynamics{
public:
  ObjectDynamicModel(): SingleRigidBodyDynamics(10,0.16,0.16,0.16,0,0,0,1)

  {}
  ~ObjectDynamicModel(){};
};
}
