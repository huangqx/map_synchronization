#ifndef incre_spectral_threshold_h_
#define incre_spectral_threshold_h_
#include <vector>
using namespace std;

struct Matrix3D {
 public:
  Matrix3D() {
   m[0][0] = m[1][1] = m[2][2] = 1.0;
   m[0][1] = m[0][2] = m[1][2] = 0.0;
   m[1][0] = m[2][0] = m[2][1] = 0.0;
  }
  ~Matrix3D() {
  }
  double FNorm();
  double m[3][3];
};

struct Rotation3D {
 public:
  Rotation3D() {
   objectId = 0;
   m[0][0] = m[1][1] = m[2][2] = 1.0;
   m[0][1] = m[0][2] = m[1][2] = 0.0;
   m[1][0] = m[2][0] = m[2][1] = 0.0;
  }
  ~Rotation3D() {
  }
  unsigned objectId;
  double m[3][3];
};

struct PairwiseRotation3D {
 public:
  // Rst*Rs = Rt
  // Rst = Rt*Rs^{T}
  PairwiseRotation3D() {
    m[0][0] = m[1][1] = m[2][2] = 1.0;
	m[0][1] = m[0][2] = m[1][2] = 0.0;
	m[1][0] = m[2][0] = m[2][1] = 0.0;
	sourceObjectId = 0;
	targetObjectId = 0;
	weight = 0;
  }
  ~PairwiseRotation3D() {
  }
  unsigned sourceObjectId;
  unsigned targetObjectId;
  double m[3][3];
  double weight; // between 0 and 1 
};

struct Parameters {
 public:
  Parameters() {
   numIterations = 30;
   delta_max = 1.0;
   delta_min = 0.02;
   minimum_weight = 1e-8;
  }
  ~Parameters() {
  }
  unsigned numIterations;
  double delta_max;
  double delta_min;
  double minimum_weight;
};

class IteractiveSpecThre3D {
 public: 
  IteractiveSpecThre3D() {
   numObjects_ = 1000;
  }
  ~IteractiveSpecThre3D() {
  }
  // IO
  void SetNumObjects(unsigned num) {
    numObjects_ = num;
  }
  void Initialize(const double *edge_data,
	  // a matrix that specifies the topology of the graph
	  const double *rot_data,
	  // a matrix that specifies the pair-wise rotations
	  const unsigned &numPairs);
  vector<Rotation3D> *GetOptimalSolution() {
    return &sol_;
  }
  const vector<Rotation3D> *GetOptimalSolution() const {
	return &sol_;
  }
  void Run(const Parameters &para);
 private:
  void GeneralizedPowerMethod();
  
  // Set the weight of a rotation to minimum_weight, if 
  // \|R_{ij}*R_j - R_i\| > \delta
  void ThresholdingAndWeighting(const double &delta,
    const double &minimum_weight);
  void CalculateResidual(const PairwiseRotation3D &Rst,
    const Rotation3D &Rs,
    const Rotation3D &Rt,
    Matrix3D *residual);
 private:
  vector<PairwiseRotation3D> observation_;
  // The input pairwise rotations
  vector<Rotation3D> sol_;
  // The optimal solution to a world coordinate system
  unsigned numObjects_;
};

#endif
