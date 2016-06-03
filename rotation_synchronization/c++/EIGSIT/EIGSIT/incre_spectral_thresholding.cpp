#include "incre_spectral_thresholding.h"

void IteractiveSpecThre3D::Initialize(
  const double *edge_data,
  const double *rot_data,
  const unsigned &numPairs) {
  observation_.resize(numPairs);
  numObjects_ = 0;
  for (unsigned pairId = 0; pairId < numPairs; ++pairId) {
	  PairwiseRotation3D *rot = &observation_[pairId];
	  rot->sourceObjectId = static_cast<unsigned> (edge_data[2*pairId]   - 0.5);
	  rot->targetObjectId = static_cast<unsigned> (edge_data[2*pairId + 1] - 0.5);
	  for (int i = 0; i < 3; ++i) {
	    for (int j = 0; j < 3; ++j) {
		    rot->m[i][j] = rot_data[9*numPairs+3*j+i];
	    }
	  }
	  if (numObjects_ <= rot->sourceObjectId)
      numObjects_ = rot->sourceObjectId;
    if (numObjects_ <= rot->targetObjectId)
      numObjects_ = rot->targetObjectId;
  }
  numObjects_++;
}

void IteractiveSpecThre3D::Run(const Parameters &para) {

}

void IteractiveSpecThre3D::GeneralizedPowerMethod() {
  vector<Rotation3D> sol_prev = sol_;
  for (unsigned iter = 0; iter < 32; ++iter) {
    // Initialize
    for (unsigned vId = 0; vId < numObjects_; ++vId) {
      for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j)
          sol_[vId].m[i][j] = 0.0;
    }

    sol_prev = sol_;
  }

}

void IteractiveSpecThre3D::ThresholdingAndWeighting(const double &delta,
  const double &minimum_weight) {
  Matrix3D residual;
  for (unsigned pairId = 0; pairId < observation_.size(); ++pairId) {
    PairwiseRotation3D *rot = &observation_[pairId];
    const Rotation3D &Rs = sol_[rot->sourceObjectId];
    const Rotation3D &Rt = sol_[rot->targetObjectId];
    CalculateResidual(*rot, Rs, Rt, &residual);
    if (residual.FNorm() < delta)
      rot->weight = 1.0;
    else
      rot->weight = minimum_weight; 
  }
  
  // Apply weighting
  vector<double> degree;
  degree.resize(numObjects_);
  for (unsigned vId = 0; vId < numObjects_; ++vId)
    degree[vId] = 0.0;

  for (unsigned pairId = 0; pairId < observation_.size(); ++pairId) {
    const PairwiseRotation3D &rot = observation_[pairId];
    degree[rot.sourceObjectId] += rot.weight;
    degree[rot.targetObjectId] += rot.weight;
  }
  for (unsigned vId = 0; vId < numObjects_; ++vId)
    degree[vId] = sqrt(1 / degree[vId]);

  for (unsigned pairId = 0; pairId < observation_.size(); ++pairId) {
    PairwiseRotation3D *rot = &observation_[pairId];
    rot->weight = rot->weight
      *degree[rot->sourceObjectId]*degree[rot->targetObjectId];
  }
}

void IteractiveSpecThre3D::CalculateResidual(
  const PairwiseRotation3D &Rst,
  const Rotation3D &Rs,
  const Rotation3D &Rt,
  Matrix3D *residual) {
  residual->m[0][0] = Rt.m[0][0] - Rst.m[0][0]*Rt.m[0][0] - Rst.m[0][1]*Rt.m[1][0] - Rst.m[0][2]*Rt.m[2][0];
  residual->m[0][1] = Rt.m[0][1] - Rst.m[0][0]*Rt.m[0][1] - Rst.m[0][1]*Rt.m[1][1] - Rst.m[0][2]*Rt.m[2][1];
  residual->m[0][2] = Rt.m[0][2] - Rst.m[0][0]*Rt.m[0][2] - Rst.m[0][1]*Rt.m[1][2] - Rst.m[0][2]*Rt.m[2][2];
  residual->m[1][0] = Rt.m[1][0] - Rst.m[1][0]*Rt.m[0][0] - Rst.m[1][1]*Rt.m[1][0] - Rst.m[1][2]*Rt.m[2][0];
  residual->m[1][1] = Rt.m[1][1] - Rst.m[1][0]*Rt.m[0][1] - Rst.m[1][1]*Rt.m[1][1] - Rst.m[1][2]*Rt.m[2][1];
  residual->m[1][2] = Rt.m[1][2] - Rst.m[1][0]*Rt.m[0][2] - Rst.m[1][1]*Rt.m[1][2] - Rst.m[1][2]*Rt.m[2][2];
  residual->m[2][0] = Rt.m[2][0] - Rst.m[2][0]*Rt.m[0][0] - Rst.m[2][1]*Rt.m[1][0] - Rst.m[2][2]*Rt.m[2][0];
  residual->m[2][1] = Rt.m[2][1] - Rst.m[2][0]*Rt.m[0][1] - Rst.m[2][1]*Rt.m[1][1] - Rst.m[2][2]*Rt.m[2][1];
  residual->m[2][2] = Rt.m[2][2] - Rst.m[2][0]*Rt.m[0][2] - Rst.m[2][1]*Rt.m[1][2] - Rst.m[2][2]*Rt.m[2][2];
}


double Matrix3D::FNorm() {
  double sqrNorm0 = m[0][0]*m[0][0] + m[0][1]*m[0][1] + m[0][2]*m[0][2];
  double sqrNorm1 = m[1][0]*m[1][0] + m[1][1]*m[1][1] + m[1][2]*m[1][2];
  double sqrNorm2 = m[2][0]*m[2][0] + m[2][1]*m[2][1] + m[2][2]*m[2][2];
  return sqrt(sqrNorm0 + sqrNorm1 + sqrNorm2);
}