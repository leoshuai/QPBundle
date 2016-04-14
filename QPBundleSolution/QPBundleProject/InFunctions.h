#pragma once
typedef Eigen::Ref<Eigen::VectorXd> Vec;
typedef const Eigen::Ref<const Eigen::VectorXd> cVec;
void CB3(cVec& x, double &y, Vec s);
void CB2(cVec& x,double &y,Vec s);
//void CB3(cVec& x,double &y,Vec s);
void DEM(cVec& x,double &y,Vec s);
void LQ(cVec& x,double &y,Vec s);
void QL(cVec& x,double &y,Vec s);
void Mifflin1(cVec& x,double &y,Vec s);
void Wolfe(cVec& x,double &y,Vec s);
void Chained_LQ(cVec& x,double &y,Vec s);

