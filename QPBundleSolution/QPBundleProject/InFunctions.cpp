#include "stdafx.h"
template<typename T>double one_norm(T& x)// One norm of a vector
	{
	double nm=0;
	for (int i=0;i<x.size();++i)
		nm=nm+abs(x[i]);
	return nm;
	}
using namespace Eigen;
using namespace std;
typedef Eigen::Ref<Eigen::VectorXd> Vec;
typedef const Eigen::Ref<const Eigen::VectorXd> cVec;

void CB2(cVec& x,double &y,Vec s)
	{
	VectorXd t(3);
	t<<pow(x(1),4)+pow(x(0),2),pow(2-x(0),2)+pow(2-x(1),2),2*exp(x(1)-x(0));
	int r,i,si;
	y=t.maxCoeff(&r);	
	vector<int> ind;
	ind.push_back(r);
	for (i=r+1;i<3;++i)
		{if (t(i)>=y)
		{
		ind.push_back(i);				
		}
		}
	si=ind.size();

	if (si==1)
		{
		switch (ind[0])
			{
			case 0:
				s<<2*x(0),4*pow(x(1),3);
				break;
			case 1:
				s<<-2*(2-x(0)),-2*(2-x(1));
				break;
			case 2:
				s<<-2*exp(-x(0)+x(1)),2*exp(-x(0)+x(1));
			}

		}
	else
		{
		if (si<=3)
			{
			MatrixXd subs(x.size(),si);
			for (i=0;i<si;++i)
				{
				switch (ind[i])
					{
					case 0:
						subs.col(i)<<2*x(0),4*pow(x(1),3);
						break;
					case 1:
						subs.col(i)<<-2*(2-x(0)),-2*(2-x(1));
						break;
					case 2:
						subs.col(i)<<-2*exp(-x(0)+x(1)),2*exp(-x(0)+x(1));

					}
				}
			vector<double> norm_selec_sub;
			for (i=0;i<si;++i)
				norm_selec_sub.push_back(one_norm(subs.col(i)));
			auto j=distance(norm_selec_sub.begin(),min_element(norm_selec_sub.begin(),norm_selec_sub.end()));
			s=subs.col(j);

			}
		else
			{
			random_device rd;
			mt19937 gen(rd());
			uniform_int_distribution<> dis(0, si-1);
			vector<int> ri(3);
			ri[0]=dis(gen);
			ri[1]=dis(gen);
			ri[2]=dis(gen);
			MatrixXd subs(x.size(),3);
			for (i=0;i<3;++i)
				{
				switch (ind[ri[i]])
					{
					case 0:
						subs.col(i)<<2*x(0),4*pow(x(1),3);
						break;
					case 1:
						subs.col(i)<<-2*(2-x(0)),-2*(2-x(1));
						break;
					case 2:
						subs.col(i)<<-2*exp(-x(0)+x(1)),2*exp(-x(0)+x(1));

					}
				}
			vector<double> norm_selec_sub;
			for (i=0;i<3;++i)
				norm_selec_sub.push_back(one_norm(subs.col(i)));
			auto j=distance(norm_selec_sub.begin(),min_element(norm_selec_sub.begin(),norm_selec_sub.end()));
			s=subs.col(j);

			}
		}
	}
void CB3(cVec& x, double &y, Vec s)
//void CB3(cVec& x,double &y,Vec s)
	{
	VectorXd t(3);
	t<<pow(x(1),2)+pow(x(0),4),pow(2-x(0),2)+pow(2-x(1),2),2*exp(x(1)-x(0));
	int r,i;
	size_t si;
	y=t.maxCoeff(&r);	
	vector<int> ind;
	ind.push_back(r);
	for (i=r+1;i<3;++i)
		{if (t(i)>=y)
		{
		ind.push_back(i);				
		}
		}
	si=ind.size();

	if (si==1)
		{
		switch (ind[0])
			{
			case 0:
				s<<4*pow(x(0),3),2*x(1);
				break;
			case 1:
				s<<-2*(2-x(0)),-2*(2-x(1));
				break;
			case 2:
				s<<-2*exp(-x(0)+x(1)),2*exp(-x(0)+x(1));
			}
		}
	else
		{
		if (si<=3)
			{
			MatrixXd subs(x.size(),si);
			for (i=0;i<si;++i)
				{
				switch (ind[i])
					{
					case 0:
						subs.col(i)<<4*pow(x(0),3),2*x(1);
						break;
					case 1:
						subs.col(i)<<-2*(2-x(0)),-2*(2-x(1));
						break;
					case 2:
						subs.col(i)<<-2*exp(-x(0)+x(1)),2*exp(-x(0)+x(1));

					}
				}
			vector<double> norm_selec_sub;
			for (i=0;i<si;++i)
				norm_selec_sub.push_back(one_norm(subs.col(i)));
			auto j=distance(norm_selec_sub.begin(),min_element(norm_selec_sub.begin(),norm_selec_sub.end()));
			s=subs.col(j);

			}
		else
			{
			random_device rd;
			mt19937 gen(rd());
			uniform_int_distribution<> dis(0, si-1);
			vector<int> ri(3);
			ri[0]=dis(gen);
			ri[1]=dis(gen);
			ri[2]=dis(gen);
			MatrixXd subs(x.size(),3);
			for (i=0;i<3;++i)
				{
				switch (ind[ri[i]])
					{
					case 0:
						subs.col(i)<<4*pow(x(0),3),2*x(1);
						break;
					case 1:
						subs.col(i)<<-2*(2-x(0)),-2*(2-x(1));
						break;
					case 2:
						subs.col(i)<<-2*exp(-x(0)+x(1)),2*exp(-x(0)+x(1));

					}
				}
			vector<double> norm_selec_sub;
			for (i=0;i<3;++i)
				norm_selec_sub.push_back(one_norm(subs.col(i)));
			auto j=distance(norm_selec_sub.begin(),min_element(norm_selec_sub.begin(),norm_selec_sub.end()));
			s=subs.col(j);

			}
		}
	}
void DEM(cVec& x,double &y,Vec s)
	{
	VectorXd t(3);
	t<<5*x(0)+x(1),-5*x(0)+x(1),pow(x(0),2)+pow(x(1),2)+4*x(1);
	int r,i,si;
	y=t.maxCoeff(&r);	
	vector<int> ind;
	ind.push_back(r);
	for (i=r+1;i<3;++i)
		{if (t(i)>=y)
		{
		ind.push_back(i);				
		}
		}
	si=ind.size();

	if (si==1)
		{
		switch (ind[0])
			{
			case 0:
				s<<5,1;
				break;
			case 1:
				s<<-5,1;
				break;
			case 2:
				s<<2*x(0),2*x(1)+4;
			}
		}
	else
		{
		if (si<=3)
			{
			MatrixXd subs(x.size(),si);
			for (i=0;i<si;++i)
				{
				switch (ind[i])
					{
					case 0:
						subs.col(i)<<5,1;
						break;
					case 1:
						subs.col(i)<<-5,1;
						break;
					case 2:
						subs.col(i)<<2*x(0),2*x(1)+4;

					}
				}
			vector<double> norm_selec_sub;
			for (i=0;i<si;++i)
				norm_selec_sub.push_back(one_norm(subs.col(i)));
			auto j=distance(norm_selec_sub.begin(),min_element(norm_selec_sub.begin(),norm_selec_sub.end()));
			s=subs.col(j);

			}
		else
			{
			random_device rd;
			mt19937 gen(rd());
			uniform_int_distribution<> dis(0, si-1);
			vector<int> ri(3);
			ri[0]=dis(gen);
			ri[1]=dis(gen);
			ri[2]=dis(gen);
			MatrixXd subs(x.size(),3);
			for (i=0;i<3;++i)
				{
				switch (ind[ri[i]])
					{
					case 0:
						subs.col(i)<<5,1;
						break;
					case 1:
						subs.col(i)<<-5,1;
						break;
					case 2:
						subs.col(i)<<2*x(0),2*x(1)+4;

					}
				}
			vector<double> norm_selec_sub;
			for (i=0;i<3;++i)
				norm_selec_sub.push_back(one_norm(subs.col(i)));
			auto j=distance(norm_selec_sub.begin(),min_element(norm_selec_sub.begin(),norm_selec_sub.end()));
			s=subs.col(j);

			}
		}
	}
void LQ(cVec& x,double &y,Vec s)
	//void f3(cVec& x,double &y,Vec s,int optional=0)// By default
	// both f(x) and a subgradient will be computed. In the future if I need to only evaluate f(x)
	// or a subgradient, I can add more options, if optional is 1 then just evaluate objective
	// if it is 2 then just evaluate subgradient.
	// But do this later because I need to get paper published asap. 
	// Still don't know how to unifying function and functors by using templates. 
	{
	double temp1 = -x(0) - x(1), temp2 = temp1 + (x(0)*x(0) + x(1)*x(1) - 1);
	//if (optional==1)
	//	y= max(temp1, temp2);
	//else if (optional==2)
	//{
	//	int n1=x.size();
	//	int n2=s.size();
	//	assert(n1==n2,  "Variable x and subgradient s should have same size"); // x and s should have same size. 1. write f3 then revise f,g,
	//	//and then revise the f and g in LPB function and the struct in Input 
	//	// and main. I searched the difference of C and advantage of C over C++. Conclusion
	//	//: just use C++!
	//	if (temp1 > temp2)
	//	{
	//		s << -1, -1;

	//	}

	//	if (temp1 < temp2)
	//	{
	//		s << 2 * x(0) - 1, 2 * x(1) - 1;

	//	}
	//	s(0) = (2 * x(0) - 1 - 1) / 2;
	//	s(1) = (2 * x(1) - 1 - 1) / 2;

	//}
	//else
	//{
	y= max(temp1, temp2);
	int n1=x.size();
	int n2=s.size();
	assert(n1==n2,  "Variable x and subgradient s should have same size"); // x and s should have same size. 1. write f3 then revise f,g,
	//and then revise the f and g in LPB function and the struct in Input 
	// and main. I searched the difference of C and advantage of C over C++. Conclusion
	//: just use C++!
	if (temp1 > temp2)
		{
		s << -1, -1;

		}

	else if (temp1 < temp2)
		{
		s << 2 * x(0) - 1, 2 * x(1) - 1;

		}
	else
		{s(0) = (2 * x(0) - 1 - 1) / 2;
	s(1) = (2 * x(1) - 1 - 1) / 2;}

	//}
	}
void QL(cVec& x,double &y,Vec s){
	vector<double> t(3);
	t[0]=pow(x[0],2)+pow(x[1],2);
	t[1]=t[0]+10*(-4*x[0]-x[1]+4);
	t[2]=t[0]+10*(-x[0]-2*x[1]+6);
	auto r=distance(t.begin(),max_element(t.begin(),t.end()));
	y=t[r];
	vector<int> ind;
	ind.push_back(r);
	for (auto i=r+1;i<3;++i)
		{if (t[i]>=y)
		{
		ind.push_back(i);				
		}
		}

	MatrixXd subs(x.size(),3);
	subs.col(0)<<2*x[0],2*x[1];
	subs.col(1)<<2*x[0]-40,2*x[1]-10;
	subs.col(2)<<2*x[0]-10,2*x[1]-20;
	if (ind.size()==1)
		{
		s= subs.col(ind[0]);
		}
	else
		{
		vector<double> norm_selec_sub;
		for (auto i=0;i<ind.size();++i)
			norm_selec_sub.push_back(one_norm(subs.col(ind[i])));
		auto j=distance(norm_selec_sub.begin(),min_element(norm_selec_sub.begin(),norm_selec_sub.end()));
		s= subs.col(ind[j]);
		}

	}
void Mifflin1(cVec& x,double &y,Vec s)
	{
	double t=pow(x(0),2)+pow(x(1),2)-1;
	y=-x(0)+20*max(t,0.0);

	if (t>0)
		{
		s(0)=-1+20*2*x(0);
		s(1)=20*2*x(1);
		}
	else 
		{
		s<<-1,0;
		}
	}
void Wolfe(cVec& x,double &y,Vec s){
	if (x(0)>abs(x(1)))
		{y=5*sqrt(9*pow(x(0),2)+16*pow(x(1),2));
	s<<18*x(0),16*2*x(1);
//std::cout<<s<<endl;
	//double temp=sqrt(9*pow(x(0),2)+16*pow(x(1),2));
//std::cout<<temp<<endl;
//double t2=5.0/2/temp;
//std::cout<<t2<<endl;
	//s=s*temp;
	//std::cout<<s<<endl;
	s*=5.0/2/sqrt(9*pow(x(0),2)+16*pow(x(1),2));}
	else if (0<x(0)&&x(0)<abs(x(1)))
		{y=9*x(0)+16*abs(x(1));
	auto dx2=0.0;
	if (x(1)>0)
		dx2=1;
	else if (x(1)<0)
		dx2=-1;
	s<<9,16*dx2;}
	else if (x(0)<0)
		{y=9*x(0)+16*abs(x(1))-pow(x(0),9);
	auto dx2=0.0;
	if (x(1)>0)
		dx2=1;
	else if (x(1)<0)
		dx2=-1;
	s<<9-9*pow(x(0),8),16*dx2;}
	else if (x(0)>0&&(x(1)>0||x(1)<0))
		{y=5*sqrt(9*pow(x(0),2)+16*pow(x(1),2));
	MatrixXd subs(s.size(),2);
	subs.col(0)<<18*x(0),16*2*x(1);
	subs.col(0)*=5.0/2/sqrt(9*pow(x(0),2)+16*pow(x(1),2));
	auto dx2=0.0;
	if (x(1)>0)
		dx2=1;
	else if (x(1)<0)
		dx2=-1;
	subs.col(1)<<9,16*dx2;
	vector<double> norm_selec_sub;	
	//norm_selec_sub.push_back(one_norm(subs.col(0)));
	norm_selec_sub.push_back(subs.col(0).norm());
	norm_selec_sub.push_back(subs.col(1).norm());
	auto j=distance(norm_selec_sub.begin(),min_element(norm_selec_sub.begin(),norm_selec_sub.end()));
	s=subs.col(j);
		}
	else if (x(1)>0||x(1)<0)
		{y=9*x(0)+16*abs(x(1))-pow(x(0),9);
	MatrixXd subs(s.size(),2);
	auto dx2=0.0;
	if (x(1)>0)
		dx2=1;
	else if (x(1)<0)
		dx2=-1;
	subs.col(1)<<9,16*dx2;
	subs.col(0)<<9-9*pow(x(0),8),16*dx2;
	vector<double> norm_selec_sub;	
	/*norm_selec_sub.push_back(one_norm(subs.col(0)));
	norm_selec_sub.push_back(one_norm(subs.col(1)));*/
	norm_selec_sub.push_back(subs.col(0).norm());
	norm_selec_sub.push_back(subs.col(1).norm());
	auto j=distance(norm_selec_sub.begin(),min_element(norm_selec_sub.begin(),norm_selec_sub.end()));
	s=subs.col(j);
		}
	else 
		{y=0.0;
	s<<9.0,0.0;
		}
	}
void Chained_LQ(cVec& x,double &y,Vec s)
	{
	int n = x.size();
	y=0.0;
	s.setZero();
	double temp;

for (auto i=0;i<n-1;i++)
    {temp=pow(x(i),2)+pow(x(i+1),2)-1;
    if (0<temp)
        {y=y+temp;
	s.segment(i,2)=s.segment(i,2)+2*x.segment(i,2);}
    
    y=y-x(i)-x(i+1);
	s(i)=s(i)-1;
	s(i+1)=s(i+1)-1;
   // s(i:i+1)=s(i:i+1)+[-1 -1]';
	}
}

