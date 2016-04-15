#pragma once
struct ProbData
	{
	std::string Prob_Name;
	double f_optimal;
	};
struct OutData
	{
	int status;// Status of output 
	double Error;// For public version, I don't need this because I will not know the optimal value of all the problems.
	double f_final;// The approximate optimal value returned by LPB
	int No_func_eval;// Number of function evaluation
	int k;// Number of all iterations
	int L;// Number of null steps
	double time;// Time spent in solving a problem
	double t_CPX;// Total time spent in CPLEX
	};
class Param
	{
	public:
		double epsilon_tol;
		double m;
		int Iter_Limi;

		Param():
			epsilon_tol(1.14e-4),
			m(1e-4),
			Iter_Limi(200)
			{}
		Param(const double eps_tol,const double mp,const unsigned int itli):
			epsilon_tol(eps_tol),
			m(mp),
			Iter_Limi(itli)
			{}
	};
template<typename T> OutData LPB(ProbData prob,const Eigen::Ref<const Eigen::VectorXd>& x,T &feval)
	//template<class X,typename T> OutData LPB(ProbData prob,X &x,T &feval)
	// Later optional input 
	// including options,pass value by reference.
	// Later I need to know the difference of explicit instantiation and implicit instantiation
	// and when to use them. // later I can set "Ref<VectorXd> x" to be "C x" where "C" is some class
	// so LPB can take vector<double> as arguments. This needs to be templated.
	// done! But inside this funciton there is still a line : VectorXd xStar(n),s(n);
	// how can I change it to X xStar,s;? I just did it and it works for Ref<Vector> 
	// I assume if clients only use vector<double> then it can also work as
	// long as it supports feval(X &x,double &y,S s) where X is type for x and S is type for 
	// subgradient vector. No! I just realized that there is inner product using eigen, and that X should 
	// be const but in the line there is a usage of X which means non-const type.
	// Conclusion: I need to move on for now and later if I can change this.
	// OLD version is template<typename T> OutData LPB(ProbData prob,const Ref<const VectorXd> x,T &feval)
	// Shuai Liu at 24/09/2015 6:00 PM
	// declaring a variable does not automatically add it to the model. I need to add it bu model.add(xv) after 
	// declaring it.
	{
	using namespace std;
	using namespace Eigen;
	ILOSTLBEGIN
		clock_t start;
	start = clock();
	OutData o;	
	int n = x.size();
	int k = 0;
	int l = 0;
	int k_l = 0;
	bool b1=false;
	double fx,ModelReduction,fxs;

	//X xStar,s;
	VectorXd xStar(n),s(n);
	feval(x,fx,s);
	/*cout<<fx<<endl;
	cout<<s<<endl;
	*/
	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);
	cplex.setParam(IloCplex::Param::Simplex::Tolerances::Optimality, 1e-9);
	cplex.setParam(IloCplex::Param::Simplex::Tolerances::Feasibility,1e-9);
	IloNumVar z(env,-IloInfinity,IloInfinity);
	model.add(IloMinimize(env,z));
	IloNumArray lb(env,n);
	IloNumArray ub(env,n);
	IloNumArray xstar(env,n);
	// set lower and upper bounds
	int i;
	for (i = 0; i < n; ++i)
		{
		lb[i]=x(i)-1;
		ub[i]=x(i)+1;
		}
	IloNumVarArray xv(env,lb,ub);	
	model.add(xv);
	IloExpr expr(env);
	expr +=fx-s.transpose()*x;
	// env.out()<<expr<<endl;
	//env.out()<<expr1<<endl;
	for (i = 0; i < n; ++i)
		expr += s(i) * xv[i];
	// env.out()<<expr<<endl;
	//expr1.setLinearCoef(xv[0],s(0));
	//expr1.setLinearCoef(xv[1],s(1));
	//env.out()<<expr1<<endl;
	//add the first cutting plane
	model.add(expr-z<=0);
	cplex.setOut(env.getNullStream());
	// int status;
	const Param p;
	o.t_CPX=0;
	/*env.out()<<xv<<endl;
	env.out()<<lb<<endl;
	env.out()<<ub<<endl;
	env.out()<<xstar<<endl;*/
	while (k<p.Iter_Limi)
		{
		cplex.resetTime();
		cplex.solve();
		// env.out()<<cplex.solve();
		o.t_CPX+=cplex.getTime();
		// env.out()<<cplex.getStatus();
		if (cplex.getStatus()!=IloAlgorithm::Status::Optimal)
			{
			env.out() << "\n Solution status = " << cplex.getStatus() << endl;
			// do something
			o.status=1;
			break;
			}
		// env.out()<<cplex.getObjValue();
		ModelReduction=fx-cplex.getObjValue();
		if (ModelReduction<=p.epsilon_tol)
			{
			b1=true;
			o.status=0;
			break;
			}
		/*env.out()<<xv<<endl;
		env.out()<<xstar<<endl;*/
		cplex.getValues(xv,xstar);
		/*cplex.getValue(xv[0],xstar[0]);
		cplex.getValue(xv[1],xstar[1]);*/
		for (i=0;i<n;++i)
			xStar(i)=xstar[i];
		feval(xStar,fxs,s);
		/*fxs=in.obj(xStar);
		s=in.subgradient(xStar);*/
		if (fx-fxs>=p.m*ModelReduction)// serious step
			{
			//x=xStar;// Do I need this? I don't think so. This way x from imput 
			// arguments can be const.
			fx=fxs;
			k_l=k+1;
			l++;
			for(i=0;i<n;i++)
				{
				lb[i]=xStar(i)-1;
				ub[i]=xStar(i)+1;
				}
			xv.setBounds(lb,ub);
			}

		IloExpr expr1(env);// I cannot just use the old expr because I don't know 
		// know how to reinitiate it.
		expr1+=fxs-s.transpose()*xStar;
		for (i = 0; i < n; ++i)
			expr1 += s(i) * xv[i];
		model.add(expr1<=z);		
		k++;
		}
	env.end();
	if (fx<fxs)
		{
		o.f_final=fx;
		o.Error=fx-prob.f_optimal;
		}
	else
		{
		o.f_final=fxs;
		o.Error=fxs-prob.f_optimal;
		}
	if(b1)
		cout<<"\n SuccesS! Result:\nk="<<k<<", f_val="<<fx<<", error="<<o.Error<<endl;
	else
		{
		cout<<"\n iter limit Result:\nk="<<k<<", f_val="<<fx<<", error="<<o.Error<<endl;
		o.status=2;
		}
	//o.status=o.status;
	o.No_func_eval=k;
	o.k=k;
	o.L=k-l;
	o.time = ( clock() - start ) / (double) CLOCKS_PER_SEC;
	o.t_CPX=o.t_CPX/(double) CLOCKS_PER_SEC;
	cout<<o.status<<" is the status\n";
	cout<<o.Error<<" is the error\n";
	cout<<o.f_final<<" is the final value\n";
	cout<<o.k<<" is the k\n";
	cout<<o.L<<" is the number of null steps\n";
	cout<<o.No_func_eval<<" is the no of func eval\n";
	cout<<o.time<<" is the time\n";
	cout<<o.t_CPX<<" is the cplex time\n";
	return o;
	}