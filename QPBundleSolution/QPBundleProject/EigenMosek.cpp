// EigenMosek.cpp : Defines the entry point for the console application.
//Shuai Liu at 24/12/2015 1:57 PM  MSK time is unknown

#include "stdafx.h"
//int testf(static const int &n)
//{
//	//double arr[n];// This does not work. I should use vector.
//	for (auto i = 0; i < n; ++i)
//	{
//		arr[i] = i*1.0;
//		//printf("arr[%d]=%lf\n", i, arr[i]);
//	}
//	return 0;
//}

int main()
{
	using namespace std;
	using namespace Eigen;
	clock_t start=clock();
	VectorXd x(2);
	//x << 2.0, 2.0;
	x << 1.0, 1.0;
	//ProbData prob = { "CB2",1.9522245 };
	ProbData prob = { "DEM",-3.0 };

	OutData o;
	const int n = x.size();
	double *xx = (double*)calloc(n + 1, sizeof(double));
	//double *xx = (double*)calloc(n + 1, sizeof(double));
	/*double  *objval;
	objval = &xx[0];*/
	int k = 0;
	int l = 0;
	int k_l = 0;
	int i,j;
	bool b1 = false;
	double fx, ModelReduction, fxs;

	//X xStar,s;
	VectorXd xStar(n), s(n), sStar(n), y(n), One;
	VectorXi rs(n);
	rs.setLinSpaced(n, 0, n - 1);
	One.setOnes(n);
	DEM(x, fx, s);
	y = x;
	

	/*double te = s.dot(y) - fx;
	cout << s << endl;
	cout << y << endl;
	cout << s.dot(y) << endl;
	cout << fx<<endl;
	cout << te << endl;*/

	MSKenv_t      env = NULL;
	MSKtask_t     task = NULL;
	MSKrescodee   r;
	MSKint32t Ind_constraint = 0;
	const Param p;
	o.t_CPX = 0;
	if (MSK_makeenv(&env, NULL) == 0)
	{
		/* Create the optimization task. */
		r = MSK_maketask(env, 1, n+1, &task);

		if (r == MSK_RES_OK)
		{
			//r = MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);

			/* Append 'NUMCON' empty constraints.
			The constraints will initially have no bounds. */
			if (r == MSK_RES_OK)
				r = MSK_appendcons(task, 1);
			/* Append 'NUMVAR' variables.
			The variables will initially be fixed at zero (x=0). */
			if (r == MSK_RES_OK)
				r = MSK_appendvars(task, n + 1);
			//if (r == MSK_RES_OK)
			//	r = MSK_putcj(task,
			//		0, //Index of the variable for which c should be changed.
			//		1.0);
			MSKboundkeye  *bkx = (MSKboundkeye  *)MSK_calloctask(task, n + 1, sizeof(MSKboundkeye));
			double *blx = (double *)MSK_calloctask(task, n + 1, sizeof(double));
			double *bux = (double *)MSK_calloctask(task, n + 1, sizeof(double));

			for (size_t i = 0; i < n+1; i++)
			{
				bkx[i] = MSK_BK_FR;
				blx[i] = -MSK_INFINITY;
				bux[i] = +MSK_INFINITY;
			}
			
			
			//if (r == MSK_RES_OK)
			//	r = MSK_putaij(task,
			//		Ind_constraint + 1,                 /* Index of the constraint in which the change should occur.*/
			//		0, /* Index of the variable in which the change should occur.*/
			//		-1.0);     /* New coecient for aij.*/

			/* Optionally add a constant term to the objective. */
			if (r == MSK_RES_OK)
				r = MSK_putcfix(task, 0.5*x.squaredNorm());

			/* Set the linear term c_j in the objective.*/
			for (j = 0; j < n && r == MSK_RES_OK; ++j)
			{
				if (r == MSK_RES_OK)
					r = MSK_putcj(task, j, -x(j));
			}
			if (r == MSK_RES_OK)
				r = MSK_putcj(task, n, 1.0);

			/*Set a row of matrix A*/
			if (r == MSK_RES_OK)
				r = MSK_putarow(task,
					Ind_constraint,                 /* Row index.*/
					n, /* Number of non-zeros in row i.*/
					rs.data(),     /* Pointer to column indexes of row i.*/
					s.data());    /* Pointer to values of row i.*/
			if (r == MSK_RES_OK)
				r = MSK_putaij(task,
					Ind_constraint,                 /* Index of the constraint in which the change should occur.*/
					n, /* Index of the variable in which the change should occur.*/
					-1.0);     /* New coecient for aij.*/
							  
			/* Set the bounds on constraints.*/
			if (r == MSK_RES_OK)
				r = MSK_putconbound(task,
					Ind_constraint,           /* Index of constraint.*/
					MSK_BK_UP,      /* Bound key.*/
					-MSK_INFINITY,      /* Numerical value of lower bound.*/
					s.dot(y)-fx);     /* Numerical value of upper bound.*/
			

			for (j = 0; j < n + 1 && r == MSK_RES_OK; ++j)
			{
				/* Set the bounds on variable j.
				blx[j] <= x_j <= bux[j] */
				if (r == MSK_RES_OK)
					r = MSK_putvarbound(task,
						j,           /* Index of variable.*/
						bkx[j],      /* Bound key.*/
						blx[j],      /* Numerical value of lower bound.*/
						bux[j]);     /* Numerical value of upper bound.*/
			}

			
			/*MSKint32t* qsubi = (MSKint32t*)calloc(n, sizeof(MSKint32t));
			MSKint32t* qsubj = (MSKint32t*)calloc(n, sizeof(MSKint32t));
			double *qval = (double*)calloc(n, sizeof(double));*/
			/*MSKint32t    qsubi[n],
				qsubj[n];
			double       qval[n];*/
			//vector<double> qvalt(n);
			/*for (i = 1; i < n+1; ++i)
			{
				qsubi[i-1] = i;
				qsubj[i-1] = i;
				qval[i-1] = 1.0;
			}*/
//The following line sets quadratic terms in a single constraint. I commented it because I am using quadratic objective 
			//function only.

			//if (r == MSK_RES_OK)
			//	r = MSK_putqconk(task,
			//		Ind_constraint,
			//		n,//Number of quadratic terms in the sense of one dimension 
			//		qsubi,
			//		qsubj,
			//		qval);
			
			if (r == MSK_RES_OK)
			{
				/* Input the Q for the objective. */

				r = MSK_putqobj(task, n, rs.data(), rs.data(), One.data());
			}

			if (r == MSK_RES_OK)
				r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
			MSKrescodee trmcode;
			MSKsolstae solsta;
			
			char symname[MSK_MAX_STR_LEN];
			char desc[MSK_MAX_STR_LEN];
			while (k < p.Iter_Limi)
			{
				if (r == MSK_RES_OK)
				{


					/* Run optimizer */
					r = MSK_optimizetrm(task, &trmcode);

					/* Print a summary containing information
					about the solution for debugging purposes*/
					MSK_solutionsummary(task, MSK_STREAM_LOG);

					if (r == MSK_RES_OK)
					{


						MSK_getsolsta(task, MSK_SOL_ITR, &solsta);

						switch (solsta)
						{
						case MSK_SOL_STA_OPTIMAL:
						case MSK_SOL_STA_NEAR_OPTIMAL:
							MSK_getxx(task,
								MSK_SOL_ITR,    /* Request the interior solution. */
								xx);

							printf("Optimal primal solution***********************\n");
							for (j = 0; j < n+1; ++j)
								printf("x[%d]: %e\n", j, xx[j]);

							break;
						case MSK_SOL_STA_DUAL_INFEAS_CER:
						case MSK_SOL_STA_PRIM_INFEAS_CER:
						case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
						case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
							printf("Primal or dual infeasibility certificate found.\n");
							break;

						case MSK_SOL_STA_UNKNOWN:
							printf("The status of the solution could not be determined.\n");
							break;
						default:
							printf("Other solution status.");
							break;
						}
					}
					else
					{
						printf("Error while optimizing.\n");
					}
				}

				if (r != MSK_RES_OK)
				{
					/* In case of an error print error code and description. */


					printf("An error occurred while optimizing.\n");
					MSK_getcodedesc(r,
						symname,
						desc);
					printf("Error %s - '%s'\n", symname, desc);
				}
				/*if (r == MSK_RES_OK)
					r = MSK_getprimalobj(task, MSK_SOL_BAS, objval);*/
				ModelReduction = fx - xx[n];
				if (ModelReduction <= p.epsilon_tol)
				{
					b1 = true;
					o.status = 0;
					//break;
				}
				//MSK_getxx(task, MSK_SOL_BAS, xx);//request the basic solution
				for (auto i = 0; i < n; ++i)
					xStar(i) = xx[i];
				DEM(xStar, fxs, sStar);// Should s be sStar? I think so 

				if (fx - fxs >= p.m*ModelReduction)
				{
					fx = fxs;
					k_l = k + 1;
					l++;
					//In serious steps I need to update c_j and cfix in the objective
					for (j = 0; j < n && r == MSK_RES_OK; ++j)
					{
						if (r == MSK_RES_OK)
							r = MSK_putcj(task, j, -xStar(j));
					}
					if (r == MSK_RES_OK)
						r = MSK_putcfix(task, 0.5*xStar.squaredNorm());
					//for (auto j = 0; j < n && r == MSK_RES_OK; ++j)
					//{
					//	if (r == MSK_RES_OK)
					//		r = MSK_putaij(task,
					//			Ind_constraint,                 /* Index of the constraint in which the change should occur.*/
					//			j+1, /* Index of the variable in which the change should occur.*/
					//			-2.0*xStar(j));
					//}
					//if (r == MSK_RES_OK)
					//	r = MSK_putconbound(task,
					//		Ind_constraint,           /* Index of constraint.*/
					//		MSK_BK_UP,      /* Bound key.*/
					//		-MSK_INFINITY,      /* Numerical value of lower bound.*/
					//		1.0 - xStar.dot(xStar));
				}
				/*if (r == MSK_RES_OK)
					r = MSK_getnumcon(task, &Ind_constraint);*/
				if (r == MSK_RES_OK)
					r = MSK_appendcons(task, 1);
				
				/*Set a row of matrix A*/
				if (r == MSK_RES_OK)
					r = MSK_putarow(task,
						Ind_constraint+1,                 /* Row index.*/
						n, /* Number of non-zeros in row i.*/
						rs.data(),     /* Pointer to column indexes of row i.*/
						sStar.data());    /* Pointer to values of row i.*/
				if (r == MSK_RES_OK)
					r = MSK_putaij(task,
						Ind_constraint+1,                 /* Index of the constraint in which the change should occur.*/
						n, /* Index of the variable in which the change should occur.*/
						-1.0);     /* New coecient for aij.*/

				/* Set the bounds on constraints.*/
				if (r == MSK_RES_OK)
					r = MSK_putconbound(task,
						Ind_constraint+1,           /* Index of constraint.*/
						MSK_BK_UP,      /* Bound key.*/
						-MSK_INFINITY,      /* Numerical value of lower bound.*/
						sStar.dot(xStar) - fxs);     /* Numerical value of upper bound.*/

				//if (r == MSK_RES_OK)
				//	r = MSK_putaij(task,
				//		Ind_constraint,                 /* Index of the constraint in which the change should occur.*/
				//		0, /* Index of the variable in which the change should occur.*/
				//		-1.0);
				//for (auto j = 0; j < n && r == MSK_RES_OK; ++j)
				//{
				//	if (r == MSK_RES_OK)
				//		r = MSK_putaij(task,
				//			Ind_constraint,                 /* Index of the constraint in which the change should occur.*/
				//			j+1, /* Index of the variable in which the change should occur.*/
				//			s(j));     /* New coecient for aij.*/

				//}

				//if (r == MSK_RES_OK)
				//	r = MSK_putconbound(task,
				//		Ind_constraint,           /* Index of constraint.*/
				//		MSK_BK_UP,      /* Bound key.*/
				//		-MSK_INFINITY,      /* Numerical value of lower bound.*/
				//		s.dot(xStar) - fx);

				k++;
				printf("\n======Now we are entering next iteration!=====\n");
			}
			MSK_deletetask(&task);

		}
		MSK_deleteenv(&env);
		if (fx < fxs)
		{
			o.f_final = fx;
			o.Error = fx - prob.f_optimal;
		}
		else
		{
			o.f_final = fxs;
			o.Error = fxs - prob.f_optimal;
		}
		if (b1)
			std::cout << "\n SuccesS! Result:\nk=" << k << ", f_val=" << fx << ", error=" << o.Error << endl;
		else
		{
			cout << "\n iter limit Result:\nk=" << k << ", f_val=" << fx << ", error=" << o.Error << endl;
			o.status = 2;
		}
		//o.status=o.status;
		o.No_func_eval = k;
		o.k = k;
		o.L = k - l;
		o.time = (clock() - start) / (double)CLOCKS_PER_SEC;
		o.t_CPX = o.t_CPX / (double)CLOCKS_PER_SEC;
		cout << o.status << " is the status\n";
		cout << o.Error << " is the error\n";
		cout << o.f_final << " is the final value\n";
		cout << o.k << " is the k\n";
		cout << o.L << " is the number of null steps\n";
		cout << o.No_func_eval << " is the no of func eval\n";
		cout << o.time << " is the time\n";
		cout << o.t_CPX << " is the cplex time\n";
		//return (r);
	}
	int ci;
	cin >> ci;
	return 0;
}