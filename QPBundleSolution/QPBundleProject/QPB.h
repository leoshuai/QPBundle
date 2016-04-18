#pragma once
//typedef Matrix<MSKrealt, Dynamic, 1> VectorXMd;

template<typename T> OutData QPB(ProbData prob, cVec& x, T &feval)
{
	using namespace std;
	using namespace Eigen;
	clock_t start = clock();
	OutData o;
	int n = x.size();
	int k = 0;
	int l = 0;
	int k_l = 0;
	bool b1 = false;
	double fx, ModelReduction, fxs, MSK_time;
	VectorXd xStar(n), sStar(n), s(n), One;
	const VectorXd rs;
	rs.setLinSpaced(n, 0, n - 1);
	One.setOnes(n);
	feval(x, fx, s);
	MSKenv_t      env = NULL;
	MSKtask_t     task = NULL;
	MSKrescodee   r;
	MSKint32t Ind_constraint = 0;
	const Param p;
	o.t_CPX = 0;
	if (MSK_makeenv(&env, NULL) == 0)
	{
		r = MSK_maketask(env, 1, n + 1, &task);
		if (r == MSK_RES_OK)
		{
			if (r == MSK_RES_OK)
				r = MSK_appendcons(task, 1);
			if (r == MSK_RES_OK)
				r = MSK_appendvars(task, n + 1);
			MSKboundkeye  *bkx = (MSKboundkeye  *)MSK_calloctask(task, n + 1, sizeof(MSKboundkeye));
			double *blx = (double *)MSK_calloctask(task, n + 1, sizeof(double));
			double *bux = (double *)MSK_calloctask(task, n + 1, sizeof(double));

			for (size_t i = 0; i < n + 1; i++)
			{
				bkx[i] = MSK_BK_FR;
				blx[i] = -MSK_INFINITY;
				bux[i] = +MSK_INFINITY;
			}
			for (size_t j = 0; j < n + 1 && r == MSK_RES_OK; ++j)
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
			if (r == MSK_RES_OK)
				r = MSK_putcfix(task, 0.5*x.squaredNorm());
			for (size_t j = 0; j < n && r == MSK_RES_OK; ++j)
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
					Ind_constraint, /* Index of the constraint in which the change should occur.*/
					n, /* Index of the variable in which the change should occur.*/
					-1.0);/* New coecient for aij.*/
						  /* Set the bounds on constraints.*/
			if (r == MSK_RES_OK)
				r = MSK_putconbound(task,
					Ind_constraint,/* Index of constraint.*/
					MSK_BK_UP,/* Bound key.*/
					-MSK_INFINITY,/* Numerical value of lower bound.*/
					s.dot(y) - fx);/* Numerical value of upper bound.*/
			if (r == MSK_RES_OK)/* Input the Q for the objective. */
				r = MSK_putqobj(task, n, rs.data(), rs.data(), One.data());
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
					if (r == MSK_RES_OK)
						r = MSK_getdouinf(task, MSK_DINF_OPTIMIZER_TIME, &MSK_time);
					o.t_CPX = o.t_CPX + MSK_time;

					/* Print a summary containing information
					about the solution for debugging purposes*/
					MSK_solutionsummary(task, MSK_STREAM_LOG);

					if (r == MSK_RES_OK)
					{


						MSK_STREAM_LOG MSK_getsolsta(task, MSK_SOL_ITR, &solsta);

						switch (solsta)
						{
						case MSK_SOL_STA_OPTIMAL:
						case MSK_SOL_STA_NEAR_OPTIMAL:
							MSK_getxx(task,
								MSK_SOL_ITR,    /* Request the interior solution. */
								xx);

							printf("Optimal primal solution***********************\n");
							/*for (size_t j = 0; j < n + 1; ++j)
								printf("x[%d]: %e\n", j, xx[j]);*/

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
				ModelReduction = fx - xx[n];
				if (ModelReduction <= p.epsilon_tol)
				{
					b1 = true;
					o.status = 0;
					//break;
				}
				for (auto i = 0; i < n; ++i)
					xStar(i) = xx[i];
				feval(xStar, fxs, sStar);
				if (fx - fxs >= p.m*ModelReduction)
				{
					fx = fxs;
					k_l = k + 1;
					l++;
					//In serious steps I need to update c_j and cfix in the objective
					for (size_t j = 0; j < n && r == MSK_RES_OK; ++j)
					{
						if (r == MSK_RES_OK)
							r = MSK_putcj(task, j, -xStar(j));
					}
					if (r == MSK_RES_OK)
						r = MSK_putcfix(task, 0.5*xStar.squaredNorm());
				}
				if (r == MSK_RES_OK)
					r = MSK_appendcons(task, 1);

				/*Set a row of matrix A*/
				if (r == MSK_RES_OK)
					r = MSK_putarow(task,
						Ind_constraint + 1,                 /* Row index.*/
						n, /* Number of non-zeros in row i.*/
						rs.data(),     /* Pointer to column indexes of row i.*/
						sStar.data());    /* Pointer to values of row i.*/
				if (r == MSK_RES_OK)
					r = MSK_putaij(task,
						Ind_constraint + 1,                 /* Index of the constraint in which the change should occur.*/
						n, /* Index of the variable in which the change should occur.*/
						-1.0);     /* New coecient for aij.*/

								   /* Set the bounds on constraints.*/
				if (r == MSK_RES_OK)
					r = MSK_putconbound(task,
						Ind_constraint + 1,           /* Index of constraint.*/
						MSK_BK_UP,      /* Bound key.*/
						-MSK_INFINITY,      /* Numerical value of lower bound.*/
						sStar.dot(xStar) - fxs);
				k++;

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
		//o.t_CPX = o.t_CPX / (double)CLOCKS_PER_SEC;
		cout << o.status << " is the status\n";
		cout << o.Error << " is the error\n";
		cout << o.f_final << " is the final value\n";
		cout << o.k << " is the k\n";
		cout << o.L << " is the number of null steps\n";
		cout << o.No_func_eval << " is the no of func eval\n";
		cout << o.time << " is the time\n";
		cout << o.t_CPX << " is the MOSEK time\n";
	}
}